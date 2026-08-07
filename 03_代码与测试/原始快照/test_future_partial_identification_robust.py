from __future__ import annotations

import itertools
import math
import random
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from future_partial_identification import identify_hypoxia_burden  # noqa: E402
from future_partial_identification_robust import (  # noqa: E402
    InsufficientReplacementBudget,
    _as_inputs,
    _compatible,
    _retained_subset_value,
    identify_hypoxia_burden_with_replacements,
    k_min,
    minimum_replacements_for_lipschitz,
)


def brute_force_endpoints(
    times: list[float],
    values: list[float],
    *,
    threshold: float,
    max_slope: float,
    max_replacements: int,
) -> tuple[float, float, float, float]:
    """Independent subset enumeration; continuous integrals reuse primitives."""

    t, y, tolerance = _as_inputs(times, values, max_slope)
    number = len(times)
    candidates: list[tuple[float, float, float, float]] = []
    for retained_count in range(max(1, number - max_replacements), number + 1):
        for retained in itertools.combinations(range(number), retained_count):
            if any(
                not _compatible(t, y, left, right, max_slope, tolerance)
                for left, right in zip(retained[:-1], retained[1:])
            ):
                continue
            candidates.append(
                (
                    _retained_subset_value(
                        t,
                        y,
                        retained,
                        max_slope,
                        threshold,
                        "upper",
                        "occupation",
                    ),
                    _retained_subset_value(
                        t,
                        y,
                        retained,
                        max_slope,
                        threshold,
                        "lower",
                        "occupation",
                    ),
                    _retained_subset_value(
                        t,
                        y,
                        retained,
                        max_slope,
                        threshold,
                        "upper",
                        "deficit",
                    ),
                    _retained_subset_value(
                        t,
                        y,
                        retained,
                        max_slope,
                        threshold,
                        "lower",
                        "deficit",
                    ),
                )
            )
    if not candidates:
        raise InsufficientReplacementBudget
    return (
        min(item[0] for item in candidates),
        max(item[1] for item in candidates),
        min(item[2] for item in candidates),
        max(item[3] for item in candidates),
    )


class RobustPartialIdentificationTests(unittest.TestCase):
    def assert_witness_valid(
        self,
        times: list[float],
        values: list[float],
        max_slope: float,
        budget: int,
        witness,
    ) -> None:
        self.assertLessEqual(len(witness.deleted_indices), budget)
        self.assertEqual(
            set(witness.retained_indices) | set(witness.deleted_indices),
            set(range(len(times))),
        )
        for left, right in zip(
            witness.retained_indices[:-1], witness.retained_indices[1:]
        ):
            self.assertLessEqual(
                abs(values[right] - values[left]),
                max_slope * (times[right] - times[left]) + 1e-10,
            )

    def test_k_zero_matches_exact_identification(self) -> None:
        times = [0.0, 1.0, 3.0]
        values = [5.5, 4.8, 5.4]
        exact = identify_hypoxia_burden(
            times, values, threshold=5.0, max_slope=1.0
        )
        robust = identify_hypoxia_burden_with_replacements(
            times,
            values,
            threshold=5.0,
            max_slope=1.0,
            max_replacements=0,
        )
        self.assertAlmostEqual(robust.occupation_lower, exact.occupation_lower)
        self.assertAlmostEqual(robust.occupation_upper, exact.occupation_upper)
        self.assertAlmostEqual(
            robust.oxygen_deficit_lower, exact.oxygen_deficit_lower
        )
        self.assertAlmostEqual(
            robust.oxygen_deficit_upper, exact.oxygen_deficit_upper
        )
        self.assertEqual(robust.occupation_lower_witness.retained_indices, (0, 1, 2))

    def test_endpoint_may_be_deleted(self) -> None:
        times = [0.0, 1.0, 2.0]
        values = [100.0, 5.0, 5.0]
        result = identify_hypoxia_burden_with_replacements(
            times,
            values,
            threshold=5.0,
            max_slope=1.0,
            max_replacements=1,
        )
        for witness in (
            result.occupation_lower_witness,
            result.occupation_upper_witness,
            result.oxygen_deficit_lower_witness,
            result.oxygen_deficit_upper_witness,
        ):
            self.assertEqual(witness.deleted_indices, (0,))
            self.assertEqual(witness.retained_indices, (1, 2))

    def test_infeasible_budget_reports_k_min_and_one_replacement_recovers(self) -> None:
        times = [0.0, 1.0, 2.0]
        values = [0.0, 100.0, 0.0]
        witness = minimum_replacements_for_lipschitz(times, values, 1.0)
        self.assertEqual(witness.minimum_replacements, 1)
        self.assertEqual(witness.retained_indices, (0, 2))
        self.assertEqual(k_min(times, values, 1.0), 1)
        with self.assertRaisesRegex(InsufficientReplacementBudget, "k_min"):
            identify_hypoxia_burden_with_replacements(
                times,
                values,
                threshold=5.0,
                max_slope=1.0,
                max_replacements=0,
            )
        recovered = identify_hypoxia_burden_with_replacements(
            times,
            values,
            threshold=5.0,
            max_slope=1.0,
            max_replacements=1,
        )
        self.assertEqual(recovered.minimum_replacements_for_feasibility, 1)

    def test_all_values_replaceable_has_correct_extended_real_bounds(self) -> None:
        result = identify_hypoxia_burden_with_replacements(
            [10.0, 12.0],
            [0.0, 100.0],
            threshold=5.0,
            max_slope=1.0,
            max_replacements=2,
        )
        self.assertEqual(result.occupation_lower, 0.0)
        self.assertEqual(result.occupation_upper, 2.0)
        self.assertEqual(result.oxygen_deficit_lower, 0.0)
        self.assertTrue(math.isinf(result.oxygen_deficit_upper))
        self.assertFalse(result.oxygen_deficit_upper_witness.attained)
        self.assertEqual(result.oxygen_deficit_upper_witness.retained_indices, ())

    def test_threshold_boundary_preserves_full_occupation_ambiguity(self) -> None:
        result = identify_hypoxia_burden_with_replacements(
            [0.0, 0.5, 1.0],
            [5.0, 5.0, 5.0],
            threshold=5.0,
            max_slope=1.0,
            max_replacements=0,
        )
        self.assertAlmostEqual(result.occupation_lower, 0.0)
        self.assertAlmostEqual(result.occupation_upper, 1.0)

    def test_dp_matches_brute_force_on_random_small_instances(self) -> None:
        generator = random.Random(20260716)
        for number in range(2, 8):
            for _ in range(30):
                increments = [generator.uniform(0.2, 1.5) for _ in range(number - 1)]
                times = [0.0]
                for increment in increments:
                    times.append(times[-1] + increment)
                values = [generator.uniform(2.0, 8.0) for _ in range(number)]
                max_slope = generator.uniform(0.0, 4.0)
                threshold = generator.uniform(3.0, 7.0)
                minimum = k_min(times, values, max_slope)
                budget = generator.randint(minimum, number - 1)
                expected = brute_force_endpoints(
                    times,
                    values,
                    threshold=threshold,
                    max_slope=max_slope,
                    max_replacements=budget,
                )
                actual = identify_hypoxia_burden_with_replacements(
                    times,
                    values,
                    threshold=threshold,
                    max_slope=max_slope,
                    max_replacements=budget,
                )
                observed = (
                    actual.occupation_lower,
                    actual.occupation_upper,
                    actual.oxygen_deficit_lower,
                    actual.oxygen_deficit_upper,
                )
                np.testing.assert_allclose(observed, expected, rtol=1e-10, atol=1e-10)
                for witness in (
                    actual.occupation_lower_witness,
                    actual.occupation_upper_witness,
                    actual.oxygen_deficit_lower_witness,
                    actual.oxygen_deficit_upper_witness,
                ):
                    self.assert_witness_valid(
                        times, values, max_slope, budget, witness
                    )
                    self.assertTrue(witness.attained)
                normalized_times, normalized_values, _ = _as_inputs(
                    times, values, max_slope
                )
                for witness, envelope, functional in (
                    (actual.occupation_lower_witness, "upper", "occupation"),
                    (actual.occupation_upper_witness, "lower", "occupation"),
                    (actual.oxygen_deficit_lower_witness, "upper", "deficit"),
                    (actual.oxygen_deficit_upper_witness, "lower", "deficit"),
                ):
                    recomputed = _retained_subset_value(
                        normalized_times,
                        normalized_values,
                        witness.retained_indices,
                        max_slope,
                        threshold,
                        envelope,
                        functional,
                    )
                    self.assertAlmostEqual(recomputed, witness.value)


if __name__ == "__main__":
    unittest.main()
