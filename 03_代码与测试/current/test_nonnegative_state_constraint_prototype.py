# -*- coding: utf-8 -*-
"""Regression tests for the independent nonnegative-state oracle."""

from __future__ import annotations

import math
import random
import unittest

from nonnegative_state_constraint_prototype import (
    NonnegativePrototypeInfeasible,
    Segment,
    clip_segments_to_floor,
    enumerate_nonnegative_bounds,
    integrate_deficit,
    minimum_replacements_for_nonnegative_lipschitz,
)
from robust_exact_current import (
    identify_hypoxia_burden_with_replacements,
    minimum_replacements_for_lipschitz,
)


class NonnegativeStateConstraintPrototypeTests(unittest.TestCase):
    def test_floor_validity_is_decided_before_value_normalization(self) -> None:
        minimum = minimum_replacements_for_lipschitz(
            [0.0, 1.0, 2.0],
            [1e16, -1.0, 0.0],
            0.0,
            state_lower_bound=0.0,
        )
        self.assertEqual(minimum.minimum_replacements, 2)

    def test_all_deleted_branch_bypasses_ill_conditioned_normalization(
        self,
    ) -> None:
        result = identify_hypoxia_burden_with_replacements(
            [0.0, 2.0],
            [1e16, 0.0],
            threshold=1.0,
            max_slope=0.0,
            max_replacements=2,
            state_lower_bound=0.0,
        )
        self.assertEqual(result.occupation_lower, 0.0)
        self.assertEqual(result.occupation_upper, 2.0)
        self.assertEqual(result.oxygen_deficit_lower, 0.0)
        self.assertEqual(result.oxygen_deficit_upper, 2.0)

    def test_all_deleted_deficit_uses_one_exact_product(self) -> None:
        tiny_horizon = 1e-308
        result = identify_hypoxia_burden_with_replacements(
            [0.0, tiny_horizon],
            [0.0, 0.0],
            threshold=1e308,
            max_slope=0.0,
            max_replacements=2,
            state_lower_bound=-1e308,
        )
        self.assertTrue(
            math.isclose(
                result.oxygen_deficit_upper,
                2.0,
                rel_tol=2e-16,
                abs_tol=0.0,
            )
        )

    def test_unrepresentable_positive_endpoint_raises_instead_of_zeroing(
        self,
    ) -> None:
        minimum_subnormal = math.ulp(0.0)
        with self.assertRaises(FloatingPointError):
            identify_hypoxia_burden_with_replacements(
                [0.0, minimum_subnormal],
                [0.0, 0.0],
                threshold=minimum_subnormal,
                max_slope=0.0,
                max_replacements=2,
                state_lower_bound=0.0,
            )

    def test_only_deficit_upper_changes_when_reports_and_threshold_are_positive(
        self,
    ) -> None:
        unconstrained = identify_hypoxia_burden_with_replacements(
            [0.0, 2.0],
            [1.0, 1.0],
            threshold=1.0,
            max_slope=2.0,
            max_replacements=0,
        )
        constrained = enumerate_nonnegative_bounds(
            [0.0, 2.0],
            [1.0, 1.0],
            threshold=1.0,
            max_slope=2.0,
            max_replacements=0,
        )
        production = identify_hypoxia_burden_with_replacements(
            [0.0, 2.0],
            [1.0, 1.0],
            threshold=1.0,
            max_slope=2.0,
            max_replacements=0,
            state_lower_bound=0.0,
        )
        self.assertAlmostEqual(
            constrained.occupation_lower,
            unconstrained.occupation_lower,
        )
        self.assertAlmostEqual(
            constrained.occupation_upper,
            unconstrained.occupation_upper,
        )
        self.assertAlmostEqual(
            constrained.oxygen_deficit_lower,
            unconstrained.oxygen_deficit_lower,
        )
        self.assertAlmostEqual(unconstrained.oxygen_deficit_upper, 2.0)
        self.assertAlmostEqual(constrained.oxygen_deficit_upper, 1.5)
        self.assertAlmostEqual(
            production.oxygen_deficit_upper,
            constrained.oxygen_deficit_upper,
        )

    def test_no_change_when_the_lower_envelope_is_already_nonnegative(self) -> None:
        unconstrained = identify_hypoxia_burden_with_replacements(
            [0.0, 1.0],
            [2.0, 2.0],
            threshold=5.0,
            max_slope=1.0,
            max_replacements=0,
        )
        constrained = enumerate_nonnegative_bounds(
            [0.0, 1.0],
            [2.0, 2.0],
            threshold=5.0,
            max_slope=1.0,
            max_replacements=0,
        )
        self.assertAlmostEqual(
            constrained.oxygen_deficit_upper,
            unconstrained.oxygen_deficit_upper,
        )

    def test_all_reports_replaceable_has_finite_attained_upper_deficit(self) -> None:
        constrained = enumerate_nonnegative_bounds(
            [0.0, 2.0],
            [1.0, 1.0],
            threshold=5.0,
            max_slope=2.0,
            max_replacements=2,
        )
        production = identify_hypoxia_burden_with_replacements(
            [0.0, 2.0],
            [1.0, 1.0],
            threshold=5.0,
            max_slope=2.0,
            max_replacements=2,
            state_lower_bound=0.0,
        )
        self.assertAlmostEqual(constrained.occupation_upper, 2.0)
        self.assertAlmostEqual(constrained.oxygen_deficit_upper, 10.0)
        self.assertTrue(constrained.oxygen_deficit_upper_witness.attained)
        self.assertEqual(
            constrained.oxygen_deficit_upper_witness.retained_indices,
            (),
        )
        self.assertAlmostEqual(production.oxygen_deficit_upper, 10.0)
        self.assertTrue(production.oxygen_deficit_upper_witness.attained)

    def test_negative_report_is_a_mandatory_replacement(self) -> None:
        self.assertEqual(
            minimum_replacements_for_nonnegative_lipschitz(
                [0.0, 1.0],
                [-1.0, 1.0],
                2.0,
            ),
            1,
        )
        production_minimum = minimum_replacements_for_lipschitz(
            [0.0, 1.0],
            [-1.0, 1.0],
            2.0,
            state_lower_bound=0.0,
        )
        self.assertEqual(production_minimum.minimum_replacements, 1)
        with self.assertRaises(NonnegativePrototypeInfeasible):
            enumerate_nonnegative_bounds(
                [0.0, 1.0],
                [-1.0, 1.0],
                threshold=1.0,
                max_slope=2.0,
                max_replacements=0,
            )

    def test_nonpositive_threshold_collapses_both_functionals(self) -> None:
        for threshold in (0.0, -1.0):
            constrained = enumerate_nonnegative_bounds(
                [0.0, 1.0],
                [0.5, 0.5],
                threshold=threshold,
                max_slope=2.0,
                max_replacements=2,
            )
            self.assertEqual(constrained.occupation_lower, 0.0)
            self.assertEqual(constrained.occupation_upper, 0.0)
            self.assertEqual(constrained.oxygen_deficit_lower, 0.0)
            self.assertEqual(constrained.oxygen_deficit_upper, 0.0)

        production = identify_hypoxia_burden_with_replacements(
            [0.0, 1.0],
            [0.5, 0.5],
            threshold=0.0,
            max_slope=2.0,
            max_replacements=0,
            state_lower_bound=0.0,
        )
        self.assertEqual(production.occupation_lower, 0.0)
        self.assertEqual(production.occupation_upper, 0.0)
        self.assertEqual(production.oxygen_deficit_lower, 0.0)
        self.assertEqual(production.oxygen_deficit_upper, 0.0)

    def test_clipping_identity_matches_hinge_difference(self) -> None:
        segments = [
            Segment(0.0, 1.0, -2.0, 1.0),
            Segment(1.0, 2.0, 2.0, -3.0),
        ]
        clipped = clip_segments_to_floor(segments)
        direct = integrate_deficit(clipped, 1.0)
        difference = integrate_deficit(segments, 1.0) - integrate_deficit(
            segments, 0.0
        )
        self.assertAlmostEqual(direct, difference)
        self.assertAlmostEqual(direct, 1.5)

    def test_random_small_cases_match_three_unconstrained_endpoints(self) -> None:
        generator = random.Random(20260727)
        for _ in range(40):
            times = [float(index) for index in range(5)]
            max_slope = 2.0
            base = [1.0]
            for _index in range(1, 5):
                increment = generator.uniform(-1.0, 1.0)
                base.append(max(0.0, base[-1] + increment))
            reports = list(base)
            reports[generator.randrange(5)] += generator.uniform(2.1, 4.0)
            threshold = 1.0
            unconstrained = identify_hypoxia_burden_with_replacements(
                times,
                reports,
                threshold=threshold,
                max_slope=max_slope,
                max_replacements=1,
            )
            production = identify_hypoxia_burden_with_replacements(
                times,
                reports,
                threshold=threshold,
                max_slope=max_slope,
                max_replacements=1,
                state_lower_bound=0.0,
            )
            constrained = enumerate_nonnegative_bounds(
                times,
                reports,
                threshold=threshold,
                max_slope=max_slope,
                max_replacements=1,
            )
            self.assertTrue(
                math.isclose(
                    constrained.occupation_lower,
                    unconstrained.occupation_lower,
                    abs_tol=1e-10,
                )
            )
            self.assertTrue(
                math.isclose(
                    constrained.occupation_upper,
                    unconstrained.occupation_upper,
                    abs_tol=1e-10,
                )
            )
            self.assertTrue(
                math.isclose(
                    constrained.oxygen_deficit_lower,
                    unconstrained.oxygen_deficit_lower,
                    abs_tol=1e-10,
                )
            )
            self.assertLessEqual(
                constrained.oxygen_deficit_upper,
                unconstrained.oxygen_deficit_upper + 1e-10,
            )
            self.assertLessEqual(
                constrained.oxygen_deficit_upper
                - constrained.oxygen_deficit_lower,
                threshold * constrained.horizon + 1e-10,
            )
            oracle_values = (
                constrained.occupation_lower,
                constrained.occupation_upper,
                constrained.oxygen_deficit_lower,
                constrained.oxygen_deficit_upper,
            )
            production_values = (
                production.occupation_lower,
                production.occupation_upper,
                production.oxygen_deficit_lower,
                production.oxygen_deficit_upper,
            )
            for oracle_value, production_value in zip(
                oracle_values, production_values
            ):
                self.assertTrue(
                    math.isclose(
                        oracle_value,
                        production_value,
                        rel_tol=1e-11,
                        abs_tol=1e-11,
                    ),
                    (oracle_values, production_values),
                )


if __name__ == "__main__":
    unittest.main()
