from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from future_partial_identification import (  # noqa: E402
    InfeasibleObservationIntervals,
    construct_envelopes,
    deficit_width_bound_exact_observations,
    deficit_width_bound_uniform_error,
    identify_hypoxia_burden,
    occupation_width_margin_bound,
)


class PartialIdentificationTests(unittest.TestCase):
    def test_short_above_threshold_gap_excludes_hypoxia(self) -> None:
        result = identify_hypoxia_burden(
            [0.0, 1.0],
            [6.0, 6.0],
            threshold=5.0,
            max_slope=1.0,
        )
        self.assertAlmostEqual(result.occupation_lower, 0.0)
        self.assertAlmostEqual(result.occupation_upper, 0.0)
        self.assertAlmostEqual(result.oxygen_deficit_lower, 0.0)
        self.assertAlmostEqual(result.oxygen_deficit_upper, 0.0)

    def test_hidden_hypoxia_has_sharp_nontrivial_upper_bounds(self) -> None:
        result = identify_hypoxia_burden(
            [0.0, 2.0],
            [6.0, 6.0],
            threshold=5.0,
            max_slope=2.0,
        )
        self.assertAlmostEqual(result.occupation_lower, 0.0)
        self.assertAlmostEqual(result.occupation_upper, 1.0)
        self.assertAlmostEqual(result.oxygen_deficit_lower, 0.0)
        self.assertAlmostEqual(result.oxygen_deficit_upper, 0.5)
        self.assertEqual(
            result.event_count_status,
            "not_bounded_in_general_by_lipschitz",
        )

    def test_certain_hypoxia_bounds_deficit_not_just_time(self) -> None:
        result = identify_hypoxia_burden(
            [0.0, 1.0],
            [4.0, 4.0],
            threshold=5.0,
            max_slope=1.0,
        )
        self.assertAlmostEqual(result.occupation_lower, 1.0)
        self.assertAlmostEqual(result.occupation_upper, 1.0)
        self.assertAlmostEqual(result.oxygen_deficit_lower, 0.75)
        self.assertAlmostEqual(result.oxygen_deficit_upper, 1.25)

    def test_measurement_intervals_use_global_feasible_envelopes(self) -> None:
        result = identify_hypoxia_burden(
            [0.0, 2.0],
            observation_lower=[5.5, 5.5],
            observation_upper=[6.5, 6.5],
            threshold=5.0,
            max_slope=1.0,
        )
        self.assertEqual(result.observation_mode, "interval")
        self.assertAlmostEqual(result.occupation_lower, 0.0)
        self.assertAlmostEqual(result.occupation_upper, 1.0)
        self.assertAlmostEqual(result.oxygen_deficit_lower, 0.0)
        self.assertAlmostEqual(result.oxygen_deficit_upper, 0.25)

    def test_measurement_error_must_propagate_across_all_observations(self) -> None:
        result = identify_hypoxia_burden(
            [0.0, 1.0, 2.0, 3.0],
            observation_lower=[0.0, 0.0, 0.0, 0.0],
            observation_upper=[0.0, 10.0, 10.0, 0.0],
            threshold=5.0,
            max_slope=1.0,
        )
        # The exact zero endpoints constrain the two broad interior intervals.
        # Treating each adjacent raw interval independently would miss this.
        self.assertAlmostEqual(result.occupation_lower, 3.0)
        self.assertAlmostEqual(result.occupation_upper, 3.0)

    def test_distant_observation_constrains_interval_envelope(self) -> None:
        lower, upper = construct_envelopes(
            [0.0, 1.0, 2.0],
            [0.0, -100.0, 10.0],
            [0.0, 100.0, 10.0],
            max_slope=5.0,
        )
        lower_at_middle = max(
            segment.value(1.0)
            for segment in lower
            if segment.start <= 1.0 <= segment.end
        )
        upper_at_middle = min(
            segment.value(1.0)
            for segment in upper
            if segment.start <= 1.0 <= segment.end
        )
        self.assertAlmostEqual(lower_at_middle, 5.0)
        self.assertAlmostEqual(upper_at_middle, 5.0)

    def test_infeasible_exact_observations_are_rejected(self) -> None:
        with self.assertRaises(InfeasibleObservationIntervals):
            identify_hypoxia_burden(
                [0.0, 1.0],
                [0.0, 2.0],
                threshold=1.0,
                max_slope=1.0,
            )

    def test_zero_slope_requires_a_common_interval(self) -> None:
        feasible = identify_hypoxia_burden(
            [0.0, 4.0],
            observation_lower=[4.0, 4.5],
            observation_upper=[5.0, 5.5],
            threshold=4.75,
            max_slope=0.0,
        )
        self.assertAlmostEqual(feasible.occupation_lower, 0.0)
        self.assertAlmostEqual(feasible.occupation_upper, 4.0)
        with self.assertRaises(InfeasibleObservationIntervals):
            identify_hypoxia_burden(
                [0.0, 4.0],
                observation_lower=[4.0, 5.0],
                observation_upper=[4.5, 5.5],
                threshold=4.75,
                max_slope=0.0,
            )

    def test_bounds_are_ordered_and_finite(self) -> None:
        result = identify_hypoxia_burden(
            [10.0, 11.5, 13.0],
            [5.5, 4.0, 5.2],
            threshold=5.0,
            max_slope=2.0,
        )
        self.assertLessEqual(result.occupation_lower, result.occupation_upper)
        self.assertLessEqual(result.oxygen_deficit_lower, result.oxygen_deficit_upper)
        self.assertTrue(math.isfinite(result.oxygen_deficit_upper))

    def test_deficit_width_is_bounded_by_irregular_grid_corridor(self) -> None:
        times = [0.0, 1.0, 3.0]
        values = [5.0, 5.5, 4.5]
        max_slope = 1.0
        result = identify_hypoxia_burden(
            times,
            values,
            threshold=5.0,
            max_slope=max_slope,
        )
        refined = deficit_width_bound_exact_observations(
            times,
            values,
            max_slope,
        )
        coarse = deficit_width_bound_uniform_error(times, max_slope, 0.0)
        self.assertLessEqual(
            result.oxygen_deficit_upper - result.oxygen_deficit_lower,
            refined + 1e-12,
        )
        self.assertLessEqual(refined, coarse + 1e-12)

    def test_threshold_samples_show_occupation_has_no_uniform_convergence(self) -> None:
        for number_of_intervals in (2, 10, 100):
            times = [index / number_of_intervals for index in range(number_of_intervals + 1)]
            result = identify_hypoxia_burden(
                times,
                [5.0] * len(times),
                threshold=5.0,
                max_slope=1.0,
            )
            self.assertAlmostEqual(result.occupation_lower, 0.0)
            self.assertAlmostEqual(result.occupation_upper, 1.0)

    def test_margin_design_bound_uses_maximum_gap_and_error(self) -> None:
        bound = occupation_width_margin_bound(
            [0.0, 1.0, 3.0],
            max_slope=0.5,
            error_halfwidth=0.1,
            margin_constant=2.0,
            margin_exponent=1.0,
        )
        self.assertAlmostEqual(bound, 2.4)


if __name__ == "__main__":
    unittest.main()
