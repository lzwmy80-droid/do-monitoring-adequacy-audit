from __future__ import annotations

import math
import unittest

from g1_stability_gate import (
    audit_exact_dp,
    audit_geometry,
    audit_regular_exact_dp,
    boundary_order_witness,
    deficit_diameter_upper_bound,
    extended_regular_sharpness_witness,
    geometry_bound,
    no_common_anchor_counterexample,
    nonconnected_range_counterexample,
    regular_grid_deficit_diameter_upper_bound,
    twice_distance_integral,
)


class G1StabilityGateTests(unittest.TestCase):
    def test_no_contamination_reduces_to_adaptive_mesh_bound(self) -> None:
        result = deficit_diameter_upper_bound(
            [0.0, 0.5, 2.0], max_slope=3.0, max_replacements=0
        )
        expected = 3.0 * 0.5 * (0.5**2 + 1.5**2)
        self.assertAlmostEqual(float(result["adaptive_bound"]), expected)

    def test_measurement_error_coefficient_uses_full_interval_width(self) -> None:
        result = deficit_diameter_upper_bound(
            [0.0, 2.0],
            max_slope=0.0,
            max_replacements=0,
            half_error=0.4,
        )
        self.assertAlmostEqual(float(result["adaptive_bound"]), 1.6)

    def test_n_must_exceed_twice_the_replacement_budget(self) -> None:
        with self.assertRaisesRegex(ValueError, "n > 2k"):
            deficit_diameter_upper_bound(
                [0.0, 1.0, 2.0, 3.0], max_slope=1.0, max_replacements=2
            )

    def test_geometry_bound_on_boundary_and_interior_missing_blocks(self) -> None:
        times = [0.0, 1.0, 2.0, 3.0, 4.0]
        for common in ((2, 3, 4), (0, 3, 4), (0, 2, 4), (2,)):
            self.assertLessEqual(
                twice_distance_integral(times, common),
                geometry_bound(times, common) + 1e-12,
            )

    def test_exhaustive_random_geometry_search(self) -> None:
        result = audit_geometry(number_of_grids=20, seed=11)
        self.assertLessEqual(float(result["maximum_actual_minus_bound"]), 1e-10)

    def test_random_exact_dp_search(self) -> None:
        result = audit_exact_dp(number_of_cases=500, seed=17)
        self.assertLessEqual(float(result["maximum_width_minus_bound"]), 1e-8)

    def test_random_regular_grid_search_uses_sharper_coefficient(self) -> None:
        result = audit_regular_exact_dp(number_of_cases=500, seed=19)
        self.assertLessEqual(float(result["maximum_width_minus_bound"]), 1e-8)

    def test_regular_grid_bound_matches_general_uniform_bound(self) -> None:
        times = [0.0, 1.0, 2.0, 3.0, 4.0]
        sharp = regular_grid_deficit_diameter_upper_bound(
            times, max_slope=2.0, max_replacements=2
        )
        general = float(
            deficit_diameter_upper_bound(
                times, max_slope=2.0, max_replacements=2
            )["adaptive_bound"]
        )
        self.assertAlmostEqual(sharp, general)

    def test_boundary_witness_has_exact_quadratic_order(self) -> None:
        rows = boundary_order_witness(max_k=5)
        for row in rows:
            self.assertTrue(
                math.isclose(
                    float(row["robust_deficit_width"]),
                    float(row["four_L_k2_h2"]),
                    rel_tol=1e-11,
                    abs_tol=1e-11,
                )
            )

    def test_extended_witness_attains_full_regular_bound(self) -> None:
        rows = extended_regular_sharpness_witness(
            max_k=3, extra_interval_counts=(0, 1, 3, 8)
        )
        self.assertEqual(len(rows), 12)
        for row in rows:
            self.assertTrue(
                math.isclose(
                    float(row["robust_deficit_width"]),
                    float(row["regular_sharp_bound"]),
                    rel_tol=1e-11,
                    abs_tol=1e-11,
                )
            )

    def test_no_common_anchor_counterexample_breaks_unconditional_formula(self) -> None:
        result = no_common_anchor_counterexample(100.0)
        self.assertEqual(float(result["robust_deficit_width"]), 100.0)
        self.assertEqual(float(result["grid_only_candidate_rhs"]), 0.0)

    def test_endpoint_hull_need_not_equal_functional_image(self) -> None:
        result = nonconnected_range_counterexample(100.0)
        self.assertEqual(result["exact_deficit_image"], [0.0, 100.0])
        self.assertFalse(bool(result["midpoint_in_exact_image"]))


if __name__ == "__main__":
    unittest.main()
