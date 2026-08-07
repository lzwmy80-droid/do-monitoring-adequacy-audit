from __future__ import annotations

import math
import sys
import unittest

from independent_endpoint_oracle import (
    OracleInfeasible,
    audit_discrete_boundary_cases,
    audit_random_cases,
    audit_vertical_translation_invariance,
    enumerate_endpoint_oracle,
)
from robust_exact_current import (
    identify_hypoxia_burden_with_replacements,
    k_min,
)


class IndependentEndpointOracleTests(unittest.TestCase):
    def test_exact_binary64_compatibility_rejects_one_ulp_violation(self) -> None:
        one_ulp_above = math.nextafter(1.0, math.inf)
        self.assertEqual(
            k_min((0.0, 1.0), (0.0, one_ulp_above), 1.0),
            1,
        )

    def test_exact_compatibility_avoids_subtraction_rounding_false_positive(
        self,
    ) -> None:
        self.assertEqual(
            k_min((0.0, 1.0), (1.0, -math.ldexp(1.0, -54)), 1.0),
            1,
        )

    def test_exact_compatibility_avoids_infinity_comparison_false_positive(
        self,
    ) -> None:
        maximum = sys.float_info.max
        almost_two = math.nextafter(2.0, 0.0)
        self.assertEqual(
            k_min((0.0, almost_two), (-maximum, maximum), maximum),
            1,
        )

    def test_k_min_does_not_alias_distinct_reports_during_normalization(
        self,
    ) -> None:
        reports = (1e12, 0.0, 1.0, math.nextafter(1.0, math.inf))
        self.assertEqual(k_min((0.0, 1.0, 2.0, 3.0), reports, 0.0), 3)
        with self.assertRaises(FloatingPointError):
            identify_hypoxia_burden_with_replacements(
                (0.0, 1.0, 2.0, 3.0),
                reports,
                threshold=0.5,
                max_slope=0.0,
                max_replacements=3,
            )

    def test_tiny_positive_slope_is_rejected_if_normalization_erases_it(
        self,
    ) -> None:
        tiny_slope = math.ulp(0.0)
        self.assertEqual(
            k_min((0.0, 1.0), (0.0, 1e308), tiny_slope),
            1,
        )
        with self.assertRaises(FloatingPointError):
            identify_hypoxia_burden_with_replacements(
                (0.0, 1.0),
                (0.0, 1e308),
                threshold=0.0,
                max_slope=tiny_slope,
                max_replacements=1,
            )

    def test_near_vertex_threshold_uses_exact_local_fallback(
        self,
    ) -> None:
        result = identify_hypoxia_burden_with_replacements(
            (0.0, 1e50),
            (0.0006197765650313158, -23.482067348105762),
            threshold=-24.18785045032845,
            max_slope=2.4894253329116173e-49,
            max_replacements=0,
        )
        self.assertEqual(result.occupation_upper, 1.493678147677149e34)
        self.assertTrue(
            math.isclose(
                result.oxygen_deficit_upper,
                1.3885232882494042e19,
                rel_tol=2e-16,
                abs_tol=0.0,
            )
        )

    def test_critical_fallback_threshold_preserves_small_endpoint_accuracy(
        self,
    ) -> None:
        result = identify_hypoxia_burden_with_replacements(
            (0.0, 1.0),
            (1685.3354141780637, 1685.2457257719755),
            threshold=1685.1749213108078,
            max_slope=0.23129732842418432,
            max_replacements=0,
        )
        self.assertEqual(
            result.occupation_upper,
            2.614790221681251e-12,
        )

    def test_exact_tangency_remains_zero_under_extreme_scaling(self) -> None:
        result = identify_hypoxia_burden_with_replacements(
            (0.0, math.ldexp(1.0, 500)),
            (3.11520354590533e-20, 6.2778872589712035e-21),
            threshold=-2.362608775210443e-21,
            max_slope=math.ldexp(4.215514026844539e-20, -500),
            max_replacements=0,
        )
        self.assertEqual(result.occupation_upper, 0.0)
        self.assertEqual(result.oxygen_deficit_upper, 0.0)

    def test_anchor_threshold_tangency_uses_exact_boundary_fallback(
        self,
    ) -> None:
        horizon = math.ldexp(1.0, 500)
        threshold = math.nextafter(1.0, 0.0)
        result = identify_hypoxia_burden_with_replacements(
            (0.0, 0.75 * horizon, horizon),
            (-2.0, threshold, -2.0),
            threshold=threshold,
            max_slope=math.ldexp(1.0, -501),
            max_replacements=2,
        )
        self.assertEqual(result.occupation_lower, 0.0)
        self.assertEqual(
            result.occupation_lower_witness.retained_indices,
            (1,),
        )

    def test_derived_cancellation_uses_exact_local_fallback(self) -> None:
        unit = math.nextafter(1.0, math.inf)
        time_scale = math.ldexp(1.0, 104)
        value_scale = math.ldexp(1.0, 106)
        horizon = time_scale * unit
        slope = (value_scale / time_scale) * unit
        threshold = -value_scale * (0.5 + (unit - 1.0))
        result = identify_hypoxia_burden_with_replacements(
            (0.0, horizon),
            (0.0, 0.0),
            threshold=threshold,
            max_slope=slope,
            max_replacements=0,
        )
        expected = 1.0 / unit
        self.assertEqual(result.occupation_upper, expected)
        self.assertEqual(result.oxygen_deficit_upper, expected)

    def test_tiny_constant_conflict_is_not_hidden_by_absolute_tolerance(self) -> None:
        self.assertEqual(k_min((0.0, 1.0), (0.0, 1e-14), 0.0), 1)
        with self.assertRaises(OracleInfeasible):
            enumerate_endpoint_oracle(
                (0.0, 1.0),
                (0.0, 1e-14),
                threshold=0.0,
                max_slope=0.0,
                max_replacements=0,
            )

    def test_small_slope_large_horizon_analytic_tent(self) -> None:
        result = identify_hypoxia_burden_with_replacements(
            (0.0, 1e16),
            (0.0, 0.0),
            threshold=-0.25,
            max_slope=1e-16,
            max_replacements=0,
        )
        self.assertEqual(result.occupation_upper, 5e15)
        self.assertEqual(result.oxygen_deficit_upper, 6.25e14)

        oracle = enumerate_endpoint_oracle(
            (0.0, 1e16),
            (0.0, 0.0),
            threshold=-0.25,
            max_slope=1e-16,
            max_replacements=0,
        )
        self.assertEqual(result.occupation_upper, oracle.occupation_upper)
        self.assertEqual(result.oxygen_deficit_upper, oracle.deficit_upper)

    def test_time_unit_scaling_transforms_endpoint_units(self) -> None:
        base = identify_hypoxia_burden_with_replacements(
            (0.0, 1.0, 2.0),
            (0.0, 0.5, 1.0),
            threshold=0.75,
            max_slope=1.0,
            max_replacements=0,
        )
        factor = 3600.0
        scaled = identify_hypoxia_burden_with_replacements(
            (0.0, factor, 2.0 * factor),
            (0.0, 0.5, 1.0),
            threshold=0.75,
            max_slope=1.0 / factor,
            max_replacements=0,
        )
        base_values = (
            base.occupation_lower,
            base.occupation_upper,
            base.oxygen_deficit_lower,
            base.oxygen_deficit_upper,
        )
        scaled_values = (
            scaled.occupation_lower,
            scaled.occupation_upper,
            scaled.oxygen_deficit_lower,
            scaled.oxygen_deficit_upper,
        )
        for original, transformed in zip(base_values, scaled_values):
            self.assertAlmostEqual(transformed, factor * original)

    def test_value_unit_scaling_preserves_occupation_and_scales_deficit(self) -> None:
        base = identify_hypoxia_burden_with_replacements(
            (0.0, 1.0, 2.0),
            (0.0, 0.5, 1.0),
            threshold=0.75,
            max_slope=1.0,
            max_replacements=0,
        )
        factor = 1000.0
        scaled = identify_hypoxia_burden_with_replacements(
            (0.0, 1.0, 2.0),
            (0.0, 500.0, 1000.0),
            threshold=750.0,
            max_slope=1000.0,
            max_replacements=0,
        )
        self.assertEqual(scaled.occupation_lower, base.occupation_lower)
        self.assertEqual(scaled.occupation_upper, base.occupation_upper)
        self.assertAlmostEqual(
            scaled.oxygen_deficit_lower,
            factor * base.oxygen_deficit_lower,
        )
        self.assertAlmostEqual(
            scaled.oxygen_deficit_upper,
            factor * base.oxygen_deficit_upper,
        )

    def test_power_of_two_metamorphic_scaling_with_one_replacement(self) -> None:
        times = (0.0, 0.75, 2.0, 4.0)
        values = (4.75, 100.0, 5.25, 4.875)
        threshold = 5.0
        max_slope = 0.75
        base = identify_hypoxia_burden_with_replacements(
            times,
            values,
            threshold=threshold,
            max_slope=max_slope,
            max_replacements=1,
        )
        base_values = (
            base.occupation_lower,
            base.occupation_upper,
            base.oxygen_deficit_lower,
            base.oxygen_deficit_upper,
        )

        for time_power, value_power in (
            (-80, 60),
            (-30, -70),
            (35, 25),
            (90, -20),
        ):
            time_factor = math.ldexp(1.0, time_power)
            value_factor = math.ldexp(1.0, value_power)
            transformed = identify_hypoxia_burden_with_replacements(
                tuple(time_factor * value for value in times),
                tuple(value_factor * value for value in values),
                threshold=value_factor * threshold,
                max_slope=value_factor * max_slope / time_factor,
                max_replacements=1,
            )
            transformed_values = (
                transformed.occupation_lower / time_factor,
                transformed.occupation_upper / time_factor,
                transformed.oxygen_deficit_lower
                / (time_factor * value_factor),
                transformed.oxygen_deficit_upper
                / (time_factor * value_factor),
            )
            for expected, observed in zip(base_values, transformed_values):
                self.assertTrue(
                    math.isclose(observed, expected, rel_tol=2e-13, abs_tol=0.0),
                    (time_power, value_power, expected, observed),
                )

    def test_vertical_translation_does_not_change_feasibility(self) -> None:
        self.assertEqual(k_min((0.0, 1.0), (0.0, 2.0), 1.0), 1)
        self.assertEqual(
            k_min((0.0, 1.0), (1e12, 1e12 + 2.0), 1.0), 1
        )
        result = audit_vertical_translation_invariance()
        self.assertLessEqual(
            float(result["maximum_endpoint_absolute_error"]), 1e-12
        )

    def test_independent_tent_integral(self) -> None:
        result = enumerate_endpoint_oracle(
            (0.0, 2.0),
            (0.0, 0.0),
            threshold=0.0,
            max_slope=1.0,
            max_replacements=0,
        )
        self.assertAlmostEqual(result.occupation_lower, 0.0)
        self.assertAlmostEqual(result.occupation_upper, 2.0)
        self.assertAlmostEqual(result.deficit_lower, 0.0)
        self.assertAlmostEqual(result.deficit_upper, 1.0)

    def test_nonconnected_image_has_correct_endpoint_hull(self) -> None:
        result = identify_hypoxia_burden_with_replacements(
            (0.0, 1.0),
            (-100.0, 100.0),
            threshold=0.0,
            max_slope=0.0,
            max_replacements=1,
        )
        self.assertEqual(result.oxygen_deficit_lower, 0.0)
        self.assertEqual(result.oxygen_deficit_upper, 100.0)
        self.assertFalse(math.isclose(50.0, 0.0) or math.isclose(50.0, 100.0))

    def test_random_four_endpoint_oracle(self) -> None:
        result = audit_random_cases(180, 314159)
        self.assertEqual(int(result["cases"]), 180)
        self.assertLessEqual(
            float(result["maximum_four_endpoint_absolute_error"]), 1e-8
        )

    def test_discrete_boundary_exhaustion(self) -> None:
        result = audit_discrete_boundary_cases()
        self.assertEqual(int(result["cases"]), 729)
        self.assertLessEqual(
            float(result["maximum_four_endpoint_absolute_error"]), 1e-10
        )


if __name__ == "__main__":
    unittest.main()
