# -*- coding: utf-8 -*-
"""Independent small-sample oracle for DIR-01 exact-inlier endpoints.

This module deliberately does not call the production feasibility,
subsequence, envelope, or integration helpers.  It enumerates retained sets,
checks all retained pairs, constructs global cone envelopes from affine-line
arrangements, and integrates threshold functionals directly.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import sys
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Literal, Sequence

from robust_exact_current import (
    InsufficientReplacementBudget,
    identify_hypoxia_burden_with_replacements,
    k_min,
)


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
Envelope = Literal["lower", "upper"]


class OracleInfeasible(ValueError):
    """Raised when no nonempty retained set is feasible within the budget."""


@dataclass(frozen=True)
class OracleBounds:
    occupation_lower: float
    occupation_upper: float
    deficit_lower: float
    deficit_upper: float
    occupation_lower_retained: tuple[int, ...]
    occupation_upper_retained: tuple[int, ...]
    deficit_lower_retained: tuple[int, ...]
    deficit_upper_retained: tuple[int, ...]


def _inputs(
    times: Iterable[float], values: Iterable[float], max_slope: float
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    t = tuple(float(value) for value in times)
    y = tuple(float(value) for value in values)
    if len(t) < 2 or len(y) != len(t):
        raise ValueError("times and values must have equal length at least two")
    if any(not math.isfinite(value) for value in (*t, *y)):
        raise ValueError("times and values must be finite")
    if any(right <= left for left, right in zip(t[:-1], t[1:])):
        raise ValueError("times must be strictly increasing")
    if not math.isfinite(max_slope) or max_slope < 0:
        raise ValueError("max_slope must be finite and non-negative")
    origin = t[0]
    return tuple(value - origin for value in t), y


def _all_pairs_feasible(
    times: Sequence[float],
    values: Sequence[float],
    retained: Sequence[int],
    max_slope: float,
) -> bool:
    """Check represented inputs with exact rational arithmetic.

    This deliberately avoids the production runtime's floating-point
    feasibility policy, so tiny but nonzero conflicts at ``L=0`` cannot be
    accepted by both implementations through a shared tolerance.
    """

    slope = Fraction.from_float(max_slope)
    for position, left in enumerate(retained):
        for right in retained[position + 1 :]:
            observed_difference = abs(
                Fraction.from_float(values[right])
                - Fraction.from_float(values[left])
            )
            allowed_difference = slope * (
                Fraction.from_float(times[right])
                - Fraction.from_float(times[left])
            )
            if observed_difference > allowed_difference:
                return False
    return True


def _cone_line(
    observation_time: float,
    observation_value: float,
    interval_midpoint: float,
    max_slope: float,
    envelope: Envelope,
) -> tuple[float, float]:
    left_of_observation = interval_midpoint < observation_time
    if envelope == "lower":
        if left_of_observation:
            return max_slope, observation_value - max_slope * observation_time
        return -max_slope, observation_value + max_slope * observation_time
    if left_of_observation:
        return -max_slope, observation_value + max_slope * observation_time
    return max_slope, observation_value - max_slope * observation_time


def _global_envelope_segments(
    times: Sequence[float],
    values: Sequence[float],
    retained: Sequence[int],
    max_slope: float,
    envelope: Envelope,
) -> list[tuple[float, float, float, float]]:
    """Return ``(left,right,slope,intercept)`` from all-cone arrangements."""

    segments: list[tuple[float, float, float, float]] = []
    choose_maximum = envelope == "lower"
    for left, right in zip(times[:-1], times[1:]):
        midpoint = 0.5 * (left + right)
        lines = [
            _cone_line(
                times[index],
                values[index],
                midpoint,
                max_slope,
                envelope,
            )
            for index in retained
        ]
        breakpoints = [left, right]
        for first_index, first in enumerate(lines):
            for second in lines[first_index + 1 :]:
                slope_difference = first[0] - second[0]
                if slope_difference == 0.0:
                    continue
                crossing = (second[1] - first[1]) / slope_difference
                if left < crossing < right:
                    breakpoints.append(crossing)
        breakpoints = sorted(set(breakpoints))
        for segment_left, segment_right in zip(
            breakpoints[:-1], breakpoints[1:]
        ):
            segment_midpoint = 0.5 * (segment_left + segment_right)
            selector = max if choose_maximum else min
            slope, intercept = selector(
                lines,
                key=lambda line: line[0] * segment_midpoint + line[1],
            )
            segments.append((segment_left, segment_right, slope, intercept))
    return segments


def _below_threshold_interval(
    left: float,
    right: float,
    slope: float,
    intercept: float,
    threshold: float,
) -> tuple[float, float] | None:
    if slope == 0.0:
        midpoint = 0.5 * (left + right)
        return (
            (left, right)
            if slope * midpoint + intercept < threshold
            else None
        )
    crossing = (threshold - intercept) / slope
    if slope > 0.0:
        below_right = min(right, crossing)
        return (left, below_right) if below_right > left else None
    below_left = max(left, crossing)
    return (below_left, right) if right > below_left else None


def _integrate_segments(
    segments: Sequence[tuple[float, float, float, float]], threshold: float
) -> tuple[float, float]:
    occupation = 0.0
    deficit = 0.0
    for left, right, slope, intercept in segments:
        below = _below_threshold_interval(
            left, right, slope, intercept, threshold
        )
        if below is None:
            continue
        active_left, active_right = below
        if active_right <= active_left:
            continue
        occupation += active_right - active_left
        deficit += (threshold - intercept) * (active_right - active_left)
        deficit -= 0.5 * slope * (
            active_right * active_right - active_left * active_left
        )
    return occupation, max(0.0, deficit)


def enumerate_endpoint_oracle(
    times: Iterable[float],
    values: Iterable[float],
    *,
    threshold: float,
    max_slope: float,
    max_replacements: int,
) -> OracleBounds:
    """Enumerate all retained sets and return four independent endpoints."""

    if not math.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if isinstance(max_replacements, bool) or not isinstance(max_replacements, int):
        raise TypeError("max_replacements must be an integer")
    if max_replacements < 0:
        raise ValueError("max_replacements must be non-negative")
    t, y = _inputs(times, values, max_slope)
    number = len(t)
    if max_replacements >= number:
        deleted = tuple(range(number))
        return OracleBounds(
            0.0,
            t[-1],
            0.0,
            math.inf,
            (),
            (),
            (),
            (),
        )

    best_values = [math.inf, -math.inf, math.inf, -math.inf]
    best_retained: list[tuple[int, ...]] = [(), (), (), ()]
    for deleted_count in range(max_replacements + 1):
        retained_count = number - deleted_count
        for retained in itertools.combinations(range(number), retained_count):
            if not _all_pairs_feasible(t, y, retained, max_slope):
                continue
            lower_segments = _global_envelope_segments(
                t, y, retained, max_slope, "lower"
            )
            upper_segments = _global_envelope_segments(
                t, y, retained, max_slope, "upper"
            )
            lower_occupation, lower_deficit = _integrate_segments(
                lower_segments, threshold
            )
            upper_occupation, upper_deficit = _integrate_segments(
                upper_segments, threshold
            )
            candidates = (
                upper_occupation,
                lower_occupation,
                upper_deficit,
                lower_deficit,
            )
            comparisons = (
                candidates[0] < best_values[0],
                candidates[1] > best_values[1],
                candidates[2] < best_values[2],
                candidates[3] > best_values[3],
            )
            for index, is_better in enumerate(comparisons):
                if is_better:
                    best_values[index] = candidates[index]
                    best_retained[index] = retained

    if not best_retained[0]:
        raise OracleInfeasible("no feasible nonempty retained set")
    return OracleBounds(
        occupation_lower=best_values[0],
        occupation_upper=best_values[1],
        deficit_lower=best_values[2],
        deficit_upper=best_values[3],
        occupation_lower_retained=best_retained[0],
        occupation_upper_retained=best_retained[1],
        deficit_lower_retained=best_retained[2],
        deficit_upper_retained=best_retained[3],
    )


def _compare_case(
    times: Sequence[float],
    values: Sequence[float],
    threshold: float,
    max_slope: float,
    max_replacements: int,
) -> tuple[bool, float]:
    oracle_error: Exception | None = None
    production_error: Exception | None = None
    oracle: OracleBounds | None = None
    production = None
    try:
        oracle = enumerate_endpoint_oracle(
            times,
            values,
            threshold=threshold,
            max_slope=max_slope,
            max_replacements=max_replacements,
        )
    except OracleInfeasible as error:
        oracle_error = error
    try:
        production = identify_hypoxia_burden_with_replacements(
            times,
            values,
            threshold=threshold,
            max_slope=max_slope,
            max_replacements=max_replacements,
        )
    except InsufficientReplacementBudget as error:
        production_error = error

    if oracle_error is not None or production_error is not None:
        if (oracle_error is None) != (production_error is None):
            raise AssertionError(
                "oracle/production feasibility disagreement: "
                f"times={times}, values={values}, L={max_slope}, "
                f"k={max_replacements}"
            )
        return False, 0.0

    assert oracle is not None and production is not None
    oracle_values = (
        oracle.occupation_lower,
        oracle.occupation_upper,
        oracle.deficit_lower,
        oracle.deficit_upper,
    )
    production_values = (
        production.occupation_lower,
        production.occupation_upper,
        production.oxygen_deficit_lower,
        production.oxygen_deficit_upper,
    )
    maximum_error = max(
        abs(left - right) for left, right in zip(oracle_values, production_values)
    )
    tolerance = 2e-9 * max(1.0, *(abs(value) for value in oracle_values))
    if maximum_error > tolerance:
        raise AssertionError(
            "oracle/production endpoint disagreement: "
            f"error={maximum_error}, oracle={oracle_values}, "
            f"production={production_values}, times={times}, values={values}, "
            f"H={threshold}, L={max_slope}, k={max_replacements}"
        )
    return True, maximum_error


def audit_random_cases(number_of_cases: int, seed: int) -> dict[str, float | int]:
    rng = random.Random(seed)
    feasible = 0
    infeasible = 0
    maximum_error = 0.0
    for _ in range(number_of_cases):
        number = rng.randint(2, 7)
        times = [0.0]
        for _index in range(1, number):
            times.append(times[-1] + rng.uniform(0.2, 1.5))
        max_slope = rng.choice((0.0, rng.uniform(0.1, 3.0)))
        budget = rng.randint(0, min(3, number - 1))
        if rng.random() < 0.6:
            values = [rng.uniform(-3.0, 3.0)]
            for left, right in zip(times[:-1], times[1:]):
                values.append(
                    values[-1]
                    + rng.uniform(-max_slope, max_slope) * (right - left)
                )
            for index in rng.sample(range(number), rng.randint(0, budget)):
                values[index] += rng.choice((-1.0, 1.0)) * rng.uniform(2.0, 12.0)
        else:
            values = [rng.uniform(-6.0, 6.0) for _ in range(number)]
        threshold = rng.uniform(min(values) - 1.0, max(values) + 1.0)
        is_feasible, error = _compare_case(
            times, values, threshold, max_slope, budget
        )
        if is_feasible:
            feasible += 1
            maximum_error = max(maximum_error, error)
        else:
            infeasible += 1
    return {
        "cases": number_of_cases,
        "feasible_cases": feasible,
        "infeasible_cases": infeasible,
        "maximum_four_endpoint_absolute_error": maximum_error,
    }


def audit_discrete_boundary_cases() -> dict[str, float | int]:
    checked = 0
    feasible = 0
    maximum_error = 0.0
    for max_slope in (0.0, 1.0, 2.0):
        for values in itertools.product((-1.0, 0.0, 1.0), repeat=3):
            for threshold in (-1.0, 0.0, 1.0):
                for budget in (0, 1, 2):
                    is_feasible, error = _compare_case(
                        (0.0, 1.0, 2.0),
                        values,
                        threshold,
                        max_slope,
                        budget,
                    )
                    checked += 1
                    feasible += int(is_feasible)
                    maximum_error = max(maximum_error, error)
    return {
        "cases": checked,
        "feasible_cases": feasible,
        "infeasible_cases": checked - feasible,
        "maximum_four_endpoint_absolute_error": maximum_error,
    }


def audit_vertical_translation_invariance() -> dict[str, object]:
    base_k_min = k_min((0.0, 1.0), (0.0, 2.0), 1.0)
    shifted_k_min = k_min((0.0, 1.0), (1e12, 1e12 + 2.0), 1.0)
    if base_k_min != 1 or shifted_k_min != 1:
        raise AssertionError(
            f"translation invariance failed: {base_k_min=} {shifted_k_min=}"
        )

    base = identify_hypoxia_burden_with_replacements(
        (0.0, 1.0, 2.0),
        (0.0, 0.5, 1.0),
        threshold=0.75,
        max_slope=1.0,
        max_replacements=0,
    )
    shift = 1e12
    translated = identify_hypoxia_burden_with_replacements(
        (0.0, 1.0, 2.0),
        (shift, shift + 0.5, shift + 1.0),
        threshold=shift + 0.75,
        max_slope=1.0,
        max_replacements=0,
    )
    base_endpoints = (
        base.occupation_lower,
        base.occupation_upper,
        base.oxygen_deficit_lower,
        base.oxygen_deficit_upper,
    )
    translated_endpoints = (
        translated.occupation_lower,
        translated.occupation_upper,
        translated.oxygen_deficit_lower,
        translated.oxygen_deficit_upper,
    )
    maximum_error = max(
        abs(left - right)
        for left, right in zip(base_endpoints, translated_endpoints)
    )
    if maximum_error > 1e-12:
        raise AssertionError(f"translated endpoint error too large: {maximum_error}")
    return {
        "base_k_min": base_k_min,
        "shifted_k_min": shifted_k_min,
        "vertical_shift": shift,
        "base_endpoints": list(base_endpoints),
        "shifted_endpoints": list(translated_endpoints),
        "maximum_endpoint_absolute_error": maximum_error,
    }


def audit_numerical_scale_regressions() -> dict[str, object]:
    """Exercise cases that absolute floating-point cutoffs previously broke."""

    production_tent = identify_hypoxia_burden_with_replacements(
        (0.0, 1e16),
        (0.0, 0.0),
        threshold=-0.25,
        max_slope=1e-16,
        max_replacements=0,
    )
    oracle_tent = enumerate_endpoint_oracle(
        (0.0, 1e16),
        (0.0, 0.0),
        threshold=-0.25,
        max_slope=1e-16,
        max_replacements=0,
    )
    expected_tent = (5e15, 6.25e14)
    observed_tent = (
        production_tent.occupation_upper,
        production_tent.oxygen_deficit_upper,
    )
    oracle_values = (
        oracle_tent.occupation_upper,
        oracle_tent.deficit_upper,
    )
    if observed_tent != expected_tent or oracle_values != expected_tent:
        raise AssertionError(
            "small-slope/large-horizon regression failed: "
            f"{observed_tent=} {oracle_values=}"
        )

    tiny_conflict_k_min = k_min((0.0, 1.0), (0.0, 1e-14), 0.0)
    oracle_rejected_tiny_conflict = False
    try:
        enumerate_endpoint_oracle(
            (0.0, 1.0),
            (0.0, 1e-14),
            threshold=0.0,
            max_slope=0.0,
            max_replacements=0,
        )
    except OracleInfeasible:
        oracle_rejected_tiny_conflict = True
    if tiny_conflict_k_min != 1 or not oracle_rejected_tiny_conflict:
        raise AssertionError(
            "tiny constant-path conflict was hidden by a numerical tolerance"
        )

    one_ulp_k_min = k_min(
        (0.0, 1.0),
        (0.0, math.nextafter(1.0, math.inf)),
        1.0,
    )
    rounded_subtraction_k_min = k_min(
        (0.0, 1.0),
        (1.0, -math.ldexp(1.0, -54)),
        1.0,
    )
    overflow_comparison_k_min = k_min(
        (0.0, math.nextafter(2.0, 0.0)),
        (-sys.float_info.max, sys.float_info.max),
        sys.float_info.max,
    )
    if (one_ulp_k_min, rounded_subtraction_k_min, overflow_comparison_k_min) != (
        1,
        1,
        1,
    ):
        raise AssertionError(
            "exact binary64 compatibility regressions failed: "
            f"{one_ulp_k_min=}, {rounded_subtraction_k_min=}, "
            f"{overflow_comparison_k_min=}"
        )

    alias_reports = (1e12, 0.0, 1.0, math.nextafter(1.0, math.inf))
    alias_k_min = k_min((0.0, 1.0, 2.0, 3.0), alias_reports, 0.0)
    alias_endpoint_rejected = False
    try:
        identify_hypoxia_burden_with_replacements(
            (0.0, 1.0, 2.0, 3.0),
            alias_reports,
            threshold=0.5,
            max_slope=0.0,
            max_replacements=3,
        )
    except FloatingPointError:
        alias_endpoint_rejected = True
    if alias_k_min != 3 or not alias_endpoint_rejected:
        raise AssertionError(
            "non-injective normalization was not separated from exact k_min"
        )

    tiny_slope_endpoint_rejected = False
    try:
        identify_hypoxia_burden_with_replacements(
            (0.0, 1.0),
            (0.0, 1e308),
            threshold=0.0,
            max_slope=math.ulp(0.0),
            max_replacements=1,
        )
    except FloatingPointError:
        tiny_slope_endpoint_rejected = True
    if not tiny_slope_endpoint_rejected:
        raise AssertionError(
            "normalization silently erased a positive Lipschitz slope"
        )

    near_vertex = identify_hypoxia_burden_with_replacements(
        (0.0, 1e50),
        (0.0006197765650313158, -23.482067348105762),
        threshold=-24.18785045032845,
        max_slope=2.4894253329116173e-49,
        max_replacements=0,
    )
    near_vertex_values = (
        near_vertex.occupation_upper,
        near_vertex.oxygen_deficit_upper,
    )
    if near_vertex_values[0] != 1.493678147677149e34 or not math.isclose(
        near_vertex_values[1],
        1.3885232882494042e19,
        rel_tol=2e-16,
        abs_tol=0.0,
    ):
        raise AssertionError(
            "exact local fallback did not preserve a near-vertex wedge"
        )

    no_anchor = identify_hypoxia_burden_with_replacements(
        (0.0, 1e-308),
        (0.0, 0.0),
        threshold=1e308,
        max_slope=0.0,
        max_replacements=2,
        state_lower_bound=-1e308,
    )
    if not math.isclose(
        no_anchor.oxygen_deficit_upper,
        2.0,
        rel_tol=2e-16,
        abs_tol=0.0,
    ):
        raise AssertionError(
            "no-anchor exact-product regression failed: "
            f"{no_anchor.oxygen_deficit_upper}"
        )

    base = identify_hypoxia_burden_with_replacements(
        (0.0, 0.75, 2.0, 4.0),
        (4.75, 100.0, 5.25, 4.875),
        threshold=5.0,
        max_slope=0.75,
        max_replacements=1,
    )
    base_values = (
        base.occupation_lower,
        base.occupation_upper,
        base.oxygen_deficit_lower,
        base.oxygen_deficit_upper,
    )
    maximum_normalized_error = 0.0
    transformations = []
    for time_power, value_power in (
        (-80, 60),
        (-30, -70),
        (35, 25),
        (90, -20),
    ):
        time_factor = math.ldexp(1.0, time_power)
        value_factor = math.ldexp(1.0, value_power)
        transformed = identify_hypoxia_burden_with_replacements(
            tuple(time_factor * value for value in (0.0, 0.75, 2.0, 4.0)),
            tuple(
                value_factor * value
                for value in (4.75, 100.0, 5.25, 4.875)
            ),
            threshold=value_factor * 5.0,
            max_slope=value_factor * 0.75 / time_factor,
            max_replacements=1,
        )
        normalized_values = (
            transformed.occupation_lower / time_factor,
            transformed.occupation_upper / time_factor,
            transformed.oxygen_deficit_lower
            / (time_factor * value_factor),
            transformed.oxygen_deficit_upper
            / (time_factor * value_factor),
        )
        error = max(
            abs(left - right)
            for left, right in zip(base_values, normalized_values)
        )
        maximum_normalized_error = max(maximum_normalized_error, error)
        transformations.append(
            {
                "time_power_of_two": time_power,
                "value_power_of_two": value_power,
                "maximum_normalized_endpoint_error": error,
            }
        )
    if maximum_normalized_error > 2e-12:
        raise AssertionError(
            f"power-of-two metamorphic error too large: "
            f"{maximum_normalized_error}"
        )
    return {
        "small_slope_large_horizon_expected": list(expected_tent),
        "small_slope_large_horizon_production": list(observed_tent),
        "small_slope_large_horizon_oracle": list(oracle_values),
        "tiny_L0_conflict_production_k_min": tiny_conflict_k_min,
        "tiny_L0_conflict_oracle_rejected": oracle_rejected_tiny_conflict,
        "one_ulp_violation_k_min": one_ulp_k_min,
        "rounded_subtraction_violation_k_min": rounded_subtraction_k_min,
        "overflow_comparison_violation_k_min": overflow_comparison_k_min,
        "normalization_alias_k_min": alias_k_min,
        "normalization_alias_endpoint_rejected": alias_endpoint_rejected,
        "tiny_positive_slope_endpoint_rejected": (
            tiny_slope_endpoint_rejected
        ),
        "near_vertex_exact_fallback_endpoints": list(near_vertex_values),
        "no_anchor_exact_product_deficit_upper": (
            no_anchor.oxygen_deficit_upper
        ),
        "power_of_two_transformations": transformations,
        "maximum_power_of_two_normalized_endpoint_error": (
            maximum_normalized_error
        ),
    }


def run_oracle_gate(random_cases: int, seed: int) -> dict[str, object]:
    tent = enumerate_endpoint_oracle(
        (0.0, 2.0),
        (0.0, 0.0),
        threshold=0.0,
        max_slope=1.0,
        max_replacements=0,
    )
    if not (
        math.isclose(tent.deficit_lower, 0.0, abs_tol=1e-12)
        and math.isclose(tent.deficit_upper, 1.0, abs_tol=1e-12)
        and math.isclose(tent.occupation_lower, 0.0, abs_tol=1e-12)
        and math.isclose(tent.occupation_upper, 2.0, abs_tol=1e-12)
    ):
        raise AssertionError(f"tent primitive failed: {tent}")
    return {
        "gate": "DIR-01-G3-independent-four-endpoint-oracle",
        "status": "pass",
        "independence": (
            "enumerated retained sets; all-pairs feasibility; all-cone affine "
            "arrangement; direct occupation/hinge integration"
        ),
        "production_helpers_not_called": [
            "_as_inputs",
            "_compatible",
            "_retained_subset_value",
            "_pair_segments",
            "integrate_occupation",
            "integrate_oxygen_deficit",
        ],
        "random_audit": audit_random_cases(random_cases, seed),
        "discrete_boundary_audit": audit_discrete_boundary_cases(),
        "vertical_translation_invariance": audit_vertical_translation_invariance(),
        "numerical_scale_regressions": audit_numerical_scale_regressions(),
        "analytic_tent_case": asdict(tent),
        "scope": "exact-inlier epsilon=0 only; endpoint range/hull, not image connectedness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--random-cases", type=int, default=1080)
    parser.add_argument("--seed", type=int, default=20260723)
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT
        / "03_代码与测试"
        / "gate_outputs"
        / "g3_independent_endpoint_oracle.json",
    )
    args = parser.parse_args()
    result = run_oracle_gate(args.random_cases, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
