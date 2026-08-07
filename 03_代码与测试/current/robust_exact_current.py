# -*- coding: utf-8 -*-
"""Conditioned numerical runtime for exact-inlier endpoint identification.

The byte-for-byte historical implementation under ``原始快照`` is retained
for provenance.  This module exposes the current public API while reusing its
proved dynamic-programming structure.

Pairwise inlier feasibility is evaluated exactly for the binary64 values
supplied by the caller. Endpoint quadrature is evaluated after a conditioned
affine normalization. Inputs whose normalization would erase distinctions are
rejected explicitly instead of returning silently incorrect endpoints.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
SNAPSHOT = HERE.parent / "原始快照"
if str(SNAPSHOT) not in sys.path:
    sys.path.insert(0, str(SNAPSHOT))

import future_partial_identification_robust as _snapshot  # noqa: E402


NUMERICAL_POLICY = (
    "exact rational compatibility for represented binary64 inputs; "
    "injective, condition-guarded affine normalization for endpoint quadrature; "
    "exact Fraction fallback for geometry-critical local costs; "
    "safe rational rescaling; "
    "direct DP ordering with deterministic first-encounter tie handling"
)

NORMALIZATION_CONDITION_LIMIT = 1 << 40
CRITICAL_GEOMETRY_FALLBACK_LIMIT = 1 << 24


@dataclass(frozen=True)
class _ExactCostContext:
    times: tuple[Fraction, ...]
    reports: tuple[Fraction, ...]
    max_slope: Fraction
    threshold: Fraction
    floor: Fraction | None
    horizon: Fraction
    value_scale: Fraction
    boundary_fallback: np.ndarray
    pair_fallback: np.ndarray


def _validated_arrays(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
) -> tuple[np.ndarray, np.ndarray]:
    raw_times = np.asarray(list(times), dtype=float)
    observations = np.asarray(list(values), dtype=float)
    if raw_times.ndim != 1 or observations.ndim != 1:
        raise ValueError("times and values must be one-dimensional")
    if raw_times.size < 2:
        raise ValueError("at least two observation times are required")
    if observations.size != raw_times.size:
        raise ValueError("times and values must have equal lengths")
    if not np.all(np.isfinite(raw_times)) or not np.all(np.isfinite(observations)):
        raise ValueError("times and values must be finite")
    if np.any(np.diff(raw_times) <= 0):
        raise ValueError("times must be strictly increasing")
    if not math.isfinite(max_slope) or max_slope < 0:
        raise ValueError("max_slope must be finite and non-negative")
    return raw_times, observations


def _as_fraction(value: float) -> Fraction:
    """Interpret a finite binary64 value as its exact dyadic rational."""

    return Fraction.from_float(float(value))


def _normalization_condition_guard(
    exact_times: Sequence[Fraction],
    exact_values: Sequence[Fraction],
    *,
    center: Fraction,
    horizon: Fraction,
    value_scale: Fraction,
    slope_span: Fraction,
) -> None:
    """Reject affine maps that cannot preserve relevant distinctions."""

    minimum_time_gap = min(
        exact_times[index + 1] - exact_times[index]
        for index in range(len(exact_times) - 1)
    )
    if horizon / minimum_time_gap > NORMALIZATION_CONDITION_LIMIT:
        raise FloatingPointError(
            "time normalization is ill-conditioned: horizon/minimum gap "
            f"exceeds {NORMALIZATION_CONDITION_LIMIT}"
        )

    normalized_times = [
        float((value - exact_times[0]) / horizon) for value in exact_times
    ]
    if any(
        right <= left
        for left, right in zip(normalized_times[:-1], normalized_times[1:])
    ):
        raise FloatingPointError(
            "time normalization is non-injective for the supplied binary64 "
            "observation times"
        )

    if (
        slope_span > 0
        and value_scale / slope_span > NORMALIZATION_CONDITION_LIMIT
    ):
        raise FloatingPointError(
            "value normalization is ill-conditioned: scale/(L*horizon) "
            f"exceeds {NORMALIZATION_CONDITION_LIMIT}"
        )

    unique_values = sorted(set(exact_values))
    if len(unique_values) < 2:
        return
    minimum_value_gap = min(
        right - left
        for left, right in zip(unique_values[:-1], unique_values[1:])
    )
    if value_scale / minimum_value_gap > NORMALIZATION_CONDITION_LIMIT:
        raise FloatingPointError(
            "value normalization is ill-conditioned: scale/minimum relevant "
            f"gap exceeds {NORMALIZATION_CONDITION_LIMIT}"
        )
    normalized_values = [
        float((value - center) / value_scale) for value in unique_values
    ]
    if any(
        right <= left
        for left, right in zip(
            normalized_values[:-1], normalized_values[1:]
        )
    ):
        raise FloatingPointError(
            "value normalization is non-injective for reports, threshold, "
            "or state lower bound"
        )


def _needs_exact_fallback(
    value_scale: Fraction,
    gap: Fraction,
) -> bool:
    """Return whether a critical event needs exact local integration."""

    absolute_gap = abs(gap)
    return bool(
        absolute_gap == 0
        or value_scale / absolute_gap > CRITICAL_GEOMETRY_FALLBACK_LIMIT
    )


def _derived_geometry_fallback_masks(
    exact_times: Sequence[Fraction],
    exact_reports: Sequence[Fraction],
    valid_reports: np.ndarray,
    compatibility: np.ndarray,
    *,
    exact_slope: Fraction,
    exact_threshold: Fraction,
    exact_floor: Fraction | None,
    value_scale: Fraction,
) -> tuple[np.ndarray, np.ndarray]:
    """Flag local costs whose critical events need exact integration."""

    number = len(exact_times)
    boundary_fallback = np.zeros((number, 2), dtype=bool)
    pair_fallback = np.zeros((number, number), dtype=bool)
    levels = [exact_threshold]
    if exact_floor is not None:
        levels.append(exact_floor)
    valid_indices = [
        index for index, is_valid in enumerate(valid_reports) if bool(is_valid)
    ]
    first_time = exact_times[0]
    last_time = exact_times[-1]

    for index in valid_indices:
        if any(
            _needs_exact_fallback(
                value_scale, level - exact_reports[index]
            )
            for level in levels
        ):
            boundary_fallback[index, :] = True
            for other in valid_indices:
                left, right = sorted((index, other))
                if left < right and bool(compatibility[left, right]):
                    pair_fallback[left, right] = True
        left_distance = exact_times[index] - first_time
        right_distance = last_time - exact_times[index]
        for side_index, distance in enumerate(
            (left_distance, right_distance)
        ):
            for sign in (-1, 1):
                boundary_value = (
                    exact_reports[index] + sign * exact_slope * distance
                )
                for level in levels:
                    if _needs_exact_fallback(
                        value_scale, level - boundary_value
                    ):
                        boundary_fallback[index, side_index] = True

    for position, left in enumerate(valid_indices):
        for right in valid_indices[position + 1 :]:
            if not bool(compatibility[left, right]):
                continue
            time_difference = exact_times[right] - exact_times[left]
            allowed = exact_slope * time_difference
            observed = abs(exact_reports[right] - exact_reports[left])
            if _needs_exact_fallback(value_scale, allowed - observed):
                pair_fallback[left, right] = True
            lower_vertex = (
                exact_reports[left] + exact_reports[right] - allowed
            ) / 2
            upper_vertex = (
                exact_reports[left] + exact_reports[right] + allowed
            ) / 2
            for level in levels:
                if _needs_exact_fallback(
                    value_scale, level - lower_vertex
                ) or _needs_exact_fallback(
                    value_scale, level - upper_vertex
                ):
                    pair_fallback[left, right] = True
    return boundary_fallback, pair_fallback


def _normalized_problem(
    times: np.ndarray,
    values: np.ndarray,
    valid_reports: np.ndarray,
    compatibility: np.ndarray,
    max_slope: float,
    *,
    threshold: float,
    state_lower_bound: float | None,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
    Fraction,
    Fraction,
    Fraction,
    float,
    float,
    float | None,
    _ExactCostContext,
]:
    """Return conditioned normalized inputs and exact inverse-scale factors."""

    exact_times = tuple(_as_fraction(value) for value in times)
    exact_reports = tuple(_as_fraction(value) for value in values)
    exact_slope = _as_fraction(max_slope)
    exact_threshold = _as_fraction(threshold)
    exact_floor = (
        None
        if state_lower_bound is None
        else _as_fraction(state_lower_bound)
    )
    horizon_exact = exact_times[-1] - exact_times[0]
    try:
        horizon = float(horizon_exact)
    except OverflowError as error:
        raise OverflowError(
            "time horizon is not representable as a finite float"
        ) from error
    if not math.isfinite(horizon) or horizon <= 0.0:
        raise OverflowError("time horizon is not representable as a finite float")

    valid_exact_reports = [
        value
        for value, is_valid in zip(exact_reports, valid_reports)
        if bool(is_valid)
    ]
    if not valid_exact_reports:
        raise ValueError("normalization requires at least one valid report")
    # Center at the functional threshold so strict crossing calculations use
    # an exactly represented normalized level of zero.
    center_exact = exact_threshold
    relevant_values = list(valid_exact_reports)
    relevant_values.append(exact_threshold)
    if exact_floor is not None:
        relevant_values.append(exact_floor)
    slope_span_exact = exact_slope * horizon_exact
    value_scale_exact = max(
        [abs(value - center_exact) for value in relevant_values]
        + [slope_span_exact]
    )
    if value_scale_exact == 0:
        value_scale_exact = Fraction(1)

    _normalization_condition_guard(
        exact_times,
        relevant_values,
        center=center_exact,
        horizon=horizon_exact,
        value_scale=value_scale_exact,
        slope_span=slope_span_exact,
    )
    boundary_fallback, pair_fallback = _derived_geometry_fallback_masks(
        exact_times,
        exact_reports,
        valid_reports,
        compatibility,
        exact_slope=exact_slope,
        exact_threshold=exact_threshold,
        exact_floor=exact_floor,
        value_scale=value_scale_exact,
    )
    exact_cost_context = _ExactCostContext(
        times=exact_times,
        reports=exact_reports,
        max_slope=exact_slope,
        threshold=exact_threshold,
        floor=exact_floor,
        horizon=horizon_exact,
        value_scale=value_scale_exact,
        boundary_fallback=boundary_fallback,
        pair_fallback=pair_fallback,
    )
    normalized_times = np.asarray(
        [
            float((value - exact_times[0]) / horizon_exact)
            for value in exact_times
        ],
        dtype=float,
    )
    normalized_values = np.asarray(
        [
            (
                float((value - center_exact) / value_scale_exact)
                if bool(is_valid)
                else 0.0
            )
            for value, is_valid in zip(exact_reports, valid_reports)
        ],
        dtype=float,
    )
    normalized_slope = float(slope_span_exact / value_scale_exact)
    if slope_span_exact > 0 and normalized_slope == 0.0:
        raise FloatingPointError(
            "positive max_slope underflows during endpoint normalization"
        )
    normalized_threshold = float(
        (exact_threshold - center_exact) / value_scale_exact
    )
    normalized_floor = (
        None
        if exact_floor is None
        else float((exact_floor - center_exact) / value_scale_exact)
    )
    return (
        normalized_times,
        normalized_values,
        horizon,
        horizon_exact,
        value_scale_exact,
        center_exact,
        normalized_slope,
        normalized_threshold,
        normalized_floor,
        exact_cost_context,
    )


@lru_cache(maxsize=32)
def _exact_compatibility_bytes(
    times: tuple[float, ...],
    values: tuple[float, ...],
    max_slope: float,
) -> bytes:
    """Cache the exact compatibility DAG for represented binary64 inputs."""

    exact_times = tuple(_as_fraction(value) for value in times)
    exact_values = tuple(_as_fraction(value) for value in values)
    exact_slope = _as_fraction(max_slope)
    number = len(times)
    matrix = bytearray(number * number)
    for index in range(number):
        matrix[index * number + index] = 1
    for left in range(number):
        for right in range(left + 1, number):
            observed = abs(exact_values[right] - exact_values[left])
            allowed = exact_slope * (exact_times[right] - exact_times[left])
            if observed <= allowed:
                matrix[left * number + right] = 1
                matrix[right * number + left] = 1
    return bytes(matrix)


def _exact_compatibility_matrix(
    times: np.ndarray,
    values: np.ndarray,
    max_slope: float,
) -> np.ndarray:
    packed = _exact_compatibility_bytes(
        tuple(float(value) for value in times),
        tuple(float(value) for value in values),
        float(max_slope),
    )
    number = int(times.size)
    return np.frombuffer(packed, dtype=np.uint8).reshape(number, number)


def _selected_two_line_segments(
    start: float,
    end: float,
    first: tuple[float, float],
    second: tuple[float, float],
    *,
    choose_maximum: bool,
) -> list[object]:
    """Return the max/min of two affine lines without absolute cutoffs."""

    slope_a, intercept_a = first
    slope_b, intercept_b = second
    points = [start, end]
    slope_difference = slope_a - slope_b
    if slope_difference != 0.0:
        crossing = (intercept_b - intercept_a) / slope_difference
        if start < crossing < end:
            points.insert(1, crossing)

    output: list[object] = []
    for left, right in zip(points[:-1], points[1:]):
        midpoint = 0.5 * (left + right)
        value_a = slope_a * midpoint + intercept_a
        value_b = slope_b * midpoint + intercept_b
        use_a = value_a >= value_b if choose_maximum else value_a <= value_b
        slope, intercept = first if use_a else second
        output.append(_snapshot.AffineSegment(left, right, slope, intercept))
    return output


def _boundary_segments(
    times: np.ndarray,
    values: np.ndarray,
    index: int,
    max_slope: float,
    envelope: str,
    side: str,
) -> list[object]:
    observation_time = float(times[index])
    observation_value = float(values[index])
    if side == "left":
        start, end = 0.0, observation_time
        slope = max_slope if envelope == "lower" else -max_slope
    else:
        start, end = observation_time, float(times[-1])
        slope = -max_slope if envelope == "lower" else max_slope
    if end <= start:
        return []
    intercept = observation_value - slope * observation_time
    return [_snapshot.AffineSegment(start, end, slope, intercept)]


def _pair_segments(
    times: np.ndarray,
    values: np.ndarray,
    left: int,
    right: int,
    max_slope: float,
    envelope: str,
) -> list[object]:
    left_time = float(times[left])
    right_time = float(times[right])
    left_value = float(values[left])
    right_value = float(values[right])
    if envelope == "lower":
        first = (-max_slope, left_value + max_slope * left_time)
        second = (max_slope, right_value - max_slope * right_time)
        choose_maximum = True
    else:
        first = (max_slope, left_value - max_slope * left_time)
        second = (-max_slope, right_value + max_slope * right_time)
        choose_maximum = False
    return _selected_two_line_segments(
        left_time,
        right_time,
        first,
        second,
        choose_maximum=choose_maximum,
    )


def _below_subinterval(segment: object, threshold: float) -> tuple[float, float] | None:
    start = float(segment.start)
    end = float(segment.end)
    slope = float(segment.slope)
    intercept = float(segment.intercept)
    if slope == 0.0:
        value = slope * (0.5 * (start + end)) + intercept
        return (start, end) if value < threshold else None

    crossing = (threshold - intercept) / slope
    if slope > 0.0:
        below_end = min(end, crossing)
        return (start, below_end) if below_end > start else None
    below_start = max(start, crossing)
    return (below_start, end) if end > below_start else None


def _integrate_occupation(segments: Sequence[object], threshold: float) -> float:
    lengths = []
    for segment in segments:
        subinterval = _below_subinterval(segment, threshold)
        if subinterval is not None:
            lengths.append(subinterval[1] - subinterval[0])
    return float(math.fsum(lengths))


def _integrate_deficit(segments: Sequence[object], threshold: float) -> float:
    contributions = []
    for segment in segments:
        subinterval = _below_subinterval(segment, threshold)
        if subinterval is None:
            continue
        start, end = subinterval
        value_start = float(segment.slope) * start + float(segment.intercept)
        value_end = float(segment.slope) * end + float(segment.intercept)
        contribution = (
            0.5
            * ((threshold - value_start) + (threshold - value_end))
            * (end - start)
        )
        contributions.append(max(0.0, contribution))
    return float(max(0.0, math.fsum(contributions)))


def _is_better(candidate: float, incumbent: float, maximize: bool) -> bool:
    """Order representable objectives directly; exact ties keep first state."""

    return candidate > incumbent if maximize else candidate < incumbent


def _clip_segments_to_floor(
    segments: Sequence[object],
    floor: float,
) -> list[object]:
    """Return affine pieces for the pointwise maximum with a constant floor."""

    clipped: list[object] = []
    for segment in segments:
        clipped.extend(
            _selected_two_line_segments(
                float(segment.start),
                float(segment.end),
                (float(segment.slope), float(segment.intercept)),
                (0.0, floor),
                choose_maximum=True,
            )
        )
    return clipped


def _envelope_integral_cost(
    segments: Sequence[object],
    threshold: float,
    floor: float | None,
    *,
    envelope: str,
    functional: str,
) -> float:
    constrained_segments = (
        _clip_segments_to_floor(segments, floor)
        if envelope == "lower" and floor is not None
        else segments
    )
    if functional == "occupation":
        return _integrate_occupation(constrained_segments, threshold)
    return _integrate_deficit(constrained_segments, threshold)


def _exact_boundary_segments(
    context: _ExactCostContext,
    index: int,
    envelope: str,
    side: str,
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    observation_time = context.times[index]
    observation_value = context.reports[index]
    if side == "left":
        start, end = context.times[0], observation_time
        slope = (
            context.max_slope
            if envelope == "lower"
            else -context.max_slope
        )
    else:
        start, end = observation_time, context.times[-1]
        slope = (
            -context.max_slope
            if envelope == "lower"
            else context.max_slope
        )
    if end <= start:
        return []
    intercept = observation_value - slope * observation_time
    return [(start, end, slope, intercept)]


def _exact_pair_segments(
    context: _ExactCostContext,
    left: int,
    right: int,
    envelope: str,
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    left_time = context.times[left]
    right_time = context.times[right]
    left_value = context.reports[left]
    right_value = context.reports[right]
    slope_bound = context.max_slope
    if slope_bound == 0:
        return [(left_time, right_time, Fraction(0), left_value)]
    if envelope == "lower":
        first = (-slope_bound, left_value + slope_bound * left_time)
        second = (slope_bound, right_value - slope_bound * right_time)
    else:
        first = (slope_bound, left_value - slope_bound * left_time)
        second = (-slope_bound, right_value + slope_bound * right_time)
    crossing = (second[1] - first[1]) / (first[0] - second[0])
    if not left_time <= crossing <= right_time:
        raise FloatingPointError(
            "exact compatibility and pair-envelope geometry disagree"
        )
    output: list[tuple[Fraction, Fraction, Fraction, Fraction]] = []
    if crossing > left_time:
        output.append((left_time, crossing, first[0], first[1]))
    if right_time > crossing:
        output.append((crossing, right_time, second[0], second[1]))
    return output


def _exact_clip_segments_to_floor(
    segments: Sequence[tuple[Fraction, Fraction, Fraction, Fraction]],
    floor: Fraction,
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    clipped: list[tuple[Fraction, Fraction, Fraction, Fraction]] = []
    for start, end, slope, intercept in segments:
        points = [start, end]
        if slope != 0:
            crossing = (floor - intercept) / slope
            if start < crossing < end:
                points.insert(1, crossing)
        for left, right in zip(points[:-1], points[1:]):
            midpoint = (left + right) / 2
            if slope * midpoint + intercept >= floor:
                clipped.append((left, right, slope, intercept))
            else:
                clipped.append((left, right, Fraction(0), floor))
    return clipped


def _exact_segment_below_interval(
    segment: tuple[Fraction, Fraction, Fraction, Fraction],
    threshold: Fraction,
) -> tuple[Fraction, Fraction] | None:
    start, end, slope, intercept = segment
    if slope == 0:
        return (start, end) if intercept < threshold else None
    crossing = (threshold - intercept) / slope
    if slope > 0:
        below_end = min(end, crossing)
        return (start, below_end) if below_end > start else None
    below_start = max(start, crossing)
    return (below_start, end) if end > below_start else None


def _exact_local_integral_cost(
    context: _ExactCostContext,
    segments: Sequence[tuple[Fraction, Fraction, Fraction, Fraction]],
    *,
    envelope: str,
    functional: str,
) -> float:
    constrained = (
        _exact_clip_segments_to_floor(segments, context.floor)
        if envelope == "lower" and context.floor is not None
        else list(segments)
    )
    total = Fraction(0)
    for segment in constrained:
        subinterval = _exact_segment_below_interval(
            segment, context.threshold
        )
        if subinterval is None:
            continue
        start, end = subinterval
        if functional == "occupation":
            total += end - start
            continue
        slope = segment[2]
        intercept = segment[3]
        value_start = slope * start + intercept
        value_end = slope * end + intercept
        total += (
            (2 * context.threshold - value_start - value_end)
            * (end - start)
            / 2
        )
    divisor = (
        context.horizon
        if functional == "occupation"
        else context.horizon * context.value_scale
    )
    return _safe_fraction_to_float(total / divisor)


def _minimum_replacements_from_matrix(
    compatibility: np.ndarray,
    valid_reports: np.ndarray,
) -> object:
    """Longest compatible subsequence among reports allowed as anchors."""

    number = int(valid_reports.size)
    lengths = [0] * number
    parents: list[int | None] = [None] * number
    for right in range(number):
        if not bool(valid_reports[right]):
            continue
        lengths[right] = 1
        for left in range(right):
            if lengths[left] == 0:
                continue
            if not bool(compatibility[left, right]):
                continue
            candidate = lengths[left] + 1
            if candidate > lengths[right]:
                lengths[right] = candidate
                parents[right] = left

    best_length = max(lengths, default=0)
    if best_length == 0:
        retained: tuple[int, ...] = ()
    else:
        best_end = next(
            index for index, length in enumerate(lengths) if length == best_length
        )
        retained_reversed: list[int] = []
        cursor: int | None = best_end
        while cursor is not None:
            retained_reversed.append(cursor)
            cursor = parents[cursor]
        retained = tuple(reversed(retained_reversed))
    retained_set = set(retained)
    deleted = tuple(index for index in range(number) if index not in retained_set)
    return MinimumReplacementResult(
        minimum_replacements=number - len(retained),
        retained_indices=retained,
        deleted_indices=deleted,
    )


def _optimize_subsequence(
    times: np.ndarray,
    values: np.ndarray,
    compatibility: np.ndarray,
    valid_reports: np.ndarray,
    max_slope: float,
    threshold: float,
    floor: float | None,
    exact_context: _ExactCostContext,
    max_replacements: int,
    *,
    envelope: str,
    functional: str,
    maximize: bool,
) -> object:
    """Conditioned ``O(k n^2)`` endpoint dynamic program."""

    number = int(times.size)
    budget = min(max_replacements, number - 1)
    scores: dict[tuple[int, int], float] = {}
    parents: dict[tuple[int, int], tuple[int, int] | None] = {}

    for right in range(number):
        if not bool(valid_reports[right]):
            continue
        if right <= budget:
            state = (right, right)
            if bool(exact_context.boundary_fallback[right, 0]):
                scores[state] = _exact_local_integral_cost(
                    exact_context,
                    _exact_boundary_segments(
                        exact_context, right, envelope, "left"
                    ),
                    envelope=envelope,
                    functional=functional,
                )
            else:
                scores[state] = _envelope_integral_cost(
                    _boundary_segments(
                        times, values, right, max_slope, envelope, "left"
                    ),
                    threshold,
                    floor,
                    envelope=envelope,
                    functional=functional,
                )
            parents[state] = None

        for left in range(right):
            if not bool(compatibility[left, right]):
                continue
            skipped = right - left - 1
            if bool(exact_context.pair_fallback[left, right]):
                edge_cost = _exact_local_integral_cost(
                    exact_context,
                    _exact_pair_segments(
                        exact_context, left, right, envelope
                    ),
                    envelope=envelope,
                    functional=functional,
                )
            else:
                edge_cost = _envelope_integral_cost(
                    _pair_segments(
                        times, values, left, right, max_slope, envelope
                    ),
                    threshold,
                    floor,
                    envelope=envelope,
                    functional=functional,
                )
            for deleted_before in range(budget - skipped + 1):
                previous_state = (deleted_before, left)
                if previous_state not in scores:
                    continue
                deleted_now = deleted_before + skipped
                state = (deleted_now, right)
                candidate = scores[previous_state] + edge_cost
                if state not in scores or _is_better(
                    candidate, scores[state], maximize
                ):
                    scores[state] = candidate
                    parents[state] = previous_state

    best_value = -math.inf if maximize else math.inf
    best_state: tuple[int, int] | None = None
    for state, prefix_value in scores.items():
        deleted_before, last = state
        if deleted_before + (number - last - 1) > budget:
            continue
        if bool(exact_context.boundary_fallback[last, 1]):
            boundary_cost = _exact_local_integral_cost(
                exact_context,
                _exact_boundary_segments(
                    exact_context, last, envelope, "right"
                ),
                envelope=envelope,
                functional=functional,
            )
        else:
            boundary_cost = _envelope_integral_cost(
                _boundary_segments(
                    times, values, last, max_slope, envelope, "right"
                ),
                threshold,
                floor,
                envelope=envelope,
                functional=functional,
            )
        value = prefix_value + boundary_cost
        if best_state is None or _is_better(value, best_value, maximize):
            best_value = value
            best_state = state

    if best_state is None:
        raise InsufficientReplacementBudget(
            "no nonempty valid Lipschitz retained subset fits the replacement "
            "budget"
        )

    retained_reversed: list[int] = []
    cursor: tuple[int, int] | None = best_state
    while cursor is not None:
        retained_reversed.append(cursor[1])
        cursor = parents[cursor]
    retained = tuple(reversed(retained_reversed))
    retained_set = set(retained)
    deleted = tuple(index for index in range(number) if index not in retained_set)
    return ExtremumWitness(
        value=float(best_value),
        retained_indices=retained,
        deleted_indices=deleted,
        attained=True,
    )


InsufficientReplacementBudget = _snapshot.InsufficientReplacementBudget
ExtremumWitness = _snapshot.ExtremumWitness
MinimumReplacementResult = _snapshot.MinimumReplacementResult
RobustBurdenBounds = _snapshot.RobustBurdenBounds


def minimum_replacements_for_lipschitz(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
    *,
    state_lower_bound: float | None = None,
) -> MinimumReplacementResult:
    """Return exact ``k_min`` for the represented binary64 input values."""

    raw_times, observations = _validated_arrays(times, values, max_slope)
    if state_lower_bound is not None and not math.isfinite(state_lower_bound):
        raise ValueError("state_lower_bound must be finite")
    compatibility = _exact_compatibility_matrix(
        raw_times, observations, max_slope
    )
    valid_reports = (
        np.ones(observations.size, dtype=bool)
        if state_lower_bound is None
        else observations >= float(state_lower_bound)
    )
    return _minimum_replacements_from_matrix(
        compatibility, valid_reports
    )


def k_min(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
    *,
    state_lower_bound: float | None = None,
) -> int:
    return minimum_replacements_for_lipschitz(
        times,
        values,
        max_slope,
        state_lower_bound=state_lower_bound,
    ).minimum_replacements


def _safe_fraction_to_float(value: Fraction) -> float:
    try:
        converted = float(value)
    except OverflowError as error:
        raise OverflowError(
            "finite endpoint is outside the representable binary64 range"
        ) from error
    if value != 0 and converted == 0.0:
        raise FloatingPointError(
            "finite nonzero endpoint underflows the binary64 return type"
        )
    return converted


def _safe_scale_value(value: float, factor: Fraction) -> float:
    """Scale once in exact rational arithmetic, avoiding ``0 * inf``."""

    if math.isnan(value):
        raise FloatingPointError("endpoint quadrature produced NaN")
    if math.isinf(value):
        raise FloatingPointError(
            "conditioned finite-anchor quadrature produced an infinite endpoint"
        )
    if value == 0.0:
        return 0.0
    return _safe_fraction_to_float(_as_fraction(value) * factor)


def _scale_witness(
    witness: ExtremumWitness,
    factor: Fraction,
) -> ExtremumWitness:
    return replace(
        witness,
        value=_safe_scale_value(float(witness.value), factor),
    )


def _all_deleted_bounds(
    raw_times: np.ndarray,
    *,
    threshold: float,
    max_slope: float,
    budget: int,
    minimum: MinimumReplacementResult,
    state_lower_bound: float | None,
) -> RobustBurdenBounds:
    """Resolve the no-anchor model in original units before normalization."""

    horizon_exact = _as_fraction(raw_times[-1]) - _as_fraction(raw_times[0])
    horizon = _safe_fraction_to_float(horizon_exact)
    if not math.isfinite(horizon):
        raise OverflowError("time horizon is not representable as a finite float")
    deleted = tuple(range(int(raw_times.size)))
    zero_witness = ExtremumWitness(0.0, (), deleted, True)

    if state_lower_bound is None:
        occupation_upper = horizon
        deficit_upper = math.inf
        deficit_attained = False
        observation_mode = "exact_k_replacement_contamination"
    else:
        if threshold > state_lower_bound:
            occupation_upper = horizon
            deficit_upper = _safe_fraction_to_float(
                (_as_fraction(threshold) - _as_fraction(state_lower_bound))
                * horizon_exact
            )
        else:
            occupation_upper = 0.0
            deficit_upper = 0.0
        deficit_attained = True
        observation_mode = (
            "exact_k_replacement_contamination_with_state_floor="
            f"{state_lower_bound:.17g}"
        )

    occupation_upper_witness = ExtremumWitness(
        occupation_upper, (), deleted, True
    )
    deficit_upper_witness = ExtremumWitness(
        deficit_upper, (), deleted, deficit_attained
    )
    return RobustBurdenBounds(
        horizon=horizon,
        threshold=float(threshold),
        max_slope=float(max_slope),
        max_replacements=budget,
        minimum_replacements_for_feasibility=minimum.minimum_replacements,
        occupation_lower=0.0,
        occupation_upper=occupation_upper,
        oxygen_deficit_lower=0.0,
        oxygen_deficit_upper=deficit_upper,
        occupation_lower_witness=zero_witness,
        occupation_upper_witness=occupation_upper_witness,
        oxygen_deficit_lower_witness=zero_witness,
        oxygen_deficit_upper_witness=deficit_upper_witness,
        observation_mode=observation_mode,
        prototype_status=(
            _snapshot.PROTOTYPE_STATUS
            + "; exact_binary64_compatibility_current_runtime"
        ),
    )


def _zero_bounds_under_state_floor(
    raw_times: np.ndarray,
    *,
    threshold: float,
    max_slope: float,
    budget: int,
    minimum: MinimumReplacementResult,
    state_lower_bound: float,
) -> RobustBurdenBounds:
    """Return the analytic zero functional when ``threshold <= floor``."""

    horizon_exact = _as_fraction(raw_times[-1]) - _as_fraction(raw_times[0])
    horizon = _safe_fraction_to_float(horizon_exact)
    witness = ExtremumWitness(
        0.0,
        minimum.retained_indices,
        minimum.deleted_indices,
        True,
    )
    return RobustBurdenBounds(
        horizon=horizon,
        threshold=float(threshold),
        max_slope=float(max_slope),
        max_replacements=budget,
        minimum_replacements_for_feasibility=minimum.minimum_replacements,
        occupation_lower=0.0,
        occupation_upper=0.0,
        oxygen_deficit_lower=0.0,
        oxygen_deficit_upper=0.0,
        occupation_lower_witness=witness,
        occupation_upper_witness=witness,
        oxygen_deficit_lower_witness=witness,
        oxygen_deficit_upper_witness=witness,
        observation_mode=(
            "exact_k_replacement_contamination_with_state_floor="
            f"{state_lower_bound:.17g}"
        ),
        prototype_status=(
            _snapshot.PROTOTYPE_STATUS
            + "; analytic_threshold_at_or_below_state_floor"
        ),
    )


def identify_hypoxia_burden_with_replacements(
    times: Iterable[float],
    values: Iterable[float],
    *,
    threshold: float,
    max_slope: float,
    max_replacements: int,
    state_lower_bound: float | None = None,
) -> RobustBurdenBounds:
    """Return sharp four endpoints for the exact-inlier replacement model.

    If ``state_lower_bound`` is supplied, every feasible latent path must stay
    above that known physical floor.  For dissolved oxygen, use
    ``state_lower_bound=0.0``.
    """

    if not math.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if state_lower_bound is not None and not math.isfinite(state_lower_bound):
        raise ValueError("state_lower_bound must be finite")
    budget = _snapshot._validate_k(max_replacements)
    raw_times, observations = _validated_arrays(times, values, max_slope)
    number = int(raw_times.size)
    compatibility = _exact_compatibility_matrix(
        raw_times, observations, max_slope
    )
    valid_reports = (
        np.ones(number, dtype=bool)
        if state_lower_bound is None
        else observations >= float(state_lower_bound)
    )
    minimum = _minimum_replacements_from_matrix(
        compatibility, valid_reports
    )
    if budget < minimum.minimum_replacements:
        floor_text = (
            ""
            if state_lower_bound is None
            else f" under state_lower_bound={state_lower_bound}"
        )
        raise InsufficientReplacementBudget(
            f"max_replacements={budget} is infeasible{floor_text}; "
            f"k_min(L)={minimum.minimum_replacements}"
        )
    if budget >= number:
        return _all_deleted_bounds(
            raw_times,
            threshold=threshold,
            max_slope=max_slope,
            budget=budget,
            minimum=minimum,
            state_lower_bound=state_lower_bound,
        )
    if (
        state_lower_bound is not None
        and threshold <= state_lower_bound
    ):
        return _zero_bounds_under_state_floor(
            raw_times,
            threshold=threshold,
            max_slope=max_slope,
            budget=budget,
            minimum=minimum,
            state_lower_bound=state_lower_bound,
        )

    (
        tau,
        z,
        horizon,
        horizon_exact,
        value_scale_exact,
        _center_exact,
        normalized_slope,
        normalized_threshold,
        normalized_floor,
        exact_cost_context,
    ) = _normalized_problem(
        raw_times,
        observations,
        valid_reports,
        compatibility,
        max_slope,
        threshold=threshold,
        state_lower_bound=state_lower_bound,
    )
    occupation_lower = _optimize_subsequence(
        tau,
        z,
        compatibility,
        valid_reports,
        normalized_slope,
        normalized_threshold,
        normalized_floor,
        exact_cost_context,
        budget,
        envelope="upper",
        functional="occupation",
        maximize=False,
    )
    occupation_upper = _optimize_subsequence(
        tau,
        z,
        compatibility,
        valid_reports,
        normalized_slope,
        normalized_threshold,
        normalized_floor,
        exact_cost_context,
        budget,
        envelope="lower",
        functional="occupation",
        maximize=True,
    )
    deficit_lower = _optimize_subsequence(
        tau,
        z,
        compatibility,
        valid_reports,
        normalized_slope,
        normalized_threshold,
        normalized_floor,
        exact_cost_context,
        budget,
        envelope="upper",
        functional="deficit",
        maximize=False,
    )
    deficit_upper = _optimize_subsequence(
        tau,
        z,
        compatibility,
        valid_reports,
        normalized_slope,
        normalized_threshold,
        normalized_floor,
        exact_cost_context,
        budget,
        envelope="lower",
        functional="deficit",
        maximize=True,
    )
    observation_mode = (
        "exact_k_replacement_contamination"
        if state_lower_bound is None
        else (
            "exact_k_replacement_contamination_with_state_floor="
            f"{state_lower_bound:.17g}"
        )
    )
    normalized = RobustBurdenBounds(
        horizon=1.0,
        threshold=normalized_threshold,
        max_slope=normalized_slope,
        max_replacements=budget,
        minimum_replacements_for_feasibility=minimum.minimum_replacements,
        occupation_lower=occupation_lower.value,
        occupation_upper=occupation_upper.value,
        oxygen_deficit_lower=deficit_lower.value,
        oxygen_deficit_upper=deficit_upper.value,
        occupation_lower_witness=occupation_lower,
        occupation_upper_witness=occupation_upper,
        oxygen_deficit_lower_witness=deficit_lower,
        oxygen_deficit_upper_witness=deficit_upper,
        observation_mode=observation_mode,
        prototype_status=_snapshot.PROTOTYPE_STATUS,
    )
    occupation_factor = horizon_exact
    deficit_factor = horizon_exact * value_scale_exact
    return RobustBurdenBounds(
        horizon=horizon,
        threshold=float(threshold),
        max_slope=float(max_slope),
        max_replacements=normalized.max_replacements,
        minimum_replacements_for_feasibility=(
            normalized.minimum_replacements_for_feasibility
        ),
        occupation_lower=_safe_scale_value(
            normalized.occupation_lower, occupation_factor
        ),
        occupation_upper=_safe_scale_value(
            normalized.occupation_upper, occupation_factor
        ),
        oxygen_deficit_lower=_safe_scale_value(
            normalized.oxygen_deficit_lower, deficit_factor
        ),
        oxygen_deficit_upper=_safe_scale_value(
            normalized.oxygen_deficit_upper, deficit_factor
        ),
        occupation_lower_witness=_scale_witness(
            normalized.occupation_lower_witness, occupation_factor
        ),
        occupation_upper_witness=_scale_witness(
            normalized.occupation_upper_witness, occupation_factor
        ),
        oxygen_deficit_lower_witness=_scale_witness(
            normalized.oxygen_deficit_lower_witness, deficit_factor
        ),
        oxygen_deficit_upper_witness=_scale_witness(
            normalized.oxygen_deficit_upper_witness, deficit_factor
        ),
        observation_mode=normalized.observation_mode,
        prototype_status=normalized.prototype_status
        + "; exact_binary64_compatibility_conditioned_current_runtime",
    )


__all__ = [
    "NUMERICAL_POLICY",
    "NORMALIZATION_CONDITION_LIMIT",
    "CRITICAL_GEOMETRY_FALLBACK_LIMIT",
    "InsufficientReplacementBudget",
    "ExtremumWitness",
    "MinimumReplacementResult",
    "RobustBurdenBounds",
    "minimum_replacements_for_lipschitz",
    "k_min",
    "identify_hypoxia_burden_with_replacements",
]
