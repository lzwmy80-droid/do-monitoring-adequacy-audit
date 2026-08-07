# -*- coding: utf-8 -*-
"""Small-n oracle for a nonnegative Lipschitz state constraint.

This file is deliberately independent of the production dynamic program.  It
enumerates every retained observation subset and is intended only for theorem
design and regression tests.  Complexity is exponential in the number of
reports, so inputs are capped at ``MAX_ENUMERATION_POINTS``.

Model
-----
The latent path is globally ``L``-Lipschitz, satisfies ``x(t) >= 0``, and
matches every retained exact report.  At most ``k`` reports may be replaced.
For each feasible retained subset, the unconstrained lower envelope is clipped
at zero and the upper envelope is unchanged.  Both constrained envelopes are
then integrated exactly.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Literal, Sequence


MAX_ENUMERATION_POINTS = 18


class NonnegativePrototypeInfeasible(ValueError):
    """No nonnegative retained subset fits the replacement budget."""


@dataclass(frozen=True)
class Segment:
    start: float
    end: float
    slope: float
    intercept: float

    def value(self, time: float) -> float:
        return self.slope * time + self.intercept


@dataclass(frozen=True)
class PrototypeWitness:
    value: float
    retained_indices: tuple[int, ...]
    deleted_indices: tuple[int, ...]
    path_rule: str
    attained: bool = True


@dataclass(frozen=True)
class NonnegativePrototypeBounds:
    horizon: float
    threshold: float
    max_slope: float
    max_replacements: int
    minimum_replacements_for_feasibility: int
    occupation_lower: float
    occupation_upper: float
    oxygen_deficit_lower: float
    oxygen_deficit_upper: float
    occupation_lower_witness: PrototypeWitness
    occupation_upper_witness: PrototypeWitness
    oxygen_deficit_lower_witness: PrototypeWitness
    oxygen_deficit_upper_witness: PrototypeWitness


Envelope = Literal["lower", "upper"]


def _validated_inputs(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
    max_replacements: int,
) -> tuple[tuple[float, ...], tuple[float, ...], int, float]:
    raw_times = tuple(float(value) for value in times)
    reports = tuple(float(value) for value in values)
    if len(raw_times) < 2:
        raise ValueError("at least two times are required")
    if len(raw_times) != len(reports):
        raise ValueError("times and values must have equal lengths")
    if len(raw_times) > MAX_ENUMERATION_POINTS:
        raise ValueError(
            f"prototype enumeration is capped at {MAX_ENUMERATION_POINTS} points"
        )
    if not all(math.isfinite(value) for value in raw_times + reports):
        raise ValueError("times and values must be finite")
    if any(right <= left for left, right in zip(raw_times[:-1], raw_times[1:])):
        raise ValueError("times must be strictly increasing")
    if not math.isfinite(max_slope) or max_slope < 0:
        raise ValueError("max_slope must be finite and non-negative")
    if isinstance(max_replacements, bool) or not isinstance(max_replacements, int):
        raise TypeError("max_replacements must be a non-negative integer")
    if max_replacements < 0:
        raise ValueError("max_replacements must be a non-negative integer")

    normalized = tuple(value - raw_times[0] for value in raw_times)
    return normalized, reports, max_replacements, 0.0


def _compatible(
    times: Sequence[float],
    values: Sequence[float],
    retained: Sequence[int],
    max_slope: float,
    _unused_tolerance: float,
) -> bool:
    # The physical state constraint is exact: a negative retained report is
    # inadmissible, rather than silently projected to zero.
    if any(values[index] < 0.0 for index in retained):
        return False
    slope = Fraction.from_float(max_slope)
    return all(
        abs(
            Fraction.from_float(values[right])
            - Fraction.from_float(values[left])
        )
        <= slope
        * (
            Fraction.from_float(times[right])
            - Fraction.from_float(times[left])
        )
        for left, right in zip(retained[:-1], retained[1:])
    )


def _selected_two_lines(
    start: float,
    end: float,
    first: tuple[float, float],
    second: tuple[float, float],
    *,
    choose_maximum: bool,
) -> list[Segment]:
    slope_a, intercept_a = first
    slope_b, intercept_b = second
    points = [start, end]
    slope_difference = slope_a - slope_b
    if slope_difference != 0.0:
        crossing = (intercept_b - intercept_a) / slope_difference
        if start < crossing < end:
            points.insert(1, crossing)

    output: list[Segment] = []
    for left, right in zip(points[:-1], points[1:]):
        midpoint = 0.5 * (left + right)
        first_value = slope_a * midpoint + intercept_a
        second_value = slope_b * midpoint + intercept_b
        use_first = (
            first_value >= second_value
            if choose_maximum
            else first_value <= second_value
        )
        slope, intercept = first if use_first else second
        output.append(Segment(left, right, slope, intercept))
    return output


def _boundary_segments(
    times: Sequence[float],
    values: Sequence[float],
    index: int,
    max_slope: float,
    envelope: Envelope,
    side: Literal["left", "right"],
) -> list[Segment]:
    anchor_time = times[index]
    anchor_value = values[index]
    if side == "left":
        start, end = 0.0, anchor_time
        slope = max_slope if envelope == "lower" else -max_slope
    else:
        start, end = anchor_time, times[-1]
        slope = -max_slope if envelope == "lower" else max_slope
    if end <= start:
        return []
    return [Segment(start, end, slope, anchor_value - slope * anchor_time)]


def _pair_segments(
    times: Sequence[float],
    values: Sequence[float],
    left: int,
    right: int,
    max_slope: float,
    envelope: Envelope,
) -> list[Segment]:
    left_time, right_time = times[left], times[right]
    left_value, right_value = values[left], values[right]
    if envelope == "lower":
        first = (-max_slope, left_value + max_slope * left_time)
        second = (max_slope, right_value - max_slope * right_time)
        choose_maximum = True
    else:
        first = (max_slope, left_value - max_slope * left_time)
        second = (-max_slope, right_value + max_slope * right_time)
        choose_maximum = False
    return _selected_two_lines(
        left_time,
        right_time,
        first,
        second,
        choose_maximum=choose_maximum,
    )


def _envelope_segments(
    times: Sequence[float],
    values: Sequence[float],
    retained: Sequence[int],
    max_slope: float,
    envelope: Envelope,
) -> list[Segment]:
    output = _boundary_segments(
        times, values, retained[0], max_slope, envelope, "left"
    )
    for left, right in zip(retained[:-1], retained[1:]):
        output.extend(
            _pair_segments(times, values, left, right, max_slope, envelope)
        )
    output.extend(
        _boundary_segments(
            times, values, retained[-1], max_slope, envelope, "right"
        )
    )
    return output


def clip_segments_to_floor(
    segments: Sequence[Segment],
    floor: float = 0.0,
) -> list[Segment]:
    """Return exact affine pieces for ``max(path, floor)``."""

    if not math.isfinite(floor):
        raise ValueError("floor must be finite")
    output: list[Segment] = []
    for segment in segments:
        points = [segment.start, segment.end]
        if segment.slope != 0.0:
            crossing = (floor - segment.intercept) / segment.slope
            if segment.start < crossing < segment.end:
                points.insert(1, crossing)
        for left, right in zip(points[:-1], points[1:]):
            midpoint = 0.5 * (left + right)
            if segment.value(midpoint) < floor:
                output.append(Segment(left, right, 0.0, floor))
            else:
                output.append(
                    Segment(
                        left,
                        right,
                        segment.slope,
                        segment.intercept,
                    )
                )
    return output


def _below_interval(
    segment: Segment,
    threshold: float,
) -> tuple[float, float] | None:
    if segment.slope == 0.0:
        midpoint = 0.5 * (segment.start + segment.end)
        return (
            (segment.start, segment.end)
            if segment.value(midpoint) < threshold
            else None
        )
    crossing = (threshold - segment.intercept) / segment.slope
    if segment.slope > 0:
        end = min(segment.end, crossing)
        return (segment.start, end) if end > segment.start else None
    start = max(segment.start, crossing)
    return (start, segment.end) if segment.end > start else None


def integrate_occupation(
    segments: Sequence[Segment],
    threshold: float,
) -> float:
    return float(
        math.fsum(
            interval[1] - interval[0]
            for segment in segments
            if (interval := _below_interval(segment, threshold)) is not None
        )
    )


def integrate_deficit(
    segments: Sequence[Segment],
    threshold: float,
) -> float:
    contributions = []
    for segment in segments:
        interval = _below_interval(segment, threshold)
        if interval is None:
            continue
        start, end = interval
        value_start = segment.value(start)
        value_end = segment.value(end)
        contributions.append(
            max(
                0.0,
                0.5
                * (
                    (threshold - value_start)
                    + (threshold - value_end)
                )
                * (end - start),
            )
        )
    return float(max(0.0, math.fsum(contributions)))


def minimum_replacements_for_nonnegative_lipschitz(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
) -> int:
    normalized, reports, _, tolerance = _validated_inputs(
        times, values, max_slope, 0
    )
    number = len(normalized)
    for retained_count in range(number, 0, -1):
        for retained in itertools.combinations(range(number), retained_count):
            if _compatible(
                normalized,
                reports,
                retained,
                max_slope,
                tolerance,
            ):
                return number - retained_count
    return number


def _witness(
    value: float,
    retained: tuple[int, ...],
    number: int,
    path_rule: str,
) -> PrototypeWitness:
    retained_set = set(retained)
    deleted = tuple(index for index in range(number) if index not in retained_set)
    return PrototypeWitness(
        value=float(value),
        retained_indices=retained,
        deleted_indices=deleted,
        path_rule=path_rule,
    )


def enumerate_nonnegative_bounds(
    times: Iterable[float],
    values: Iterable[float],
    *,
    threshold: float,
    max_slope: float,
    max_replacements: int,
) -> NonnegativePrototypeBounds:
    """Enumerate sharp four endpoints under ``x(t) >= 0`` for small ``n``."""

    if not math.isfinite(threshold):
        raise ValueError("threshold must be finite")
    t, y, budget, tolerance = _validated_inputs(
        times, values, max_slope, max_replacements
    )
    number = len(t)
    horizon = t[-1]
    minimum = minimum_replacements_for_nonnegative_lipschitz(
        t, y, max_slope
    )
    if budget < minimum:
        raise NonnegativePrototypeInfeasible(
            f"max_replacements={budget} is infeasible; "
            f"nonnegative k_min(L)={minimum}"
        )

    candidates: list[
        tuple[
            tuple[int, ...],
            float,
            float,
            float,
            float,
        ]
    ] = []
    minimum_retained = max(1, number - budget)
    for retained_count in range(minimum_retained, number + 1):
        for retained in itertools.combinations(range(number), retained_count):
            if not _compatible(t, y, retained, max_slope, tolerance):
                continue
            lower = clip_segments_to_floor(
                _envelope_segments(
                    t, y, retained, max_slope, envelope="lower"
                )
            )
            upper = clip_segments_to_floor(
                _envelope_segments(
                    t, y, retained, max_slope, envelope="upper"
                )
            )
            candidates.append(
                (
                    retained,
                    integrate_occupation(upper, threshold),
                    integrate_occupation(lower, threshold),
                    integrate_deficit(upper, threshold),
                    integrate_deficit(lower, threshold),
                )
            )

    # The empty retained set is feasible only when every report may be
    # replaced.  Nonnegativity makes its deficit supremum finite and attained.
    if budget >= number:
        upper_occ = horizon if threshold > 0.0 else 0.0
        upper_deficit = max(threshold, 0.0) * horizon
        candidates.append(((), 0.0, upper_occ, 0.0, upper_deficit))

    if not candidates:
        raise NonnegativePrototypeInfeasible(
            "no nonnegative retained subset fits the replacement budget"
        )

    occupation_lower_row = min(candidates, key=lambda row: row[1])
    occupation_upper_row = max(candidates, key=lambda row: row[2])
    deficit_lower_row = min(candidates, key=lambda row: row[3])
    deficit_upper_row = max(candidates, key=lambda row: row[4])

    return NonnegativePrototypeBounds(
        horizon=horizon,
        threshold=threshold,
        max_slope=max_slope,
        max_replacements=budget,
        minimum_replacements_for_feasibility=minimum,
        occupation_lower=occupation_lower_row[1],
        occupation_upper=occupation_upper_row[2],
        oxygen_deficit_lower=deficit_lower_row[3],
        oxygen_deficit_upper=deficit_upper_row[4],
        occupation_lower_witness=_witness(
            occupation_lower_row[1],
            occupation_lower_row[0],
            number,
            "upper_envelope",
        ),
        occupation_upper_witness=_witness(
            occupation_upper_row[2],
            occupation_upper_row[0],
            number,
            "clip(lower_envelope, 0)",
        ),
        oxygen_deficit_lower_witness=_witness(
            deficit_lower_row[3],
            deficit_lower_row[0],
            number,
            "upper_envelope",
        ),
        oxygen_deficit_upper_witness=_witness(
            deficit_upper_row[4],
            deficit_upper_row[0],
            number,
            "clip(lower_envelope, 0)",
        ),
    )


__all__ = [
    "MAX_ENUMERATION_POINTS",
    "NonnegativePrototypeBounds",
    "NonnegativePrototypeInfeasible",
    "PrototypeWitness",
    "Segment",
    "clip_segments_to_floor",
    "enumerate_nonnegative_bounds",
    "integrate_deficit",
    "integrate_occupation",
    "minimum_replacements_for_nonnegative_lipschitz",
]
