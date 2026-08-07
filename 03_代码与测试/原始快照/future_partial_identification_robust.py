# -*- coding: utf-8 -*-
"""Prototype: sharp burden bounds with up to ``k`` replaced observations.

The observation times and the analysis domain ``[t[0], t[-1]]`` are fixed.
At most ``k`` reported values may be arbitrary replacements (gross errors).
Every non-replaced value must lie exactly on one globally ``L``-Lipschitz
trajectory.  For this finite-sample model, this module computes the sharp
identified lower and upper bounds for threshold occupation and cumulative
threshold deficit.

For a retained feasible subsequence, the McShane--Whitney lower and upper
envelopes depend only on the first retained point, adjacent retained pairs,
and the last retained point.  Their integrals are therefore additive.  Four
shortest/longest-path dynamic programs find the extremizing subsequences in
``O(k n^2)`` time and ``O(k n)`` memory.

This is an algorithmically exact research prototype, not a formally published
theorem or peer-reviewed implementation.  Lipschitz-extension algorithms have
substantial precedent; see Kyng, Rao, Sachdeva, and Spielman (2015),
"Algorithms for Lipschitz Learning on Graphs", COLT/PMLR 40.  The ordered
maximum-feasible-subsequence calculation used for ``k_min`` is a standard DAG
longest-path dynamic program and is not claimed as a new algorithm.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Literal, Sequence

import numpy as np

from future_partial_identification import (
    AffineSegment,
    integrate_occupation,
    integrate_oxygen_deficit,
)


PROTOTYPE_STATUS = (
    "algorithmically_exact_finite_sample_prototype; "
    "theorem_not_formally_peer_reviewed"
)


class InsufficientReplacementBudget(ValueError):
    """Raised when no nonempty feasible retained set exists within ``k``."""


@dataclass(frozen=True)
class ExtremumWitness:
    """One retained/deleted-index witness for a sharp functional endpoint."""

    value: float
    retained_indices: tuple[int, ...]
    deleted_indices: tuple[int, ...]
    attained: bool


@dataclass(frozen=True)
class MinimumReplacementResult:
    """Minimum replacements needed for exact ``L``-Lipschitz consistency."""

    minimum_replacements: int
    retained_indices: tuple[int, ...]
    deleted_indices: tuple[int, ...]


@dataclass(frozen=True)
class RobustBurdenBounds:
    """Sharp extended-real identified interval under replacement contamination."""

    horizon: float
    threshold: float
    max_slope: float
    max_replacements: int
    minimum_replacements_for_feasibility: int
    occupation_lower: float
    occupation_upper: float
    oxygen_deficit_lower: float
    oxygen_deficit_upper: float
    occupation_lower_witness: ExtremumWitness
    occupation_upper_witness: ExtremumWitness
    oxygen_deficit_lower_witness: ExtremumWitness
    oxygen_deficit_upper_witness: ExtremumWitness
    observation_mode: str
    prototype_status: str


Envelope = Literal["lower", "upper"]
Functional = Literal["occupation", "deficit"]


def _as_inputs(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
) -> tuple[np.ndarray, np.ndarray, float]:
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
    normalized = raw_times - raw_times[0]
    scale = max(
        1.0,
        float(np.max(np.abs(observations))),
        float(max_slope * normalized[-1]),
    )
    return normalized, observations, 1e-11 * scale


def _validate_k(max_replacements: int) -> int:
    if isinstance(max_replacements, bool) or not isinstance(
        max_replacements, (int, np.integer)
    ):
        raise TypeError("max_replacements must be a non-negative integer")
    if max_replacements < 0:
        raise ValueError("max_replacements must be a non-negative integer")
    return int(max_replacements)


def _compatible(
    times: np.ndarray,
    values: np.ndarray,
    left: int,
    right: int,
    max_slope: float,
    tolerance: float,
) -> bool:
    return bool(
        abs(float(values[right] - values[left]))
        <= max_slope * float(times[right] - times[left]) + tolerance
    )


def minimum_replacements_for_lipschitz(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
) -> MinimumReplacementResult:
    """Return ``k_min(L)`` and a maximum-cardinality feasible subsequence.

    Pairwise consistency of adjacent retained points is sufficient on the
    ordered line.  Thus ``k_min = n - longest_path_length`` in the DAG whose
    edge ``i -> j`` exists when ``|y[j]-y[i]| <= L (t[j]-t[i])``.

    This combines standard ordered-subsequence dynamic programming with the
    well-established Lipschitz-extension feasibility primitive; it is not a
    novelty claim (cf. Kyng et al., 2015).
    """

    t, y, tolerance = _as_inputs(times, values, max_slope)
    number = int(t.size)
    lengths = [1] * number
    parents: list[int | None] = [None] * number
    for right in range(number):
        for left in range(right):
            if not _compatible(t, y, left, right, max_slope, tolerance):
                continue
            candidate = lengths[left] + 1
            if candidate > lengths[right]:
                lengths[right] = candidate
                parents[right] = left

    best_end = max(range(number), key=lambda index: lengths[index])
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


def k_min(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
) -> int:
    """Return only the scalar ``k_min(L)``; see the witness-returning function."""

    return minimum_replacements_for_lipschitz(
        times,
        values,
        max_slope,
    ).minimum_replacements


def _selected_two_line_segments(
    start: float,
    end: float,
    first: tuple[float, float],
    second: tuple[float, float],
    *,
    choose_maximum: bool,
) -> list[AffineSegment]:
    """Return the max/min of two affine lines on one closed interval."""

    slope_a, intercept_a = first
    slope_b, intercept_b = second
    points = [start, end]
    slope_difference = slope_a - slope_b
    if abs(slope_difference) > 1e-15:
        crossing = (intercept_b - intercept_a) / slope_difference
        tolerance = 1e-12 * max(1.0, abs(start), abs(end))
        if start + tolerance < crossing < end - tolerance:
            points.insert(1, crossing)

    output: list[AffineSegment] = []
    for left, right in zip(points[:-1], points[1:]):
        midpoint = 0.5 * (left + right)
        value_a = slope_a * midpoint + intercept_a
        value_b = slope_b * midpoint + intercept_b
        use_a = value_a >= value_b if choose_maximum else value_a <= value_b
        slope, intercept = first if use_a else second
        output.append(AffineSegment(left, right, slope, intercept))
    return output


def _boundary_segments(
    times: np.ndarray,
    values: np.ndarray,
    index: int,
    max_slope: float,
    envelope: Envelope,
    side: Literal["left", "right"],
) -> list[AffineSegment]:
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
    return [AffineSegment(start, end, slope, intercept)]


def _pair_segments(
    times: np.ndarray,
    values: np.ndarray,
    left: int,
    right: int,
    max_slope: float,
    envelope: Envelope,
) -> list[AffineSegment]:
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


def _integral_cost(
    segments: Sequence[AffineSegment],
    threshold: float,
    functional: Functional,
) -> float:
    if functional == "occupation":
        return integrate_occupation(segments, threshold)
    return integrate_oxygen_deficit(segments, threshold)


def _retained_subset_value(
    times: np.ndarray,
    values: np.ndarray,
    retained: Sequence[int],
    max_slope: float,
    threshold: float,
    envelope: Envelope,
    functional: Functional,
) -> float:
    """Exact additive integral for one nonempty feasible retained subset."""

    if not retained:
        raise ValueError("retained must be nonempty")
    segments = _boundary_segments(
        times, values, retained[0], max_slope, envelope, "left"
    )
    for left, right in zip(retained[:-1], retained[1:]):
        segments.extend(
            _pair_segments(times, values, left, right, max_slope, envelope)
        )
    segments.extend(
        _boundary_segments(
            times, values, retained[-1], max_slope, envelope, "right"
        )
    )
    return _integral_cost(segments, threshold, functional)


def _is_better(candidate: float, incumbent: float, maximize: bool) -> bool:
    scale = max(1.0, abs(candidate), abs(incumbent))
    tolerance = 1e-12 * scale
    return candidate > incumbent + tolerance if maximize else candidate < incumbent - tolerance


def _optimize_subsequence(
    times: np.ndarray,
    values: np.ndarray,
    max_slope: float,
    threshold: float,
    max_replacements: int,
    tolerance: float,
    *,
    envelope: Envelope,
    functional: Functional,
    maximize: bool,
) -> ExtremumWitness:
    """Optimize one envelope integral by an ``O(k n^2)`` subsequence DP."""

    number = int(times.size)
    budget = min(max_replacements, number - 1)
    scores: dict[tuple[int, int], float] = {}
    parents: dict[tuple[int, int], tuple[int, int] | None] = {}

    for right in range(number):
        if right <= budget:
            state = (right, right)
            scores[state] = _integral_cost(
                _boundary_segments(
                    times, values, right, max_slope, envelope, "left"
                ),
                threshold,
                functional,
            )
            parents[state] = None

        for left in range(right):
            if not _compatible(
                times, values, left, right, max_slope, tolerance
            ):
                continue
            skipped = right - left - 1
            edge_cost = _integral_cost(
                _pair_segments(
                    times, values, left, right, max_slope, envelope
                ),
                threshold,
                functional,
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
        total_deleted = deleted_before + (number - last - 1)
        if total_deleted > budget:
            continue
        value = prefix_value + _integral_cost(
            _boundary_segments(
                times, values, last, max_slope, envelope, "right"
            ),
            threshold,
            functional,
        )
        if best_state is None or _is_better(value, best_value, maximize):
            best_value = value
            best_state = state

    if best_state is None:
        raise InsufficientReplacementBudget(
            "no nonempty L-Lipschitz retained subset fits the replacement budget"
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


def identify_hypoxia_burden_with_replacements(
    times: Iterable[float],
    values: Iterable[float],
    *,
    threshold: float,
    max_slope: float,
    max_replacements: int,
) -> RobustBurdenBounds:
    """Compute sharp bounds with at most ``max_replacements`` gross errors.

    Separate retained/deleted witnesses are returned for all four endpoints;
    the same contaminated-index set need not attain every endpoint.

    If every observation may be replaced, occupation remains sharply bounded
    by ``[0, T]`` and deficit by ``[0, +inf)``.  The infinite upper deficit is
    a supremum and has no finite attaining trajectory.
    """

    if not math.isfinite(threshold):
        raise ValueError("threshold must be finite")
    budget = _validate_k(max_replacements)
    t, y, tolerance = _as_inputs(times, values, max_slope)
    number = int(t.size)
    horizon = float(t[-1])
    minimum = minimum_replacements_for_lipschitz(t, y, max_slope)

    if budget >= number:
        deleted = tuple(range(number))
        no_anchor_low = ExtremumWitness(0.0, (), deleted, True)
        no_anchor_high_occupation = ExtremumWitness(horizon, (), deleted, True)
        no_anchor_high_deficit = ExtremumWitness(math.inf, (), deleted, False)
        return RobustBurdenBounds(
            horizon=horizon,
            threshold=float(threshold),
            max_slope=float(max_slope),
            max_replacements=budget,
            minimum_replacements_for_feasibility=minimum.minimum_replacements,
            occupation_lower=0.0,
            occupation_upper=horizon,
            oxygen_deficit_lower=0.0,
            oxygen_deficit_upper=math.inf,
            occupation_lower_witness=no_anchor_low,
            occupation_upper_witness=no_anchor_high_occupation,
            oxygen_deficit_lower_witness=no_anchor_low,
            oxygen_deficit_upper_witness=no_anchor_high_deficit,
            observation_mode="exact_k_replacement_contamination",
            prototype_status=PROTOTYPE_STATUS,
        )

    if budget < minimum.minimum_replacements:
        raise InsufficientReplacementBudget(
            f"max_replacements={budget} is infeasible; "
            f"k_min(L)={minimum.minimum_replacements}"
        )

    occupation_lower = _optimize_subsequence(
        t,
        y,
        max_slope,
        threshold,
        budget,
        tolerance,
        envelope="upper",
        functional="occupation",
        maximize=False,
    )
    occupation_upper = _optimize_subsequence(
        t,
        y,
        max_slope,
        threshold,
        budget,
        tolerance,
        envelope="lower",
        functional="occupation",
        maximize=True,
    )
    deficit_lower = _optimize_subsequence(
        t,
        y,
        max_slope,
        threshold,
        budget,
        tolerance,
        envelope="upper",
        functional="deficit",
        maximize=False,
    )
    deficit_upper = _optimize_subsequence(
        t,
        y,
        max_slope,
        threshold,
        budget,
        tolerance,
        envelope="lower",
        functional="deficit",
        maximize=True,
    )
    return RobustBurdenBounds(
        horizon=horizon,
        threshold=float(threshold),
        max_slope=float(max_slope),
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
        observation_mode="exact_k_replacement_contamination",
        prototype_status=PROTOTYPE_STATUS,
    )

