# -*- coding: utf-8 -*-
"""Finite-sample partial identification of hypoxia burden.

This module treats the unobserved dissolved-oxygen trajectory as an unknown
globally L-Lipschitz function.  Observations may be exact or interval-valued.
It constructs the pointwise smallest and largest feasible trajectories and
integrates them exactly to obtain sharp bounds for total time below a threshold
and cumulative oxygen deficit.

Event counts are deliberately not estimated: a bounded slope alone does not
give a finite upper bound on the number of threshold crossings.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    PROJECT
    / "research_branch_outputs"
    / "future_agenda_20260716"
    / "partial_identification_prototype"
)


class InfeasibleObservationIntervals(ValueError):
    """Raised when no L-Lipschitz path can meet all observation intervals."""


@dataclass(frozen=True)
class AffineSegment:
    """An affine envelope piece ``value(t) = slope * t + intercept``."""

    start: float
    end: float
    slope: float
    intercept: float

    def value(self, time: float) -> float:
        return self.slope * time + self.intercept


@dataclass(frozen=True)
class BurdenBounds:
    """Sharp identified interval for two monotone trajectory functionals."""

    horizon: float
    threshold: float
    max_slope: float
    occupation_lower: float
    occupation_upper: float
    occupation_fraction_lower: float
    occupation_fraction_upper: float
    oxygen_deficit_lower: float
    oxygen_deficit_upper: float
    event_count_status: str
    observation_mode: str


def deficit_width_bound_exact_observations(
    times: Iterable[float],
    values: Iterable[float],
    max_slope: float,
) -> float:
    """Upper-bound deficit identification width for exact observations.

    The returned irregular-grid bound is
    ``sum((L^2 Delta_i^2 - |Delta y_i|^2) / (2L))``.  It is the
    integrated feasible-corridor width; the actual deficit width can be smaller
    because the hinge functional may be flat above the threshold.
    """

    exact = list(values)
    normalized_times, _, _ = _prepare_inputs(times, exact, exact, max_slope)
    if max_slope == 0:
        return 0.0
    delta_t = np.diff(normalized_times)
    delta_y = np.diff(np.asarray(exact, dtype=float))
    terms = (
        max_slope * max_slope * delta_t * delta_t - delta_y * delta_y
    ) / (2.0 * max_slope)
    return float(np.sum(np.maximum(0.0, terms)))


def deficit_width_bound_uniform_error(
    times: Iterable[float],
    max_slope: float,
    error_halfwidth: float,
) -> float:
    """Uniform bound with deterministic symmetric endpoint error ``+/-epsilon``.

    This uses ``2 epsilon T + (L/2) sum Delta_i^2`` and does not require
    observed values.  It is a design bound, not an empirical guarantee that the
    chosen L or error half-width is scientifically valid.
    """

    time_vector = _as_finite_vector(times, "times")
    if time_vector.size < 2 or np.any(np.diff(time_vector) <= 0):
        raise ValueError("times must contain at least two strictly increasing values")
    if not math.isfinite(max_slope) or max_slope < 0:
        raise ValueError("max_slope must be finite and non-negative")
    if not math.isfinite(error_halfwidth) or error_halfwidth < 0:
        raise ValueError("error_halfwidth must be finite and non-negative")
    delta = np.diff(time_vector)
    horizon = float(time_vector[-1] - time_vector[0])
    return float(
        2.0 * error_halfwidth * horizon
        + 0.5 * max_slope * np.sum(delta * delta)
    )


def occupation_width_margin_bound(
    times: Iterable[float],
    max_slope: float,
    error_halfwidth: float,
    margin_constant: float,
    margin_exponent: float,
) -> float:
    """Return ``min(T, C (L h + 2 epsilon)^alpha)`` under a margin condition."""

    time_vector = _as_finite_vector(times, "times")
    if time_vector.size < 2 or np.any(np.diff(time_vector) <= 0):
        raise ValueError("times must contain at least two strictly increasing values")
    parameters = {
        "max_slope": max_slope,
        "error_halfwidth": error_halfwidth,
        "margin_constant": margin_constant,
        "margin_exponent": margin_exponent,
    }
    if any(not math.isfinite(value) or value < 0 for value in parameters.values()):
        raise ValueError("margin-bound parameters must be finite and non-negative")
    horizon = float(time_vector[-1] - time_vector[0])
    max_gap = float(np.max(np.diff(time_vector)))
    radius = max_slope * max_gap + 2.0 * error_halfwidth
    return float(min(horizon, margin_constant * radius**margin_exponent))


def _as_finite_vector(values: Iterable[float], name: str) -> np.ndarray:
    vector = np.asarray(list(values), dtype=float)
    if vector.ndim != 1 or vector.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} contains non-finite values")
    return vector


def _prepare_inputs(
    times: Iterable[float],
    observation_lower: Iterable[float],
    observation_upper: Iterable[float],
    max_slope: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw_times = _as_finite_vector(times, "times")
    lower = _as_finite_vector(observation_lower, "observation_lower")
    upper = _as_finite_vector(observation_upper, "observation_upper")
    if raw_times.size < 2:
        raise ValueError("at least two observation times are required")
    if lower.size != raw_times.size or upper.size != raw_times.size:
        raise ValueError("times and observation intervals must have equal lengths")
    if np.any(np.diff(raw_times) <= 0):
        raise ValueError("times must be strictly increasing")
    if np.any(lower > upper):
        raise ValueError("every observation lower bound must be <= its upper bound")
    if not math.isfinite(max_slope) or max_slope < 0:
        raise ValueError("max_slope must be finite and non-negative")

    normalized_times = raw_times - raw_times[0]
    distance = np.abs(normalized_times[:, None] - normalized_times[None, :])
    required_lower_at_observations = np.max(
        lower[None, :] - max_slope * distance,
        axis=1,
    )
    scale = max(1.0, float(np.max(np.abs(np.concatenate((lower, upper))))))
    tolerance = 1e-11 * scale
    violations = np.flatnonzero(required_lower_at_observations > upper + tolerance)
    if violations.size:
        first = int(violations[0])
        raise InfeasibleObservationIntervals(
            "no globally L-Lipschitz path satisfies all intervals; "
            f"at observation {first}, required lower value "
            f"{required_lower_at_observations[first]:.12g} exceeds upper "
            f"{upper[first]:.12g}"
        )
    return normalized_times, lower, upper


def _append_selected_segments(
    output: list[AffineSegment],
    start: float,
    end: float,
    first: tuple[float, float],
    second: tuple[float, float],
    choose_maximum: bool,
) -> None:
    """Append the max/min of two affine lines over one observation interval."""

    slope_a, intercept_a = first
    slope_b, intercept_b = second
    split_points = [start, end]
    slope_difference = slope_a - slope_b
    if abs(slope_difference) > 1e-15:
        crossing = (intercept_b - intercept_a) / slope_difference
        tolerance = 1e-12 * max(1.0, abs(start), abs(end))
        if start + tolerance < crossing < end - tolerance:
            split_points.insert(1, crossing)

    for left, right in zip(split_points[:-1], split_points[1:]):
        midpoint = (left + right) / 2.0
        value_a = slope_a * midpoint + intercept_a
        value_b = slope_b * midpoint + intercept_b
        select_a = value_a >= value_b if choose_maximum else value_a <= value_b
        slope, intercept = first if select_a else second
        output.append(AffineSegment(left, right, slope, intercept))


def construct_envelopes(
    times: Iterable[float],
    observation_lower: Iterable[float],
    observation_upper: Iterable[float],
    max_slope: float,
) -> tuple[list[AffineSegment], list[AffineSegment]]:
    """Return sharp pointwise lower and upper feasible-path envelopes.

    The feasible class is all globally ``max_slope``-Lipschitz functions whose
    value at every observation time lies in the corresponding closed interval.
    Segment times are normalized so the first observation is at zero.
    """

    t, lower, upper = _prepare_inputs(
        times,
        observation_lower,
        observation_upper,
        max_slope,
    )
    if max_slope == 0:
        return (
            [AffineSegment(0.0, float(t[-1]), 0.0, float(np.max(lower)))],
            [AffineSegment(0.0, float(t[-1]), 0.0, float(np.min(upper)))],
        )

    lower_left = np.maximum.accumulate(lower + max_slope * t)
    lower_right = np.maximum.accumulate((lower - max_slope * t)[::-1])[::-1]
    upper_left = np.minimum.accumulate(upper - max_slope * t)
    upper_right = np.minimum.accumulate((upper + max_slope * t)[::-1])[::-1]

    lower_segments: list[AffineSegment] = []
    upper_segments: list[AffineSegment] = []
    for index in range(t.size - 1):
        start = float(t[index])
        end = float(t[index + 1])
        _append_selected_segments(
            lower_segments,
            start,
            end,
            (-max_slope, float(lower_left[index])),
            (max_slope, float(lower_right[index + 1])),
            choose_maximum=True,
        )
        _append_selected_segments(
            upper_segments,
            start,
            end,
            (max_slope, float(upper_left[index])),
            (-max_slope, float(upper_right[index + 1])),
            choose_maximum=False,
        )
    return lower_segments, upper_segments


def _below_subinterval(
    segment: AffineSegment,
    threshold: float,
) -> tuple[float, float] | None:
    """Return the subinterval where an affine segment is strictly below H."""

    start, end = segment.start, segment.end
    slope, intercept = segment.slope, segment.intercept
    if abs(slope) <= 1e-15:
        return (start, end) if segment.value((start + end) / 2.0) < threshold else None

    crossing = (threshold - intercept) / slope
    if slope > 0:
        below_end = min(end, crossing)
        return (start, below_end) if below_end > start else None
    below_start = max(start, crossing)
    return (below_start, end) if end > below_start else None


def integrate_occupation(
    segments: Sequence[AffineSegment],
    threshold: float,
) -> float:
    """Exactly integrate time for which a piecewise-affine path is below H."""

    total = 0.0
    for segment in segments:
        subinterval = _below_subinterval(segment, threshold)
        if subinterval is not None:
            total += subinterval[1] - subinterval[0]
    return float(total)


def integrate_oxygen_deficit(
    segments: Sequence[AffineSegment],
    threshold: float,
) -> float:
    """Exactly integrate ``max(threshold - path(t), 0)``."""

    total = 0.0
    for segment in segments:
        subinterval = _below_subinterval(segment, threshold)
        if subinterval is None:
            continue
        start, end = subinterval
        total += (
            (threshold - segment.intercept) * (end - start)
            - 0.5 * segment.slope * (end * end - start * start)
        )
    return float(max(0.0, total))


def identify_hypoxia_burden(
    times: Iterable[float],
    values: Iterable[float] | None = None,
    *,
    observation_lower: Iterable[float] | None = None,
    observation_upper: Iterable[float] | None = None,
    threshold: float,
    max_slope: float,
) -> BurdenBounds:
    """Compute sharp bounds for hypoxia occupation and cumulative deficit.

    Pass either exact ``values`` or both ``observation_lower`` and
    ``observation_upper``.  Time units determine the units of occupation and
    oxygen-deficit outputs; ``max_slope`` must use the matching time unit.
    """

    if not math.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if values is not None:
        if observation_lower is not None or observation_upper is not None:
            raise ValueError("pass exact values or interval bounds, not both")
        exact = list(values)
        lower_values = exact
        upper_values = exact
        mode = "exact"
    else:
        if observation_lower is None or observation_upper is None:
            raise ValueError("both observation interval bounds are required")
        lower_values = list(observation_lower)
        upper_values = list(observation_upper)
        mode = "interval"

    normalized_times, _, _ = _prepare_inputs(
        times,
        lower_values,
        upper_values,
        max_slope,
    )
    lower_envelope, upper_envelope = construct_envelopes(
        normalized_times,
        lower_values,
        upper_values,
        max_slope,
    )
    horizon = float(normalized_times[-1])
    occupation_lower = integrate_occupation(upper_envelope, threshold)
    occupation_upper = integrate_occupation(lower_envelope, threshold)
    deficit_lower = integrate_oxygen_deficit(upper_envelope, threshold)
    deficit_upper = integrate_oxygen_deficit(lower_envelope, threshold)
    return BurdenBounds(
        horizon=horizon,
        threshold=float(threshold),
        max_slope=float(max_slope),
        occupation_lower=occupation_lower,
        occupation_upper=occupation_upper,
        occupation_fraction_lower=occupation_lower / horizon,
        occupation_fraction_upper=occupation_upper / horizon,
        oxygen_deficit_lower=deficit_lower,
        oxygen_deficit_upper=deficit_upper,
        event_count_status="not_bounded_in_general_by_lipschitz",
        observation_mode=mode,
    )


def _prototype_scenarios() -> list[dict[str, object]]:
    scenarios = [
        {
            "name": "above_threshold_short_gap",
            "times": [0.0, 1.0],
            "values": [6.0, 6.0],
            "threshold": 5.0,
            "max_slope": 1.0,
        },
        {
            "name": "hidden_hypoxia_possible",
            "times": [0.0, 2.0],
            "values": [6.0, 6.0],
            "threshold": 5.0,
            "max_slope": 2.0,
        },
        {
            "name": "certain_hypoxia",
            "times": [0.0, 1.0],
            "values": [4.0, 4.0],
            "threshold": 5.0,
            "max_slope": 1.0,
        },
        {
            "name": "measurement_interval",
            "times": [0.0, 2.0],
            "observation_lower": [5.5, 5.5],
            "observation_upper": [6.5, 6.5],
            "threshold": 5.0,
            "max_slope": 1.0,
        },
    ]
    results: list[dict[str, object]] = []
    for scenario in scenarios:
        arguments = dict(scenario)
        name = str(arguments.pop("name"))
        results.append({"name": name, **asdict(identify_hypoxia_burden(**arguments))})
    return results


def write_prototype_outputs(output_directory: Path = DEFAULT_OUT) -> Path:
    output_directory.mkdir(parents=True, exist_ok=True)
    results = _prototype_scenarios()
    json_path = output_directory / "partial_identification_scenarios.json"
    json_path.write_text(
        json.dumps(
            {
                "status": "theory_prototype_only",
                "guarantee_scope": "globally_L_lipschitz_feasible_paths",
                "event_count": "not_finitely_bounded_without_extra_assumptions",
                "scenarios": results,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    report_lines = [
        "# Partial-identification prototype",
        "",
        "Status: `theory_prototype_only`.",
        "",
        "The outputs are exact for the stated global Lipschitz feasible-path class. "
        "They are not empirical estimates and do not validate a value of L.",
        "",
        "| Scenario | Occupation interval | Oxygen-deficit interval |",
        "|---|---:|---:|",
    ]
    for row in results:
        report_lines.append(
            f"| {row['name']} | "
            f"[{row['occupation_lower']:.6g}, {row['occupation_upper']:.6g}] | "
            f"[{row['oxygen_deficit_lower']:.6g}, {row['oxygen_deficit_upper']:.6g}] |"
        )
    report_lines.extend(
        [
            "",
            "Event-count upper bounds remain unbounded under a Lipschitz condition alone.",
        ]
    )
    (output_directory / "PARTIAL_IDENTIFICATION_PROTOTYPE_REPORT.md").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )
    return json_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    output = write_prototype_outputs(args.output_dir)
    print(output)


if __name__ == "__main__":
    main()
