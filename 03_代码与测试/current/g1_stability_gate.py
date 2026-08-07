# -*- coding: utf-8 -*-
"""Falsification gate for the DIR-01 robust deficit-stability theorem.

The imported dynamic program is a frozen source snapshot.  This file adds no
new estimator: it tests a deterministic upper bound implied by the overlap of
two admissible retained-index sets.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import sys
from pathlib import Path
from typing import Iterable, Sequence


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SNAPSHOT = HERE.parent / "原始快照"
if str(SNAPSHOT) not in sys.path:
    sys.path.insert(0, str(SNAPSHOT))

from robust_exact_current import (  # noqa: E402
    identify_hypoxia_burden_with_replacements,
)


def _validated_times(times: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in times)
    if len(values) < 2 or any(not math.isfinite(value) for value in values):
        raise ValueError("times must contain at least two finite values")
    if any(right <= left for left, right in zip(values[:-1], values[1:])):
        raise ValueError("times must be strictly increasing")
    return values


def deficit_diameter_upper_bound(
    times: Iterable[float],
    *,
    max_slope: float,
    max_replacements: int,
    half_error: float = 0.0,
) -> dict[str, float | int]:
    """Return the G1 bound for symmetric deterministic ``+/- half_error`` limits.

    The bound applies only when ``n > 2k``.  This guarantees that any two
    admissible deletion sets leave at least one common retained observation.
    """

    t = _validated_times(times)
    if not math.isfinite(max_slope) or max_slope < 0:
        raise ValueError("max_slope must be finite and non-negative")
    if isinstance(max_replacements, bool) or not isinstance(max_replacements, int):
        raise TypeError("max_replacements must be an integer")
    if max_replacements < 0:
        raise ValueError("max_replacements must be non-negative")
    if not math.isfinite(half_error) or half_error < 0:
        raise ValueError("half_error must be finite and non-negative")
    if len(t) <= 2 * max_replacements:
        raise ValueError("the finite bound requires n > 2k")

    gaps = tuple(right - left for left, right in zip(t[:-1], t[1:]))
    horizon = t[-1] - t[0]
    mesh = max(gaps)
    measurement_term = 2.0 * half_error * horizon
    adaptive_mesh_term = 0.5 * sum(gap * gap for gap in gaps)
    contamination_term = (
        4 * max_replacements**2 - max_replacements
    ) * mesh**2
    adaptive_bound = measurement_term + max_slope * (
        adaptive_mesh_term + contamination_term
    )
    coarse_bound = measurement_term + max_slope * (
        0.5 * horizon * mesh
        + (4.0 * max_replacements**2 - max_replacements) * mesh**2
    )
    return {
        "n": len(t),
        "k": max_replacements,
        "horizon": horizon,
        "mesh": mesh,
        "measurement_term": measurement_term,
        "adaptive_mesh_term_inside_L": adaptive_mesh_term,
        "contamination_term_inside_L": contamination_term,
        "adaptive_bound": adaptive_bound,
        "coarse_bound": coarse_bound,
    }


def regular_grid_deficit_diameter_upper_bound(
    times: Iterable[float],
    *,
    max_slope: float,
    max_replacements: int,
    half_error: float = 0.0,
) -> float:
    """Return the sharper regular-grid bound.

    With ``n > 2k``, deleting the first ``2k`` common anchors is the worst
    missing-anchor pattern (ties occur at ``k=1``).  Relative to the complete
    grid its exact geometry increment is ``(4k^2-k)h^2``.
    """

    t = _validated_times(times)
    base = deficit_diameter_upper_bound(
        t,
        max_slope=max_slope,
        max_replacements=max_replacements,
        half_error=half_error,
    )
    gaps = tuple(right - left for left, right in zip(t[:-1], t[1:]))
    mesh = gaps[0]
    tolerance = 1e-11 * max(1.0, abs(mesh))
    if any(abs(gap - mesh) > tolerance for gap in gaps[1:]):
        raise ValueError("times must form a regular grid")
    horizon = t[-1] - t[0]
    k = max_replacements
    return float(base["measurement_term"]) + max_slope * (
        0.5 * horizon * mesh + (4 * k * k - k) * mesh * mesh
    )


def twice_distance_integral(times: Sequence[float], common_indices: Sequence[int]) -> float:
    """Compute ``2 * integral distance(t, common anchors) dt`` exactly."""

    t = _validated_times(times)
    common = tuple(sorted(set(int(index) for index in common_indices)))
    if not common or common[0] < 0 or common[-1] >= len(t):
        raise ValueError("common_indices must be a nonempty subset of grid indices")
    left = t[common[0]] - t[0]
    right = t[-1] - t[common[-1]]
    value = left * left + right * right
    value += 0.5 * sum(
        (t[b] - t[a]) ** 2 for a, b in zip(common[:-1], common[1:])
    )
    return value


def geometry_bound(times: Sequence[float], common_indices: Sequence[int]) -> float:
    """Bound the twice-distance integral using the number of missing anchors."""

    t = _validated_times(times)
    common = tuple(sorted(set(int(index) for index in common_indices)))
    if not common:
        raise ValueError("at least one common anchor is required")
    gaps = tuple(right - left for left, right in zip(t[:-1], t[1:]))
    missing = len(t) - len(common)
    if missing == 0:
        penalty = 0.0
    elif missing == 1:
        penalty = max(gaps) ** 2
    else:
        penalty = (missing**2 - 0.5 * missing) * max(gaps) ** 2
    return 0.5 * sum(gap * gap for gap in gaps) + penalty


def audit_geometry(number_of_grids: int, seed: int) -> dict[str, float | int]:
    """Enumerate all nonempty common-anchor sets on random small grids."""

    rng = random.Random(seed)
    checked = 0
    maximum_ratio = 0.0
    maximum_excess = -math.inf
    for _ in range(number_of_grids):
        n = rng.randint(2, 10)
        times = [0.0]
        for _index in range(1, n):
            times.append(times[-1] + rng.uniform(0.05, 1.5))
        for retained_count in range(1, n + 1):
            for common in itertools.combinations(range(n), retained_count):
                actual = twice_distance_integral(times, common)
                bound = geometry_bound(times, common)
                maximum_excess = max(maximum_excess, actual - bound)
                if bound > 0:
                    maximum_ratio = max(maximum_ratio, actual / bound)
                if actual > bound + 1e-10 * max(1.0, bound):
                    raise AssertionError(
                        f"geometry counterexample: actual={actual}, bound={bound}, "
                        f"times={times}, common={common}"
                    )
                checked += 1
    return {
        "random_grids": number_of_grids,
        "common_anchor_subsets_checked": checked,
        "maximum_actual_to_bound_ratio": maximum_ratio,
        "maximum_actual_minus_bound": maximum_excess,
    }


def _random_lipschitz_path(
    times: Sequence[float], max_slope: float, rng: random.Random
) -> list[float]:
    values = [rng.uniform(-3.0, 3.0)]
    for left, right in zip(times[:-1], times[1:]):
        values.append(
            values[-1] + rng.uniform(-max_slope, max_slope) * (right - left)
        )
    return values


def audit_exact_dp(number_of_cases: int, seed: int) -> dict[str, float | int]:
    """Search exact-observation DP outputs for a violation of the G1 bound."""

    rng = random.Random(seed)
    maximum_ratio = 0.0
    maximum_excess = -math.inf
    for case_index in range(number_of_cases):
        n = rng.randint(2, 12)
        k = rng.randint(0, (n - 1) // 2)
        times = [0.0]
        for _index in range(1, n):
            times.append(times[-1] + rng.uniform(0.05, 1.5))
        max_slope = rng.uniform(0.05, 5.0)
        reports = _random_lipschitz_path(times, max_slope, rng)
        contaminated_count = rng.randint(0, k)
        for index in rng.sample(range(n), contaminated_count):
            reports[index] += rng.choice((-1.0, 1.0)) * rng.uniform(0.0, 30.0)
        threshold = rng.uniform(min(reports) - 2.0, max(reports) + 2.0)

        result = identify_hypoxia_burden_with_replacements(
            times,
            reports,
            threshold=threshold,
            max_slope=max_slope,
            max_replacements=k,
        )
        width = result.oxygen_deficit_upper - result.oxygen_deficit_lower
        bound = float(
            deficit_diameter_upper_bound(
                times,
                max_slope=max_slope,
                max_replacements=k,
            )["adaptive_bound"]
        )
        maximum_excess = max(maximum_excess, width - bound)
        if bound > 0:
            maximum_ratio = max(maximum_ratio, width / bound)
        if width > bound + 2e-9 * max(1.0, bound):
            raise AssertionError(
                f"DP counterexample at case {case_index}: width={width}, "
                f"bound={bound}, n={n}, k={k}, L={max_slope}"
            )
    return {
        "exact_dp_cases_checked": number_of_cases,
        "maximum_width_to_bound_ratio": maximum_ratio,
        "maximum_width_minus_bound": maximum_excess,
    }


def audit_regular_exact_dp(number_of_cases: int, seed: int) -> dict[str, float | int]:
    """Search regular grids for a violation of the sharper coefficient."""

    rng = random.Random(seed)
    maximum_ratio = 0.0
    maximum_excess = -math.inf
    for case_index in range(number_of_cases):
        n = rng.randint(2, 14)
        k = rng.randint(0, (n - 1) // 2)
        mesh = rng.uniform(0.05, 1.5)
        times = [index * mesh for index in range(n)]
        max_slope = rng.uniform(0.05, 5.0)
        reports = _random_lipschitz_path(times, max_slope, rng)
        for index in rng.sample(range(n), rng.randint(0, k)):
            reports[index] += rng.choice((-1.0, 1.0)) * rng.uniform(0.0, 30.0)
        threshold = rng.uniform(min(reports) - 2.0, max(reports) + 2.0)
        result = identify_hypoxia_burden_with_replacements(
            times,
            reports,
            threshold=threshold,
            max_slope=max_slope,
            max_replacements=k,
        )
        width = result.oxygen_deficit_upper - result.oxygen_deficit_lower
        bound = regular_grid_deficit_diameter_upper_bound(
            times, max_slope=max_slope, max_replacements=k
        )
        maximum_excess = max(maximum_excess, width - bound)
        if bound > 0:
            maximum_ratio = max(maximum_ratio, width / bound)
        if width > bound + 2e-9 * max(1.0, bound):
            raise AssertionError(
                f"regular-grid counterexample at case {case_index}: width={width}, "
                f"bound={bound}, n={n}, k={k}, L={max_slope}"
            )
    return {
        "regular_exact_dp_cases_checked": number_of_cases,
        "maximum_width_to_bound_ratio": maximum_ratio,
        "maximum_width_minus_bound": maximum_excess,
    }


def boundary_order_witness(max_k: int = 8) -> list[dict[str, float | int]]:
    """Construct witnesses showing that the contamination term needs order k^2 h^2."""

    rows: list[dict[str, float | int]] = []
    max_slope = 2.0
    mesh = 1.0
    for k in range(1, max_k + 1):
        anchor = 2.0 * k * mesh
        times = [index * mesh for index in range(2 * k + 1)]
        reports = []
        for index, time in enumerate(times):
            if index < k:
                reports.append(max_slope * (anchor - time))
            elif index < 2 * k:
                reports.append(-max_slope * (anchor - time))
            else:
                reports.append(0.0)
        threshold = max_slope * anchor + 1.0
        result = identify_hypoxia_burden_with_replacements(
            times,
            reports,
            threshold=threshold,
            max_slope=max_slope,
            max_replacements=k,
        )
        width = result.oxygen_deficit_upper - result.oxygen_deficit_lower
        order_term = 4.0 * max_slope * k**2 * mesh**2
        sharp_bound = regular_grid_deficit_diameter_upper_bound(
            times, max_slope=max_slope, max_replacements=k
        )
        if not math.isclose(width, order_term, rel_tol=1e-11, abs_tol=1e-11):
            raise AssertionError(
                f"boundary witness failed for k={k}: width={width}, expected={order_term}"
            )
        rows.append(
            {
                "k": k,
                "n": len(times),
                "robust_deficit_width": width,
                "four_L_k2_h2": order_term,
                "regular_sharp_bound": sharp_bound,
            }
        )
    return rows


def extended_regular_sharpness_witness(
    max_k: int = 5,
    extra_interval_counts: Sequence[int] = (0, 1, 3, 8),
) -> list[dict[str, float | int]]:
    """Attain the full regular-grid bound for several horizons ``N >= 2k``.

    The first ``2k`` reports split into two mutually incompatible Lipschitz
    ramps.  Each admissible explanation deletes one half of that split.  All
    later reports equal the shared anchor value, so the exact endpoint solver
    supplies opposite tents on every remaining interval.  With a threshold
    above both extremal paths, deficit differences equal signed-integral
    differences and attain the complete regular-grid diameter formula.
    """

    if isinstance(max_k, bool) or not isinstance(max_k, int) or max_k < 1:
        raise ValueError("max_k must be a positive integer")
    extras = tuple(int(value) for value in extra_interval_counts)
    if any(value < 0 for value in extras):
        raise ValueError("extra_interval_counts must be non-negative")

    rows: list[dict[str, float | int]] = []
    max_slope = 2.0
    mesh = 1.0
    for k in range(1, max_k + 1):
        anchor = 2.0 * k * mesh
        for extra_intervals in extras:
            interval_count = 2 * k + extra_intervals
            times = [index * mesh for index in range(interval_count + 1)]
            reports = []
            for index, time in enumerate(times):
                if index < k:
                    reports.append(max_slope * (anchor - time))
                elif index < 2 * k:
                    reports.append(-max_slope * (anchor - time))
                else:
                    reports.append(0.0)
            threshold = max_slope * anchor + 1.0
            result = identify_hypoxia_burden_with_replacements(
                times,
                reports,
                threshold=threshold,
                max_slope=max_slope,
                max_replacements=k,
            )
            width = result.oxygen_deficit_upper - result.oxygen_deficit_lower
            sharp_bound = regular_grid_deficit_diameter_upper_bound(
                times, max_slope=max_slope, max_replacements=k
            )
            expected = max_slope * mesh**2 * (
                interval_count / 2.0 + 4.0 * k**2 - k
            )
            if not math.isclose(width, expected, rel_tol=1e-11, abs_tol=1e-11):
                raise AssertionError(
                    "extended sharpness witness failed: "
                    f"k={k}, N={interval_count}, width={width}, expected={expected}"
                )
            if not math.isclose(width, sharp_bound, rel_tol=1e-11, abs_tol=1e-11):
                raise AssertionError(
                    "witness did not attain the regular-grid bound: "
                    f"k={k}, N={interval_count}, width={width}, bound={sharp_bound}"
                )
            rows.append(
                {
                    "k": k,
                    "point_count": len(times),
                    "interval_count_N": interval_count,
                    "extra_intervals_after_first_common_anchor": extra_intervals,
                    "robust_deficit_width": width,
                    "regular_sharp_bound": sharp_bound,
                    "width_minus_bound": width - sharp_bound,
                }
            )
    return rows


def no_common_anchor_counterexample(magnitude: float = 100.0) -> dict[str, float | int]:
    """Show that no grid-only finite bound exists when ``n <= 2k``."""

    result = identify_hypoxia_burden_with_replacements(
        [0.0, 1.0],
        [-magnitude, magnitude],
        threshold=0.0,
        max_slope=0.0,
        max_replacements=1,
    )
    width = result.oxygen_deficit_upper - result.oxygen_deficit_lower
    if not math.isclose(width, magnitude, rel_tol=0.0, abs_tol=1e-12):
        raise AssertionError(f"unexpected no-anchor counterexample width: {width}")
    return {
        "n": 2,
        "k": 1,
        "L": 0.0,
        "epsilon": 0.0,
        "report_magnitude": magnitude,
        "robust_deficit_width": width,
        "grid_only_candidate_rhs": 0.0,
    }


def nonconnected_range_counterexample(magnitude: float = 100.0) -> dict[str, object]:
    """Record an exact two-point functional image whose hull is an interval.

    For two conflicting reports, ``L=0`` and ``k=1``, every feasible path is
    either the constant ``-magnitude`` or ``+magnitude`` path.  At threshold
    zero the deficit image is therefore ``{0, magnitude}``, not its convex
    hull.  Endpoint algorithms remain sharp for extrema but do not establish
    connectedness of the functional image.
    """

    if not math.isfinite(magnitude) or magnitude <= 0:
        raise ValueError("magnitude must be finite and positive")
    return {
        "times": [0.0, 1.0],
        "reports": [-magnitude, magnitude],
        "L": 0.0,
        "epsilon": 0.0,
        "k": 1,
        "threshold": 0.0,
        "exact_deficit_image": [0.0, magnitude],
        "endpoint_hull": [0.0, magnitude],
        "midpoint_in_exact_image": False,
        "interpretation": "sharp endpoints do not imply a connected identified set",
    }


def run_gate(dp_cases: int, geometry_grids: int, seed: int) -> dict[str, object]:
    return {
        "gate": "DIR-01-G1-deficit-stability",
        "status": "pass_no_counterexample_in_finite_search",
        "scope": "exact-observation DP plus deterministic grid geometry; not a proof by testing",
        "required_condition": "n > 2k",
        "candidate_adaptive_bound": (
            "2*epsilon*T + L*(0.5*sum(Delta_i^2) + (4*k^2-k)*h^2)"
        ),
        "candidate_coarse_bound": (
            "2*epsilon*T + L*(T*h/2 + (4*k^2-k)*h^2)"
        ),
        "regular_grid_sharp_bound": (
            "2*epsilon*T + L*(T*h/2 + (4*k^2-k)*h^2)"
        ),
        "geometry_audit": audit_geometry(geometry_grids, seed),
        "exact_dp_audit": audit_exact_dp(dp_cases, seed + 1),
        "regular_exact_dp_audit": audit_regular_exact_dp(dp_cases, seed + 2),
        "boundary_order_witnesses": boundary_order_witness(),
        "extended_regular_sharpness_witnesses": (
            extended_regular_sharpness_witness()
        ),
        "no_common_anchor_counterexample": no_common_anchor_counterexample(),
        "nonconnected_range_counterexample": nonconnected_range_counterexample(),
        "interpretation": (
            "Finite search did not falsify either valid bound. The boundary "
            "construction attains the regular-grid formula, while n<=2k gives "
            "an explicit counterexample to any bound lacking a shared anchor."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dp-cases", type=int, default=10_000)
    parser.add_argument("--geometry-grids", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260723)
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT / "03_代码与测试" / "gate_outputs" / "g1_stability_gate.json",
    )
    args = parser.parse_args()
    result = run_gate(args.dp_cases, args.geometry_grids, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
