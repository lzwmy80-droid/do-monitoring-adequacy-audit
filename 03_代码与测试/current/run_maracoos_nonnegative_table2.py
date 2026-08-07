# -*- coding: utf-8 -*-
"""Recompute the frozen 7-day Table 2 slice with the physical DO floor.

This script preserves the frozen station split, collision-aware loader,
calibration-only q99.9 slope, window rule, cadences, and all index phases.  It
changes only the latent state space by comparing the unconstrained model with
the physically constrained model ``x(t) >= 0 mg/L``.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from robust_exact_current import (
    InsufficientReplacementBudget,
    identify_hypoxia_burden_with_replacements,
)


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SNAPSHOT = HERE.parent / "原始快照"
if str(SNAPSHOT) not in sys.path:
    sys.path.insert(0, str(SNAPSHOT))

import future_partial_identification_maracoos_pilot as _pilot  # noqa: E402


DATA = PROJECT / "04_真实数据" / "MARACOOS_MDDNR_discovery_v3"
MANIFEST = DATA / "download_manifest.json"
V3_METADATA = (
    PROJECT
    / "05_既有结果"
    / "partial_identification_maracoos_pilot_v3_collision_identity"
    / "pilot_metadata.json"
)
OUTPUT = (
    PROJECT
    / "05_既有结果"
    / "current_maracoos_7day_nonnegative_table2"
)

EXPECTED_CALIBRATION_IDS = (
    "mddnr_Arundel_on_the_Bay",
    "mddnr_Budds_Landing",
    "mddnr_Dares_Beach",
    "mddnr_Harris_Creek_Downstream",
)
EXPECTED_VALIDATION_IDS = (
    "mddnr_Bishopville_Prong",
    "mddnr_Camp_Tockwogh",
    "mddnr_Greys_Creek",
    "mddnr_Harris_Creek_Upstream",
)

THRESHOLD_MGL = 5.0
STATE_LOWER_BOUND_MGL = 0.0
CADENCE_HOURS = (0.5, 1.0, 2.0, 4.0)
BASE_CADENCE_HOURS = 0.25
WINDOW_POINTS = 7 * 24 * 4 + 1
MIN_WINDOW_HOURS = 167.5
MAX_WINDOW_HOURS = 168.5
MAX_GAP_HOURS = 18.0 / 60.0
MAX_BLOCKS_PER_STATION = 4


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_frozen_inputs() -> tuple[list[dict[str, object]], float]:
    metadata = json.loads(V3_METADATA.read_text(encoding="utf-8"))
    calibration_ids = tuple(metadata["calibration_dataset_ids"])
    validation_ids = tuple(metadata["validation_dataset_ids"])
    if calibration_ids != EXPECTED_CALIBRATION_IDS:
        raise ValueError("calibration split differs from the frozen v3 split")
    if validation_ids != EXPECTED_VALIDATION_IDS:
        raise ValueError("validation split differs from the frozen v3 split")

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if _sha256(MANIFEST) != metadata["source_manifest_sha256"]:
        raise ValueError("manifest hash differs from frozen v3 metadata")
    selected = sorted(
        manifest["selected"],
        key=lambda item: str(item["dataset_id"]).casefold(),
    )
    selected_by_id = {str(item["dataset_id"]): item for item in selected}
    items = [selected_by_id[dataset_id] for dataset_id in validation_ids]
    slope = float(
        metadata["slope_parameters_mgL_per_hour"]["calibration_q999"]
    )
    if not math.isfinite(slope) or slope < 0.0:
        raise ValueError("frozen calibration q99.9 slope is invalid")
    return items, slope


def _find_blocks(frame: pd.DataFrame) -> list[pd.DataFrame]:
    ordered = frame.sort_values("time", kind="stable").reset_index(drop=True)
    blocks: list[pd.DataFrame] = []
    cursor = 0
    while (
        cursor + WINDOW_POINTS <= len(ordered)
        and len(blocks) < MAX_BLOCKS_PER_STATION
    ):
        candidate = ordered.iloc[cursor : cursor + WINDOW_POINTS].copy()
        elapsed = (
            candidate["time"] - candidate["time"].iloc[0]
        ).dt.total_seconds().to_numpy(dtype=float) / 3600.0
        adjacent = np.diff(elapsed)
        valid = (
            len(candidate) == WINDOW_POINTS
            and MIN_WINDOW_HOURS <= elapsed[-1] <= MAX_WINDOW_HOURS
            and np.all(adjacent > 0.0)
            and np.all(adjacent <= MAX_GAP_HOURS)
        )
        if valid:
            candidate["elapsed_hours"] = elapsed
            blocks.append(candidate)
            cursor += WINDOW_POINTS
        else:
            cursor += 1
    return blocks


def _reference_metrics(
    times: np.ndarray,
    values: np.ndarray,
    threshold: float,
) -> tuple[float, float, float]:
    occupation_terms: list[float] = []
    deficit_terms: list[float] = []
    slopes: list[float] = []
    for left_time, right_time, left_value, right_value in zip(
        times[:-1],
        times[1:],
        values[:-1],
        values[1:],
    ):
        duration = float(right_time - left_time)
        slopes.append(abs(float((right_value - left_value) / duration)))
        left_below = bool(left_value < threshold)
        right_below = bool(right_value < threshold)
        if left_below and right_below:
            occupation_terms.append(duration)
            deficit_terms.append(
                0.5
                * (
                    (threshold - float(left_value))
                    + (threshold - float(right_value))
                )
                * duration
            )
        elif not left_below and not right_below:
            occupation_terms.append(0.0)
            deficit_terms.append(0.0)
        else:
            crossing_fraction = float(
                (threshold - left_value) / (right_value - left_value)
            )
            crossing_fraction = min(1.0, max(0.0, crossing_fraction))
            if left_below:
                active_fraction = crossing_fraction
                active_height = threshold - float(left_value)
            else:
                active_fraction = 1.0 - crossing_fraction
                active_height = threshold - float(right_value)
            occupation_terms.append(duration * active_fraction)
            deficit_terms.append(
                0.5 * active_height * duration * active_fraction
            )
    return (
        float(math.fsum(occupation_terms)),
        float(math.fsum(deficit_terms)),
        max(slopes),
    )


def run() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    # The frozen loader resolves this module global at call time.
    _pilot.DATA = DATA
    _pilot.MANIFEST = MANIFEST
    validation_items, max_slope = _load_frozen_inputs()

    blocks: list[dict[str, object]] = []
    for item in validation_items:
        dataset_id = str(item["dataset_id"])
        frame, _audit = _pilot.load_station(item)
        for block_index, block in enumerate(_find_blocks(frame)):
            blocks.append(
                {
                    "dataset_id": dataset_id,
                    "block_index": block_index,
                    "block_start": str(block["time"].iloc[0]),
                    "block_end": str(block["time"].iloc[-1]),
                    "times": block["elapsed_hours"].to_numpy(dtype=float),
                    "values": block["do_mgL"].to_numpy(dtype=float),
                }
            )

    rows: list[dict[str, object]] = []
    started = time.perf_counter()
    for block in blocks:
        full_times = np.asarray(block["times"], dtype=float)
        full_values = np.asarray(block["values"], dtype=float)
        reference_block = _reference_metrics(
            full_times, full_values, THRESHOLD_MGL
        )
        block_reference_max_slope = reference_block[2]
        for cadence in CADENCE_HOURS:
            step = int(round(cadence / BASE_CADENCE_HOURS))
            for phase in range(step):
                sparse_indices = np.arange(
                    phase, full_times.size, step, dtype=int
                )
                if sparse_indices.size < 2:
                    continue
                first = int(sparse_indices[0])
                last = int(sparse_indices[-1])
                reference_times = full_times[first : last + 1]
                reference_values = full_values[first : last + 1]
                reference_occupation, reference_deficit, _phase_slope = (
                    _reference_metrics(
                        reference_times,
                        reference_values,
                        THRESHOLD_MGL,
                    )
                )
                sparse_times = full_times[sparse_indices]
                sparse_values = full_values[sparse_indices]
                horizon = float(reference_times[-1] - reference_times[0])
                base = {
                    "dataset_id": block["dataset_id"],
                    "block_index": block["block_index"],
                    "block_start": block["block_start"],
                    "block_end": block["block_end"],
                    "threshold_mgL": THRESHOLD_MGL,
                    "cadence_hours": cadence,
                    "phase_index": phase,
                    "max_slope_mgL_per_hour": max_slope,
                    "reference_L_compatible": (
                        block_reference_max_slope <= max_slope + 1e-12
                    ),
                    "horizon_hours": horizon,
                    "n_sparse_reports": int(sparse_indices.size),
                    "reference_occupation_hours": reference_occupation,
                    "reference_deficit_mgL_hours": reference_deficit,
                }
                case_started = time.perf_counter()
                try:
                    unconstrained = (
                        identify_hypoxia_burden_with_replacements(
                            sparse_times,
                            sparse_values,
                            threshold=THRESHOLD_MGL,
                            max_slope=max_slope,
                            max_replacements=0,
                        )
                    )
                    constrained = (
                        identify_hypoxia_burden_with_replacements(
                            sparse_times,
                            sparse_values,
                            threshold=THRESHOLD_MGL,
                            max_slope=max_slope,
                            max_replacements=0,
                            state_lower_bound=STATE_LOWER_BOUND_MGL,
                        )
                    )
                except InsufficientReplacementBudget as error:
                    rows.append(
                        {
                            **base,
                            "status": "infeasible",
                            "error": str(error),
                            "runtime_seconds": time.perf_counter()
                            - case_started,
                        }
                    )
                    continue

                unconstrained_width = (
                    unconstrained.oxygen_deficit_upper
                    - unconstrained.oxygen_deficit_lower
                ) / (THRESHOLD_MGL * horizon)
                constrained_width = (
                    constrained.oxygen_deficit_upper
                    - constrained.oxygen_deficit_lower
                ) / (THRESHOLD_MGL * horizon)
                rows.append(
                    {
                        **base,
                        "status": "feasible",
                        "error": "",
                        "occupation_lower": constrained.occupation_lower,
                        "occupation_upper": constrained.occupation_upper,
                        "occupation_width_fraction": (
                            constrained.occupation_upper
                            - constrained.occupation_lower
                        )
                        / horizon,
                        "deficit_lower_mgL_hours": (
                            constrained.oxygen_deficit_lower
                        ),
                        "deficit_upper_unconstrained_mgL_hours": (
                            unconstrained.oxygen_deficit_upper
                        ),
                        "deficit_upper_nonnegative_mgL_hours": (
                            constrained.oxygen_deficit_upper
                        ),
                        "deficit_width_fraction_unconstrained_H_times_T": (
                            unconstrained_width
                        ),
                        "deficit_width_fraction_nonnegative_H_times_T": (
                            constrained_width
                        ),
                        "deficit_width_fraction_reduction": (
                            unconstrained_width - constrained_width
                        ),
                        "occupation_reference_covered": (
                            constrained.occupation_lower - 1e-8
                            <= reference_occupation
                            <= constrained.occupation_upper + 1e-8
                        ),
                        "deficit_reference_covered": (
                            constrained.oxygen_deficit_lower - 1e-8
                            <= reference_deficit
                            <= constrained.oxygen_deficit_upper + 1e-8
                        ),
                        "runtime_seconds": time.perf_counter() - case_started,
                    }
                )

    detail = pd.DataFrame(rows)
    compatible = detail.loc[
        detail["status"].eq("feasible")
        & detail["reference_L_compatible"].eq(True)
    ]
    block_diagnostics = detail[
        [
            "dataset_id",
            "block_index",
            "block_start",
            "block_end",
            "reference_L_compatible",
        ]
    ].drop_duplicates()
    compatible_blocks = block_diagnostics.loc[
        block_diagnostics["reference_L_compatible"].eq(True)
    ]
    compatible_station_ids = sorted(
        compatible_blocks["dataset_id"].unique().tolist()
    )
    summary = (
        compatible.groupby("cadence_hours", dropna=False)
        .agg(
            compatible_feasible_cases=("status", "size"),
            median_occupation_width_fraction=(
                "occupation_width_fraction",
                "median",
            ),
            median_deficit_width_fraction_unconstrained=(
                "deficit_width_fraction_unconstrained_H_times_T",
                "median",
            ),
            median_deficit_width_fraction_nonnegative=(
                "deficit_width_fraction_nonnegative_H_times_T",
                "median",
            ),
            median_deficit_width_fraction_reduction=(
                "deficit_width_fraction_reduction",
                "median",
            ),
            maximum_deficit_width_fraction_nonnegative=(
                "deficit_width_fraction_nonnegative_H_times_T",
                "max",
            ),
            containment_failures=(
                "deficit_reference_covered",
                lambda values: int(values.eq(False).sum()),
            ),
            median_runtime_seconds=("runtime_seconds", "median"),
        )
        .reset_index()
    )
    metadata = {
        "status": "complete",
        "estimand": "exact_inlier_endpoints_with_known_DO_state_floor",
        "state_lower_bound_mgL": STATE_LOWER_BOUND_MGL,
        "threshold_mgL": THRESHOLD_MGL,
        "max_replacements": 0,
        "calibration_q999_max_slope_mgL_per_hour": max_slope,
        "validation_data_used_to_select_or_repair_L": False,
        "frozen_validation_dataset_ids": list(EXPECTED_VALIDATION_IDS),
        "window_rule": {
            "points": WINDOW_POINTS,
            "span_hours_inclusive": [
                MIN_WINDOW_HOURS,
                MAX_WINDOW_HOURS,
            ],
            "maximum_adjacent_gap_minutes": MAX_GAP_HOURS * 60.0,
            "maximum_blocks_per_station": MAX_BLOCKS_PER_STATION,
            "selection": "earliest qualifying row-disjoint windows",
        },
        "cadence_hours_all_phases": list(CADENCE_HOURS),
        "validation_block_count": len(blocks),
        "compatible_validation_block_count": len(compatible_blocks),
        "validation_station_count": len(EXPECTED_VALIDATION_IDS),
        "compatible_validation_station_count": len(
            compatible_station_ids
        ),
        "compatible_validation_station_ids": compatible_station_ids,
        "source_manifest": str(MANIFEST),
        "source_manifest_sha256": _sha256(MANIFEST),
        "source_v3_metadata": str(V3_METADATA),
        "source_v3_metadata_sha256": _sha256(V3_METADATA),
        "runtime_source": str(HERE / "robust_exact_current.py"),
        "runtime_source_sha256": _sha256(HERE / "robust_exact_current.py"),
        "analysis_source_sha256": _sha256(Path(__file__)),
        "total_runtime_seconds": time.perf_counter() - started,
        "case_count": len(detail),
        "feasible_case_count": int(detail["status"].eq("feasible").sum()),
        "compatible_feasible_case_count": len(compatible),
        "compatible_containment_failures": int(
            compatible["occupation_reference_covered"].eq(False).sum()
            + compatible["deficit_reference_covered"].eq(False).sum()
        ),
    }
    return detail, summary, metadata


def write_outputs() -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    detail, summary, metadata = run()
    detail.to_csv(OUTPUT / "table2_nonnegative_detail.csv", index=False)
    summary.to_csv(OUTPUT / "table2_nonnegative_summary.csv", index=False)
    (OUTPUT / "table2_nonnegative_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# MARACOOS/MDDNR 7-day physical-state reanalysis",
        "",
        "This is the frozen Table 2 slice recomputed with "
        "`x(t) >= 0 mg/L`. No validation observation was used to select or "
        "repair `L`.",
        "",
        f"Compatibility was met by "
        f"{metadata['compatible_validation_block_count']}/"
        f"{metadata['validation_block_count']} frozen blocks from "
        f"{metadata['compatible_validation_station_count']}/"
        f"{metadata['validation_station_count']} validation stations "
        f"({', '.join(metadata['compatible_validation_station_ids'])}), "
        f"yielding {metadata['compatible_feasible_case_count']}/"
        f"{metadata['case_count']} pooled station-window-cadence-phase cases. "
        "These pooled cases overlap across phases and are not independent "
        "windows. Width summaries below are conditional on those compatible "
        "feasible cases.",
        "",
        "Containment is a deterministic check against the retained 15-minute "
        "operational reference, not statistical coverage. All compatible "
        "occupation and deficit containment checks passed.",
        "",
        "| Cadence h | Compatible feasible cases | Occupation width/T | "
        "Unconstrained deficit width/(HT) | Nonnegative deficit width/(HT) | "
        "Median reduction | Maximum nonnegative width | "
        "Deficit containment failures |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.cadence_hours:g} | "
            f"{int(row.compatible_feasible_cases)} | "
            f"{row.median_occupation_width_fraction:.3f} | "
            f"{row.median_deficit_width_fraction_unconstrained:.3f} | "
            f"{row.median_deficit_width_fraction_nonnegative:.3f} | "
            f"{row.median_deficit_width_fraction_reduction:.3f} | "
            f"{row.maximum_deficit_width_fraction_nonnegative:.3f} | "
            f"{int(row.containment_failures)} |"
        )
    lines.extend(
        [
            "",
            f"Runtime: {metadata['total_runtime_seconds']:.2f} s for "
            f"{metadata['case_count']} pooled station-window-cadence-phase cases.",
        ]
    )
    report = OUTPUT / "MARACOOS_7DAY_NONNEGATIVE_TABLE2_REPORT.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    print(write_outputs())
