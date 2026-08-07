# -*- coding: utf-8 -*-
"""Recompute the frozen 24-hour, H=5 slice with ``DO >= 0``."""

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
    identify_hypoxia_burden_with_replacements,
    minimum_replacements_for_lipschitz,
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
    / "current_maracoos_24h_nonnegative_trilemma"
)

EXPECTED_VALIDATION_IDS = (
    "mddnr_Bishopville_Prong",
    "mddnr_Camp_Tockwogh",
    "mddnr_Greys_Creek",
    "mddnr_Harris_Creek_Upstream",
)
WINDOW_POINTS = 24 * 4 + 1
MAX_GAP_HOURS = 18.0 / 60.0
THRESHOLD_MGL = 5.0
STATE_LOWER_BOUND_MGL = 0.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _frozen_items_and_slopes() -> tuple[
    list[dict[str, object]],
    list[tuple[str, float]],
    dict[str, object],
]:
    metadata = json.loads(V3_METADATA.read_text(encoding="utf-8"))
    validation_ids = tuple(metadata["validation_dataset_ids"])
    if validation_ids != EXPECTED_VALIDATION_IDS:
        raise ValueError("validation split differs from frozen v3 metadata")
    if _sha256(MANIFEST) != metadata["source_manifest_sha256"]:
        raise ValueError("manifest hash differs from frozen v3 metadata")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    selected = sorted(
        manifest["selected"],
        key=lambda item: str(item["dataset_id"]).casefold(),
    )
    selected_by_id = {str(item["dataset_id"]): item for item in selected}
    items = [selected_by_id[dataset_id] for dataset_id in validation_ids]

    inherited = metadata["slope_parameters_mgL_per_hour"]
    labelled_slopes = [
        ("grid_4", 4.0),
        ("grid_8", 8.0),
        ("calibration_q999", float(inherited["calibration_q999"])),
        ("grid_16", 16.0),
        ("calibration_max", float(inherited["calibration_max"])),
        ("grid_32", 32.0),
        ("grid_64", 64.0),
    ]
    labelled_slopes.sort(key=lambda item: item[1])
    return items, labelled_slopes, metadata


def _first_complete_day(frame: pd.DataFrame) -> pd.DataFrame:
    ordered = frame.sort_values("time", kind="stable").reset_index(drop=True)
    for start in range(0, len(ordered) - WINDOW_POINTS + 1):
        candidate = ordered.iloc[start : start + WINDOW_POINTS].copy()
        elapsed = (
            candidate["time"] - candidate["time"].iloc[0]
        ).dt.total_seconds().to_numpy(dtype=float) / 3600.0
        gaps = np.diff(elapsed)
        if (
            23.5 <= elapsed[-1] <= 24.5
            and np.all(gaps > 0.0)
            and np.all(gaps <= MAX_GAP_HOURS)
        ):
            candidate["elapsed_hours"] = elapsed
            return candidate
    raise ValueError("no complete 24-hour block")


def run() -> tuple[pd.DataFrame, dict[str, object]]:
    _pilot.DATA = DATA
    _pilot.MANIFEST = MANIFEST
    validation_items, labelled_slopes, source_metadata = (
        _frozen_items_and_slopes()
    )
    rows: list[dict[str, object]] = []
    started = time.perf_counter()
    for item in validation_items:
        frame, _audit = _pilot.load_station(item)
        block = _first_complete_day(frame)
        dataset_id = str(item["dataset_id"])
        times = block["elapsed_hours"].to_numpy(dtype=float)
        values = block["do_mgL"].to_numpy(dtype=float)
        previous_k: int | None = None
        for slope_kind, max_slope in labelled_slopes:
            case_started = time.perf_counter()
            feasibility = minimum_replacements_for_lipschitz(
                times,
                values,
                max_slope,
            )
            floor_feasibility = minimum_replacements_for_lipschitz(
                times,
                values,
                max_slope,
                state_lower_bound=STATE_LOWER_BOUND_MGL,
            )
            if (
                floor_feasibility.minimum_replacements
                != feasibility.minimum_replacements
            ):
                raise AssertionError(
                    "nonnegative reports must have the same floor-aware k_min"
                )
            if (
                previous_k is not None
                and feasibility.minimum_replacements > previous_k
            ):
                raise AssertionError("k_min(L) must be non-increasing")
            previous_k = feasibility.minimum_replacements
            budget = feasibility.minimum_replacements

            unconstrained = identify_hypoxia_burden_with_replacements(
                times,
                values,
                threshold=THRESHOLD_MGL,
                max_slope=max_slope,
                max_replacements=budget,
            )
            constrained = identify_hypoxia_burden_with_replacements(
                times,
                values,
                threshold=THRESHOLD_MGL,
                max_slope=max_slope,
                max_replacements=budget,
                state_lower_bound=STATE_LOWER_BOUND_MGL,
            )
            unchanged_pairs = (
                (
                    unconstrained.occupation_lower,
                    constrained.occupation_lower,
                ),
                (
                    unconstrained.occupation_upper,
                    constrained.occupation_upper,
                ),
                (
                    unconstrained.oxygen_deficit_lower,
                    constrained.oxygen_deficit_lower,
                ),
            )
            if any(
                not math.isclose(
                    left, right, rel_tol=1e-11, abs_tol=1e-11
                )
                for left, right in unchanged_pairs
            ):
                raise AssertionError(
                    "positive-threshold floor changed an invariant endpoint"
                )
            horizon = constrained.horizon
            unconstrained_width = (
                unconstrained.oxygen_deficit_upper
                - unconstrained.oxygen_deficit_lower
            ) / (THRESHOLD_MGL * horizon)
            constrained_width = (
                constrained.oxygen_deficit_upper
                - constrained.oxygen_deficit_lower
            ) / (THRESHOLD_MGL * horizon)
            if constrained_width > 1.0 + 1e-12:
                raise AssertionError("nonnegative deficit width exceeds H*T")

            rows.append(
                {
                    "dataset_id": dataset_id,
                    "block_start": str(block["time"].iloc[0]),
                    "block_end": str(block["time"].iloc[-1]),
                    "n_observations": len(block),
                    "threshold_mgL": THRESHOLD_MGL,
                    "state_lower_bound_mgL": STATE_LOWER_BOUND_MGL,
                    "L_kind": slope_kind,
                    "L_mgL_per_hour": max_slope,
                    "k_min": budget,
                    "k_min_fraction": budget / len(block),
                    "occupation_lower_hours": (
                        constrained.occupation_lower
                    ),
                    "occupation_upper_hours": (
                        constrained.occupation_upper
                    ),
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
                    "feasibility_deleted_indices": json.dumps(
                        feasibility.deleted_indices
                    ),
                    "deficit_upper_deleted_indices_nonnegative": json.dumps(
                        constrained.oxygen_deficit_upper_witness.deleted_indices
                    ),
                    "runtime_seconds": time.perf_counter() - case_started,
                }
            )

    detail = pd.DataFrame(rows)
    metadata = {
        "status": "complete_nonconfirmatory_mechanism_demonstration",
        "design": (
            "H=5 slice; earliest_complete_24h_block_per_validation_station"
        ),
        "threshold_mgL": THRESHOLD_MGL,
        "state_lower_bound_mgL": STATE_LOWER_BOUND_MGL,
        "validation_ids": list(EXPECTED_VALIDATION_IDS),
        "L_grid_mgL_per_hour": [
            {"label": label, "value": value}
            for label, value in labelled_slopes
        ],
        "source_manifest_sha256": source_metadata[
            "source_manifest_sha256"
        ],
        "source_v3_metadata_sha256": _sha256(V3_METADATA),
        "validation_data_used_to_select_or_repair_L": False,
        "interpretation": (
            "For each L, k is the minimum discrete feasibility budget. "
            "Deleted reports are witnesses, not verified errors; empirical "
            "slopes do not validate a continuous-path physical L."
        ),
        "runtime_source_sha256": _sha256(
            HERE / "robust_exact_current.py"
        ),
        "analysis_source_sha256": _sha256(Path(__file__)),
        "total_runtime_seconds": time.perf_counter() - started,
        "case_count": len(detail),
        "maximum_nonnegative_deficit_width_fraction": float(
            detail[
                "deficit_width_fraction_nonnegative_H_times_T"
            ].max()
        ),
    }
    return detail, metadata


def write_outputs() -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    detail, metadata = run()
    detail.to_csv(OUTPUT / "trilemma_nonnegative_detail.csv", index=False)
    (OUTPUT / "trilemma_nonnegative_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# MARACOOS/MDDNR 24-hour physical-state L-k-width surface: H=5 slice",
        "",
        "Status: `nonconfirmatory mechanism demonstration`.",
        "",
        "This output recomputes only the frozen `H=5 mg/L` slice. The station, "
        "day, L grid and k-min rule are unchanged, and the new estimand "
        "imposes `x(t) >= 0 mg/L`. The inherited `H=2 mg/L` slice is not "
        "recomputed here.",
        "",
        "| Station | L source | L | k/n | Occupation width/T | "
        "Unconstrained deficit width/(HT) | Nonnegative deficit width/(HT) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in detail.itertuples(index=False):
        lines.append(
            f"| {row.dataset_id} | {row.L_kind} | "
            f"{row.L_mgL_per_hour:.4g} | {row.k_min_fraction:.4f} | "
            f"{row.occupation_width_fraction:.3f} | "
            f"{row.deficit_width_fraction_unconstrained_H_times_T:.3f} | "
            f"{row.deficit_width_fraction_nonnegative_H_times_T:.3f} |"
        )
    lines.extend(
        [
            "",
            "This surface is conditional sensitivity evidence, not a claim "
            "that the selected L or deleted reports are physically true.",
            "",
            f"Runtime: {metadata['total_runtime_seconds']:.2f} s.",
        ]
    )
    report = OUTPUT / "MARACOOS_24H_NONNEGATIVE_TRILEMMA_REPORT.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    print(write_outputs())
