# -*- coding: utf-8 -*-
"""Run the frozen seven-day fixed-k=1 sensitivity analysis.

The analysis reuses exactly the 240 feasible, reference-L-compatible cases
from the frozen k=0 Table 2 analysis and changes only the replacement budget
to the smallest positive fixed value, k=1. It evaluates the physical
nonnegative DO model only.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

import run_maracoos_nonnegative_table2 as baseline
from robust_exact_current import (
    InsufficientReplacementBudget,
    identify_hypoxia_burden_with_replacements,
)


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
BASELINE_OUTPUT = (
    PROJECT / "05_既有结果" / "current_maracoos_7day_nonnegative_table2"
)
OUTPUT = (
    PROJECT / "05_既有结果" / "current_maracoos_7day_fixed_k1_sensitivity"
)
FREEZE_NOTE = (
    PROJECT
    / "09_论文写作输出"
    / "K1_SENSITIVITY_FREEZE_NOTE_20260727.md"
)
CORE_RUNTIME = HERE / "robust_exact_current.py"
BASELINE_DETAIL = BASELINE_OUTPUT / "table2_nonnegative_detail.csv"
MAX_REPLACEMENTS = 1
EXPECTED_CASES = 240
TOLERANCE = 1e-8


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _case_key(row: pd.Series) -> tuple[str, int, float, int]:
    return (
        str(row["dataset_id"]),
        int(row["block_index"]),
        float(row["cadence_hours"]),
        int(row["phase_index"]),
    )


def _deleted_count(result: object, name: str) -> int:
    witness = getattr(result, name)
    return len(witness.deleted_indices)


def _deleted_indices(result: object, name: str) -> str:
    witness = getattr(result, name)
    return json.dumps(list(witness.deleted_indices), separators=(",", ":"))


def run() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not FREEZE_NOTE.exists():
        raise FileNotFoundError("the pre-computation analysis freeze note is missing")

    baseline_detail = pd.read_csv(BASELINE_DETAIL)
    selected = baseline_detail.loc[
        baseline_detail["status"].eq("feasible")
        & baseline_detail["reference_L_compatible"].eq(True)
    ].copy()
    if len(selected) != EXPECTED_CASES:
        raise ValueError(
            f"expected {EXPECTED_CASES} frozen compatible cases, found {len(selected)}"
        )
    lookup = {_case_key(row): row for _, row in selected.iterrows()}
    if len(lookup) != EXPECTED_CASES:
        raise ValueError("frozen case keys are not unique")

    baseline._pilot.DATA = baseline.DATA
    baseline._pilot.MANIFEST = baseline.MANIFEST
    validation_items, max_slope = baseline._load_frozen_inputs()

    blocks: list[dict[str, object]] = []
    for item in validation_items:
        dataset_id = str(item["dataset_id"])
        frame, _audit = baseline._pilot.load_station(item)
        for block_index, block in enumerate(baseline._find_blocks(frame)):
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
        block_reference = baseline._reference_metrics(
            full_times, full_values, baseline.THRESHOLD_MGL
        )
        block_reference_max_slope = block_reference[2]
        for cadence in baseline.CADENCE_HOURS:
            step = int(round(cadence / baseline.BASE_CADENCE_HOURS))
            for phase in range(step):
                key = (
                    str(block["dataset_id"]),
                    int(block["block_index"]),
                    float(cadence),
                    int(phase),
                )
                if key not in lookup:
                    continue
                k0 = lookup[key]
                sparse_indices = np.arange(
                    phase, full_times.size, step, dtype=int
                )
                first = int(sparse_indices[0])
                last = int(sparse_indices[-1])
                reference_times = full_times[first : last + 1]
                reference_values = full_values[first : last + 1]
                reference_occupation, reference_deficit, _phase_slope = (
                    baseline._reference_metrics(
                        reference_times,
                        reference_values,
                        baseline.THRESHOLD_MGL,
                    )
                )
                sparse_times = full_times[sparse_indices]
                sparse_values = full_values[sparse_indices]
                horizon = float(reference_times[-1] - reference_times[0])
                case_started = time.perf_counter()
                try:
                    result = identify_hypoxia_burden_with_replacements(
                        sparse_times,
                        sparse_values,
                        threshold=baseline.THRESHOLD_MGL,
                        max_slope=max_slope,
                        max_replacements=MAX_REPLACEMENTS,
                        state_lower_bound=baseline.STATE_LOWER_BOUND_MGL,
                    )
                except InsufficientReplacementBudget as error:
                    rows.append(
                        {
                            "dataset_id": key[0],
                            "block_index": key[1],
                            "cadence_hours": key[2],
                            "phase_index": key[3],
                            "status": "infeasible",
                            "error": str(error),
                            "runtime_seconds": time.perf_counter()
                            - case_started,
                        }
                    )
                    continue

                if result.minimum_replacements_for_feasibility != 0:
                    raise AssertionError(
                        "a frozen k=0-feasible case no longer has k_min=0"
                    )
                if sparse_indices.size <= 2 * MAX_REPLACEMENTS:
                    raise AssertionError("a case violates the M>2k condition")
                occupation_width_k1 = (
                    result.occupation_upper - result.occupation_lower
                ) / horizon
                deficit_width_k1 = (
                    result.oxygen_deficit_upper
                    - result.oxygen_deficit_lower
                ) / (baseline.THRESHOLD_MGL * horizon)
                if not (
                    -TOLERANCE <= occupation_width_k1 <= 1.0 + TOLERANCE
                ):
                    raise AssertionError(
                        "normalized occupation width violates its physical range"
                    )
                if not (
                    -TOLERANCE <= deficit_width_k1 <= 1.0 + TOLERANCE
                ):
                    raise AssertionError(
                        "normalized nonnegative deficit width violates its "
                        "physical range"
                    )
                occupation_width_k0 = float(
                    k0["occupation_width_fraction"]
                )
                deficit_width_k0 = float(
                    k0[
                        "deficit_width_fraction_nonnegative_H_times_T"
                    ]
                )
                deleted_counts = {
                    name: _deleted_count(result, name)
                    for name in (
                        "occupation_lower_witness",
                        "occupation_upper_witness",
                        "oxygen_deficit_lower_witness",
                        "oxygen_deficit_upper_witness",
                    )
                }
                if any(count > MAX_REPLACEMENTS for count in deleted_counts.values()):
                    raise AssertionError("an endpoint witness exceeded the fixed budget")

                rows.append(
                    {
                        "dataset_id": key[0],
                        "block_index": key[1],
                        "block_start": block["block_start"],
                        "block_end": block["block_end"],
                        "cadence_hours": key[2],
                        "phase_index": key[3],
                        "status": "feasible",
                        "error": "",
                        "reference_L_compatible": (
                            block_reference_max_slope <= max_slope + 1e-12
                        ),
                        "n_sparse_reports": int(sparse_indices.size),
                        "permitted_replacement_fraction": (
                            MAX_REPLACEMENTS / int(sparse_indices.size)
                        ),
                        "minimum_replacements_for_feasibility": (
                            result.minimum_replacements_for_feasibility
                        ),
                        "common_anchor_condition_M_gt_2k": (
                            int(sparse_indices.size)
                            > 2 * MAX_REPLACEMENTS
                        ),
                        "horizon_hours": horizon,
                        "threshold_mgL": baseline.THRESHOLD_MGL,
                        "state_lower_bound_mgL": (
                            baseline.STATE_LOWER_BOUND_MGL
                        ),
                        "max_slope_mgL_per_hour": max_slope,
                        "max_replacements": MAX_REPLACEMENTS,
                        "reference_occupation_hours": reference_occupation,
                        "reference_deficit_mgL_hours": reference_deficit,
                        "occupation_lower_k0": float(k0["occupation_lower"]),
                        "occupation_upper_k0": float(k0["occupation_upper"]),
                        "occupation_lower_k1": result.occupation_lower,
                        "occupation_upper_k1": result.occupation_upper,
                        "occupation_width_fraction_k0": occupation_width_k0,
                        "occupation_width_fraction_k1": occupation_width_k1,
                        "occupation_width_fraction_increase": (
                            occupation_width_k1 - occupation_width_k0
                        ),
                        "deficit_lower_k0_mgL_hours": float(
                            k0["deficit_lower_mgL_hours"]
                        ),
                        "deficit_upper_k0_mgL_hours": float(
                            k0["deficit_upper_nonnegative_mgL_hours"]
                        ),
                        "deficit_lower_k1_mgL_hours": (
                            result.oxygen_deficit_lower
                        ),
                        "deficit_upper_k1_mgL_hours": (
                            result.oxygen_deficit_upper
                        ),
                        "deficit_width_fraction_k0": deficit_width_k0,
                        "deficit_width_fraction_k1": deficit_width_k1,
                        "deficit_width_fraction_increase": (
                            deficit_width_k1 - deficit_width_k0
                        ),
                        "occupation_reference_covered_k1": (
                            result.occupation_lower - TOLERANCE
                            <= reference_occupation
                            <= result.occupation_upper + TOLERANCE
                        ),
                        "deficit_reference_covered_k1": (
                            result.oxygen_deficit_lower - TOLERANCE
                            <= reference_deficit
                            <= result.oxygen_deficit_upper + TOLERANCE
                        ),
                        "k0_nested_in_k1": (
                            result.occupation_lower
                            <= float(k0["occupation_lower"]) + TOLERANCE
                            and result.occupation_upper
                            >= float(k0["occupation_upper"]) - TOLERANCE
                            and result.oxygen_deficit_lower
                            <= float(k0["deficit_lower_mgL_hours"])
                            + TOLERANCE
                            and result.oxygen_deficit_upper
                            >= float(
                                k0[
                                    "deficit_upper_nonnegative_mgL_hours"
                                ]
                            )
                            - TOLERANCE
                        ),
                        "any_endpoint_uses_deletion": any(
                            count > 0 for count in deleted_counts.values()
                        ),
                        **{
                            f"{name}_deleted_count": count
                            for name, count in deleted_counts.items()
                        },
                        **{
                            f"{name}_deleted_indices": _deleted_indices(
                                result, name
                            )
                            for name in deleted_counts
                        },
                        "runtime_seconds": time.perf_counter() - case_started,
                    }
                )

    detail = pd.DataFrame(rows)
    feasible = detail.loc[detail["status"].eq("feasible")].copy()
    if len(detail) != EXPECTED_CASES or len(feasible) != EXPECTED_CASES:
        raise RuntimeError(
            "fixed-k=1 run did not return all frozen compatible cases"
        )
    if not feasible["reference_L_compatible"].all():
        raise RuntimeError("a frozen selected case lost reference-L compatibility")

    summary_rows: list[dict[str, object]] = []
    for cadence, group in feasible.groupby("cadence_hours", sort=True):
        summary_rows.append(
            {
                "cadence_hours": float(cadence),
                "cases": int(len(group)),
                "min_sparse_reports": int(group["n_sparse_reports"].min()),
                "max_sparse_reports": int(group["n_sparse_reports"].max()),
                "min_permitted_replacement_fraction": float(
                    group["permitted_replacement_fraction"].min()
                ),
                "max_permitted_replacement_fraction": float(
                    group["permitted_replacement_fraction"].max()
                ),
                "median_occupation_width_fraction_k0": float(
                    group["occupation_width_fraction_k0"].median()
                ),
                "median_occupation_width_fraction_k1": float(
                    group["occupation_width_fraction_k1"].median()
                ),
                "median_occupation_width_fraction_increase": float(
                    group["occupation_width_fraction_increase"].median()
                ),
                "median_deficit_width_fraction_k0": float(
                    group["deficit_width_fraction_k0"].median()
                ),
                "median_deficit_width_fraction_k1": float(
                    group["deficit_width_fraction_k1"].median()
                ),
                "median_deficit_width_fraction_increase": float(
                    group["deficit_width_fraction_increase"].median()
                ),
                "p90_occupation_width_fraction_k1": float(
                    group["occupation_width_fraction_k1"].quantile(0.9)
                ),
                "p90_deficit_width_fraction_k1": float(
                    group["deficit_width_fraction_k1"].quantile(0.9)
                ),
                "max_occupation_width_fraction_k1": float(
                    group["occupation_width_fraction_k1"].max()
                ),
                "max_deficit_width_fraction_k1": float(
                    group["deficit_width_fraction_k1"].max()
                ),
                "cases_occupation_width_increased": int(
                    (
                        group["occupation_width_fraction_increase"]
                        > TOLERANCE
                    ).sum()
                ),
                "cases_deficit_width_increased": int(
                    (
                        group["deficit_width_fraction_increase"]
                        > TOLERANCE
                    ).sum()
                ),
                "cases_any_endpoint_uses_deletion": int(
                    group["any_endpoint_uses_deletion"].sum()
                ),
                "occupation_containment_failures_k1": int(
                    (~group["occupation_reference_covered_k1"]).sum()
                ),
                "deficit_containment_failures_k1": int(
                    (~group["deficit_reference_covered_k1"]).sum()
                ),
                "nesting_failures": int((~group["k0_nested_in_k1"]).sum()),
                "median_runtime_seconds": float(
                    group["runtime_seconds"].median()
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)

    metadata: dict[str, object] = {
        "analysis": "fixed_k1_sensitivity_on_frozen_compatible_7day_cases",
        "interpretation": (
            "smallest positive fixed replacement budget; not an estimated "
            "corruption count and not k_min(L)"
        ),
        "analysis_freeze_note": str(FREEZE_NOTE),
        "analysis_freeze_note_sha256": _sha256(FREEZE_NOTE),
        "max_replacements": MAX_REPLACEMENTS,
        "threshold_mgL": baseline.THRESHOLD_MGL,
        "state_lower_bound_mgL": baseline.STATE_LOWER_BOUND_MGL,
        "calibration_q999_max_slope_mgL_per_hour": max_slope,
        "case_selection": (
            "all and only frozen k=0 feasible reference-L-compatible cases"
        ),
        "expected_case_count": EXPECTED_CASES,
        "case_count": int(len(detail)),
        "feasible_case_count": int(len(feasible)),
        "compatible_block_count": int(
            feasible[["dataset_id", "block_index"]]
            .drop_duplicates()
            .shape[0]
        ),
        "compatible_station_count": int(feasible["dataset_id"].nunique()),
        "all_minimum_replacements_equal_zero": bool(
            feasible["minimum_replacements_for_feasibility"].eq(0).all()
        ),
        "all_cases_satisfy_M_gt_2k": bool(
            feasible["common_anchor_condition_M_gt_2k"].all()
        ),
        "all_physical_width_checks_passed": bool(
            feasible["occupation_width_fraction_k1"]
            .between(-TOLERANCE, 1.0 + TOLERANCE)
            .all()
            and feasible["deficit_width_fraction_k1"]
            .between(-TOLERANCE, 1.0 + TOLERANCE)
            .all()
        ),
        "source_baseline_detail": str(BASELINE_DETAIL),
        "source_baseline_detail_sha256": _sha256(BASELINE_DETAIL),
        "source_manifest": str(baseline.MANIFEST),
        "source_manifest_sha256": _sha256(baseline.MANIFEST),
        "source_v3_metadata": str(baseline.V3_METADATA),
        "source_v3_metadata_sha256": _sha256(baseline.V3_METADATA),
        "core_runtime": str(CORE_RUNTIME),
        "core_runtime_sha256": _sha256(CORE_RUNTIME),
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": _sha256(Path(__file__).resolve()),
        "total_runtime_seconds": time.perf_counter() - started,
        "containment_failures": int(
            (~feasible["occupation_reference_covered_k1"]).sum()
            + (~feasible["deficit_reference_covered_k1"]).sum()
        ),
        "nesting_failures": int((~feasible["k0_nested_in_k1"]).sum()),
        "summary": json.loads(summary.to_json(orient="records")),
    }
    return detail, summary, metadata


def write_outputs() -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    detail, summary, metadata = run()
    detail_path = OUTPUT / "fixed_k1_detail.csv"
    summary_path = OUTPUT / "fixed_k1_summary.csv"
    metadata_path = OUTPUT / "fixed_k1_metadata.json"
    detail.to_csv(detail_path, index=False)
    summary.to_csv(summary_path, index=False)
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# Fixed-k=1 seven-day sensitivity report",
        "",
        "This analysis reuses exactly the frozen reference-compatible Table 2 "
        "cases and changes only the maximum replacement budget from k=0 to "
        "the fixed value k=1, which permits at most one arbitrary replacement. "
        "It is a sensitivity scenario, not an estimated error prevalence.",
        "",
        "| Cadence h | Cases | Median occupation width k=0 | k=1 | "
        "Median within-case increase | Median deficit width k=0 | k=1 | "
        "Median within-case increase | Cases with an endpoint witness using "
        "a permitted deletion | Failures |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        failures = (
            int(row.occupation_containment_failures_k1)
            + int(row.deficit_containment_failures_k1)
            + int(row.nesting_failures)
        )
        lines.append(
            f"| {row.cadence_hours:g} | {int(row.cases)} | "
            f"{row.median_occupation_width_fraction_k0:.3f} | "
            f"{row.median_occupation_width_fraction_k1:.3f} | "
            f"{row.median_occupation_width_fraction_increase:.3f} | "
            f"{row.median_deficit_width_fraction_k0:.3f} | "
            f"{row.median_deficit_width_fraction_k1:.3f} | "
            f"{row.median_deficit_width_fraction_increase:.3f} | "
            f"{int(row.cases_any_endpoint_uses_deletion)} | {failures} |"
        )
    lines.extend(
        [
            "",
            f"Total runtime: {metadata['total_runtime_seconds']:.2f} s.",
            "All widths are normalized by T for occupation and HT for deficit.",
        ]
    )
    report = OUTPUT / "FIXED_K1_SENSITIVITY_REPORT.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    print(write_outputs())
