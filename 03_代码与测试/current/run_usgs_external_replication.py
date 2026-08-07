# -*- coding: utf-8 -*-
"""Run the frozen USGS cross-network dissolved-oxygen stress test.

This adapter does not recalibrate the Maryland scenario labels.  It applies
the previously frozen H, B, L, cadence, phase, seven-day window and k rules to
all nine metadata-eligible USGS Maryland series.  Main summaries use the
station-block as the descriptive unit; overlapping cadence phases are retained
only as within-block numerical scenarios.
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
DATA = PROJECT / "04_真实数据" / "USGS_WDFN_external_2024"
MANIFEST = DATA / "download_manifest.json"
FREEZE_NOTE = (
    PROJECT
    / "09_论文写作输出"
    / "USGS_EXTERNAL_REPLICATION_FREEZE_NOTE_20260802.md"
)
OUTPUT = PROJECT / "05_既有结果" / "usgs_external_7day_replication"
CORE_RUNTIME = HERE / "robust_exact_current.py"

EXPECTED_SERIES = (
    ("USGS-01491000", "79480fea74964d81bf3822a6a66ee1c7"),
    ("USGS-01579550", "abaf8aaaf7bd495483198ad33d277024"),
    ("USGS-01594441", "979a9c0c3a4842e8bd56bd05b649ab28"),
    ("USGS-01638500", "7329a5d59c164090aee93904299d4a7a"),
    ("USGS-01643580", "4d259ae612cb4d82bd36c7ccd1051434"),
    ("USGS-01646500", "bb18f3558934465194f8da1d5f1c3df3"),
    ("USGS-01649190", "42075f9554bd4ee0b1909d52959482ae"),
    ("USGS-01649500", "80262e6c200045eb8755f9f1780b799b"),
    ("USGS-01650800", "d4a0194d99aa4663b9c935003d7ed3f0"),
)

ANALYSIS_START = pd.Timestamp("2024-06-01T00:00:00Z")
ANALYSIS_END = pd.Timestamp("2024-09-01T00:00:00Z")
MAX_SLOPE_MGL_PER_HOUR = 13.203959999999936
MAX_REPLACEMENTS_SENSITIVITY = 1
TOLERANCE = 1e-8


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _deleted_count(result: object, name: str) -> int:
    return len(getattr(result, name).deleted_indices)


def _deleted_indices(result: object, name: str) -> str:
    return json.dumps(
        list(getattr(result, name).deleted_indices), separators=(",", ":")
    )


def _load_inputs() -> tuple[list[dict[str, object]], pd.DataFrame]:
    if not FREEZE_NOTE.exists():
        raise FileNotFoundError("pre-computation freeze note is missing")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    observed = tuple(
        (
            str(item["monitoring_location_id"]),
            str(item["time_series_id"]),
        )
        for item in manifest["series"]
    )
    if observed != EXPECTED_SERIES:
        raise ValueError("download manifest differs from the frozen series list")
    if manifest["parameter_code"] != "00300" or manifest["unit"] != "mg/l":
        raise ValueError("unexpected parameter identity or unit")
    if tuple(manifest["frozen_analysis_interval_utc_half_open"]) != (
        "2024-06-01T00:00:00Z",
        "2024-09-01T00:00:00Z",
    ):
        raise ValueError("download manifest differs from the frozen interval")

    blocks: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    for item in manifest["series"]:
        location_id = str(item["monitoring_location_id"])
        source_path = DATA / str(item["normalized_file"])
        if _sha256(source_path) != str(item["normalized_sha256"]):
            raise ValueError(f"{location_id}: normalized source hash mismatch")
        frame = pd.read_csv(source_path, dtype=str, keep_default_na=False)
        required = {
            "monitoring_location_id",
            "time_series_id",
            "site_name",
            "agency_code",
            "parameter_code",
            "unit",
            "statistic_code",
            "time_utc",
            "value_mgL",
            "qualifiers",
        }
        if not required.issubset(frame.columns):
            raise ValueError(f"{location_id}: normalized columns are incomplete")
        identity_checks = {
            "monitoring_location_id": {location_id},
            "time_series_id": {str(item["time_series_id"])},
            "agency_code": {"USGS"},
            "parameter_code": {"00300"},
            "unit": {"mg/l"},
            "statistic_code": {"00000"},
        }
        for column, expected in identity_checks.items():
            observed_values = set(frame[column].astype(str).unique())
            if observed_values != expected:
                raise ValueError(
                    f"{location_id}: unexpected {column}: {observed_values}"
                )

        frame["time"] = pd.to_datetime(frame["time_utc"], utc=True)
        frame["do_mgL"] = pd.to_numeric(frame["value_mgL"], errors="coerce")
        in_interval = frame["time"].ge(ANALYSIS_START) & frame["time"].lt(
            ANALYSIS_END
        )
        interval_frame = frame.loc[in_interval].copy()
        approved_only = interval_frame["qualifiers"].astype(str).eq("A")
        finite = np.isfinite(interval_frame["do_mgL"].to_numpy(dtype=float))
        analysis_candidates = interval_frame.loc[approved_only & finite].copy()
        on_quarter_hour = (
            analysis_candidates["time"].dt.minute.isin([0, 15, 30, 45])
            & analysis_candidates["time"].dt.second.eq(0)
            & analysis_candidates["time"].dt.microsecond.eq(0)
        )
        reference = analysis_candidates.loc[on_quarter_hour, ["time", "do_mgL"]]
        reference = reference.sort_values("time", kind="stable").reset_index(
            drop=True
        )
        duplicate_count = int(reference["time"].duplicated(keep=False).sum())
        if duplicate_count:
            raise ValueError(
                f"{location_id}: duplicate UTC timestamps in reference series"
            )
        if (reference["do_mgL"] < 0.0).any():
            raise ValueError(
                f"{location_id}: an approved value violates the frozen state floor"
            )

        raw_times = interval_frame["time"].sort_values(kind="stable")
        raw_gaps = raw_times.diff().dt.total_seconds().div(60.0).dropna()
        reference_gaps = (
            reference["time"].diff().dt.total_seconds().div(60.0).dropna()
        )
        station_blocks = baseline._find_blocks(reference)
        qualifier_counts = {
            str(key): int(value)
            for key, value in interval_frame["qualifiers"]
            .astype(str)
            .value_counts(dropna=False)
            .sort_index()
            .items()
        }
        audit_rows.append(
            {
                "monitoring_location_id": location_id,
                "time_series_id": str(item["time_series_id"]),
                "site_name": str(item["site_name"]),
                "latitude": item.get("latitude"),
                "longitude": item.get("longitude"),
                "huc_code": item.get("huc_code"),
                "raw_response_rows": int(item["row_count"]),
                "rows_in_frozen_interval": int(len(interval_frame)),
                "exact_A_rows_in_interval": int(approved_only.sum()),
                "nonfinite_rows_excluded": int((~finite).sum()),
                "other_qualifier_rows_excluded": int((~approved_only).sum()),
                "quarter_hour_reference_rows": int(len(reference)),
                "qualifying_7day_blocks": int(len(station_blocks)),
                "duplicate_reference_timestamps": duplicate_count,
                "native_median_gap_minutes": (
                    float(raw_gaps.median()) if len(raw_gaps) else math.nan
                ),
                "native_modal_gap_minutes": (
                    float(raw_gaps.mode().iloc[0]) if len(raw_gaps) else math.nan
                ),
                "reference_median_gap_minutes": (
                    float(reference_gaps.median())
                    if len(reference_gaps)
                    else math.nan
                ),
                "minimum_DO_mgL_in_reference": (
                    float(reference["do_mgL"].min())
                    if len(reference)
                    else math.nan
                ),
                "maximum_DO_mgL_in_reference": (
                    float(reference["do_mgL"].max())
                    if len(reference)
                    else math.nan
                ),
                "qualifier_counts_json": json.dumps(
                    qualifier_counts, sort_keys=True, separators=(",", ":")
                ),
                "source_query_url": str(item["query_url"]),
                "raw_sha256": str(item["raw_sha256"]),
                "normalized_sha256": str(item["normalized_sha256"]),
            }
        )
        for block_index, block in enumerate(station_blocks):
            blocks.append(
                {
                    "monitoring_location_id": location_id,
                    "time_series_id": str(item["time_series_id"]),
                    "site_name": str(item["site_name"]),
                    "block_index": block_index,
                    "block_start": str(block["time"].iloc[0]),
                    "block_end": str(block["time"].iloc[-1]),
                    "times": block["elapsed_hours"].to_numpy(dtype=float),
                    "values": block["do_mgL"].to_numpy(dtype=float),
                }
            )
    return blocks, pd.DataFrame(audit_rows)


def run() -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, object],
]:
    blocks, data_audit = _load_inputs()
    k0_rows: list[dict[str, object]] = []
    k1_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    started = time.perf_counter()

    for number, block in enumerate(blocks, start=1):
        full_times = np.asarray(block["times"], dtype=float)
        full_values = np.asarray(block["values"], dtype=float)
        full_occ, full_def, block_max_slope = baseline._reference_metrics(
            full_times, full_values, baseline.THRESHOLD_MGL
        )
        compatible = bool(
            block_max_slope <= MAX_SLOPE_MGL_PER_HOUR + TOLERANCE
        )
        block_rows.append(
            {
                "monitoring_location_id": block["monitoring_location_id"],
                "time_series_id": block["time_series_id"],
                "site_name": block["site_name"],
                "block_index": block["block_index"],
                "block_start": block["block_start"],
                "block_end": block["block_end"],
                "reference_points": int(full_times.size),
                "minimum_DO_mgL": float(full_values.min()),
                "maximum_DO_mgL": float(full_values.max()),
                "reference_occupation_hours": full_occ,
                "reference_deficit_mgL_hours": full_def,
                "reference_max_slope_mgL_per_hour": block_max_slope,
                "reference_L_compatible": compatible,
            }
        )
        print(
            f"block {number}/{len(blocks)} "
            f"{block['monitoring_location_id']}:{block['block_index']} "
            f"compatible={compatible}",
            flush=True,
        )

        for cadence in baseline.CADENCE_HOURS:
            step = int(round(cadence / baseline.BASE_CADENCE_HOURS))
            for phase in range(step):
                sparse_indices = np.arange(phase, full_times.size, step, dtype=int)
                if sparse_indices.size < 2:
                    continue
                first = int(sparse_indices[0])
                last = int(sparse_indices[-1])
                reference_times = full_times[first : last + 1]
                reference_values = full_values[first : last + 1]
                reference_occ, reference_def, _ = baseline._reference_metrics(
                    reference_times,
                    reference_values,
                    baseline.THRESHOLD_MGL,
                )
                sparse_times = full_times[sparse_indices]
                sparse_values = full_values[sparse_indices]
                horizon = float(reference_times[-1] - reference_times[0])
                base = {
                    "monitoring_location_id": block["monitoring_location_id"],
                    "time_series_id": block["time_series_id"],
                    "site_name": block["site_name"],
                    "block_index": block["block_index"],
                    "block_start": block["block_start"],
                    "block_end": block["block_end"],
                    "cadence_hours": cadence,
                    "phase_index": phase,
                    "n_sparse_reports": int(sparse_indices.size),
                    "horizon_hours": horizon,
                    "threshold_mgL": baseline.THRESHOLD_MGL,
                    "state_lower_bound_mgL": baseline.STATE_LOWER_BOUND_MGL,
                    "max_slope_mgL_per_hour": MAX_SLOPE_MGL_PER_HOUR,
                    "reference_L_compatible": compatible,
                    "reference_occupation_hours": reference_occ,
                    "reference_deficit_mgL_hours": reference_def,
                }
                case_started = time.perf_counter()
                try:
                    unconstrained = identify_hypoxia_burden_with_replacements(
                        sparse_times,
                        sparse_values,
                        threshold=baseline.THRESHOLD_MGL,
                        max_slope=MAX_SLOPE_MGL_PER_HOUR,
                        max_replacements=0,
                    )
                    constrained = identify_hypoxia_burden_with_replacements(
                        sparse_times,
                        sparse_values,
                        threshold=baseline.THRESHOLD_MGL,
                        max_slope=MAX_SLOPE_MGL_PER_HOUR,
                        max_replacements=0,
                        state_lower_bound=baseline.STATE_LOWER_BOUND_MGL,
                    )
                except InsufficientReplacementBudget as error:
                    k0_rows.append(
                        {
                            **base,
                            "status": "infeasible",
                            "error": str(error),
                            "runtime_seconds": time.perf_counter() - case_started,
                        }
                    )
                    continue

                occupation_width = (
                    constrained.occupation_upper - constrained.occupation_lower
                ) / horizon
                deficit_width_unconstrained = (
                    unconstrained.oxygen_deficit_upper
                    - unconstrained.oxygen_deficit_lower
                ) / (baseline.THRESHOLD_MGL * horizon)
                deficit_width = (
                    constrained.oxygen_deficit_upper
                    - constrained.oxygen_deficit_lower
                ) / (baseline.THRESHOLD_MGL * horizon)
                k0_row = {
                    **base,
                    "status": "feasible",
                    "error": "",
                    "occupation_lower": constrained.occupation_lower,
                    "occupation_upper": constrained.occupation_upper,
                    "occupation_width_fraction": occupation_width,
                    "deficit_lower_mgL_hours": constrained.oxygen_deficit_lower,
                    "deficit_upper_unconstrained_mgL_hours": (
                        unconstrained.oxygen_deficit_upper
                    ),
                    "deficit_upper_nonnegative_mgL_hours": (
                        constrained.oxygen_deficit_upper
                    ),
                    "deficit_width_fraction_unconstrained_H_times_T": (
                        deficit_width_unconstrained
                    ),
                    "deficit_width_fraction_nonnegative_H_times_T": deficit_width,
                    "deficit_width_fraction_reduction": (
                        deficit_width_unconstrained - deficit_width
                    ),
                    "occupation_reference_covered": (
                        constrained.occupation_lower - TOLERANCE
                        <= reference_occ
                        <= constrained.occupation_upper + TOLERANCE
                    ),
                    "deficit_reference_covered": (
                        constrained.oxygen_deficit_lower - TOLERANCE
                        <= reference_def
                        <= constrained.oxygen_deficit_upper + TOLERANCE
                    ),
                    "runtime_seconds": time.perf_counter() - case_started,
                }
                k0_rows.append(k0_row)

                if not compatible:
                    continue
                k1_started = time.perf_counter()
                result = identify_hypoxia_burden_with_replacements(
                    sparse_times,
                    sparse_values,
                    threshold=baseline.THRESHOLD_MGL,
                    max_slope=MAX_SLOPE_MGL_PER_HOUR,
                    max_replacements=MAX_REPLACEMENTS_SENSITIVITY,
                    state_lower_bound=baseline.STATE_LOWER_BOUND_MGL,
                )
                if result.minimum_replacements_for_feasibility != 0:
                    raise AssertionError("a k=0-feasible case lost k_min=0")
                if sparse_indices.size <= 2 * MAX_REPLACEMENTS_SENSITIVITY:
                    raise AssertionError("M > 2k failed")
                occ_width_k1 = (
                    result.occupation_upper - result.occupation_lower
                ) / horizon
                def_width_k1 = (
                    result.oxygen_deficit_upper - result.oxygen_deficit_lower
                ) / (baseline.THRESHOLD_MGL * horizon)
                witness_names = (
                    "occupation_lower_witness",
                    "occupation_upper_witness",
                    "oxygen_deficit_lower_witness",
                    "oxygen_deficit_upper_witness",
                )
                deleted_counts = {
                    name: _deleted_count(result, name) for name in witness_names
                }
                k1_rows.append(
                    {
                        **base,
                        "status": "feasible",
                        "error": "",
                        "max_replacements": MAX_REPLACEMENTS_SENSITIVITY,
                        "minimum_replacements_for_feasibility": (
                            result.minimum_replacements_for_feasibility
                        ),
                        "common_anchor_condition_M_gt_2k": bool(
                            sparse_indices.size
                            > 2 * MAX_REPLACEMENTS_SENSITIVITY
                        ),
                        "occupation_lower_k0": constrained.occupation_lower,
                        "occupation_upper_k0": constrained.occupation_upper,
                        "occupation_lower_k1": result.occupation_lower,
                        "occupation_upper_k1": result.occupation_upper,
                        "occupation_width_fraction_k0": occupation_width,
                        "occupation_width_fraction_k1": occ_width_k1,
                        "occupation_width_fraction_increase": (
                            occ_width_k1 - occupation_width
                        ),
                        "deficit_lower_k0_mgL_hours": (
                            constrained.oxygen_deficit_lower
                        ),
                        "deficit_upper_k0_mgL_hours": (
                            constrained.oxygen_deficit_upper
                        ),
                        "deficit_lower_k1_mgL_hours": result.oxygen_deficit_lower,
                        "deficit_upper_k1_mgL_hours": result.oxygen_deficit_upper,
                        "deficit_width_fraction_k0": deficit_width,
                        "deficit_width_fraction_k1": def_width_k1,
                        "deficit_width_fraction_increase": (
                            def_width_k1 - deficit_width
                        ),
                        "occupation_reference_covered_k1": (
                            result.occupation_lower - TOLERANCE
                            <= reference_occ
                            <= result.occupation_upper + TOLERANCE
                        ),
                        "deficit_reference_covered_k1": (
                            result.oxygen_deficit_lower - TOLERANCE
                            <= reference_def
                            <= result.oxygen_deficit_upper + TOLERANCE
                        ),
                        "k0_nested_in_k1": bool(
                            result.occupation_lower
                            <= constrained.occupation_lower + TOLERANCE
                            and result.occupation_upper
                            >= constrained.occupation_upper - TOLERANCE
                            and result.oxygen_deficit_lower
                            <= constrained.oxygen_deficit_lower + TOLERANCE
                            and result.oxygen_deficit_upper
                            >= constrained.oxygen_deficit_upper - TOLERANCE
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
                            for name in witness_names
                        },
                        "runtime_seconds": time.perf_counter() - k1_started,
                    }
                )

    k0_detail = pd.DataFrame(k0_rows)
    k1_detail = pd.DataFrame(k1_rows)
    station_blocks = pd.DataFrame(block_rows)
    compatible_k0 = k0_detail.loc[
        k0_detail["status"].eq("feasible")
        & k0_detail["reference_L_compatible"].eq(True)
    ].copy()
    block_cadence = (
        compatible_k0.groupby(
            [
                "monitoring_location_id",
                "time_series_id",
                "site_name",
                "block_index",
                "block_start",
                "block_end",
                "cadence_hours",
            ],
            as_index=False,
        )
        .agg(
            phase_cases=("phase_index", "size"),
            median_occupation_width_fraction=(
                "occupation_width_fraction",
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
            occupation_containment_failures=(
                "occupation_reference_covered",
                lambda values: int((~values).sum()),
            ),
            deficit_containment_failures=(
                "deficit_reference_covered",
                lambda values: int((~values).sum()),
            ),
        )
    )
    if len(k1_detail):
        k1_block = (
            k1_detail.groupby(
                ["monitoring_location_id", "block_index", "cadence_hours"],
                as_index=False,
            )
            .agg(
                median_occupation_width_fraction_k1=(
                    "occupation_width_fraction_k1",
                    "median",
                ),
                median_occupation_width_fraction_increase=(
                    "occupation_width_fraction_increase",
                    "median",
                ),
                median_deficit_width_fraction_k1=(
                    "deficit_width_fraction_k1",
                    "median",
                ),
                median_deficit_width_fraction_increase=(
                    "deficit_width_fraction_increase",
                    "median",
                ),
                k1_nesting_failures=(
                    "k0_nested_in_k1",
                    lambda values: int((~values).sum()),
                ),
            )
        )
        block_cadence = block_cadence.merge(
            k1_block,
            on=["monitoring_location_id", "block_index", "cadence_hours"],
            how="left",
            validate="one_to_one",
        )

    cadence_summary = (
        block_cadence.groupby("cadence_hours", as_index=False)
        .agg(
            compatible_station_blocks=("block_index", "size"),
            compatible_stations=("monitoring_location_id", "nunique"),
            within_block_phase_cases=("phase_cases", "sum"),
            median_block_occupation_width_fraction=(
                "median_occupation_width_fraction",
                "median",
            ),
            q25_block_occupation_width_fraction=(
                "median_occupation_width_fraction",
                lambda values: float(values.quantile(0.25)),
            ),
            q75_block_occupation_width_fraction=(
                "median_occupation_width_fraction",
                lambda values: float(values.quantile(0.75)),
            ),
            median_block_deficit_width_fraction=(
                "median_deficit_width_fraction_nonnegative",
                "median",
            ),
            q25_block_deficit_width_fraction=(
                "median_deficit_width_fraction_nonnegative",
                lambda values: float(values.quantile(0.25)),
            ),
            q75_block_deficit_width_fraction=(
                "median_deficit_width_fraction_nonnegative",
                lambda values: float(values.quantile(0.75)),
            ),
            median_block_occupation_k1_increase=(
                "median_occupation_width_fraction_increase",
                "median",
            ),
            median_block_deficit_k1_increase=(
                "median_deficit_width_fraction_increase",
                "median",
            ),
            containment_failures=(
                "occupation_containment_failures",
                lambda values: int(values.sum()),
            ),
            deficit_containment_failures=(
                "deficit_containment_failures",
                lambda values: int(values.sum()),
            ),
            k1_nesting_failures=(
                "k1_nesting_failures",
                lambda values: int(values.sum()),
            ),
        )
    )

    compatible_blocks = station_blocks.loc[
        station_blocks["reference_L_compatible"].eq(True)
    ]
    metadata = {
        "status": "complete",
        "interpretation": (
            "cross-network transportability stress test; not external "
            "validation or population inference"
        ),
        "selected_series_count": len(EXPECTED_SERIES),
        "series_with_qualifying_blocks": int(
            data_audit["qualifying_7day_blocks"].gt(0).sum()
        ),
        "qualifying_block_count": int(len(station_blocks)),
        "compatible_block_count": int(len(compatible_blocks)),
        "compatible_station_count": int(
            compatible_blocks["monitoring_location_id"].nunique()
        ),
        "all_phase_case_count": int(len(k0_detail)),
        "compatible_feasible_phase_case_count": int(len(compatible_k0)),
        "fixed_k1_case_count": int(len(k1_detail)),
        "threshold_mgL": baseline.THRESHOLD_MGL,
        "state_lower_bound_mgL": baseline.STATE_LOWER_BOUND_MGL,
        "transferred_calibration_q999_L_mgL_per_hour": (
            MAX_SLOPE_MGL_PER_HOUR
        ),
        "cadence_hours_all_phases": list(baseline.CADENCE_HOURS),
        "max_replacements_main": 0,
        "max_replacements_sensitivity": MAX_REPLACEMENTS_SENSITIVITY,
        "main_summary_unit": (
            "median across station-block phase medians; phases are not "
            "independent replicates"
        ),
        "analysis_interval_utc_half_open": [
            str(ANALYSIS_START),
            str(ANALYSIS_END),
        ],
        "approved_rule": (
            "legacy NWIS qualifier exactly A; other/nonfinite observations "
            "retained in source files but excluded before window selection"
        ),
        "quarter_hour_rule": (
            "exact UTC minutes 00/15/30/45, second and microsecond zero; "
            "no interpolation or averaging"
        ),
        "window_rule": {
            "points": baseline.WINDOW_POINTS,
            "span_hours_inclusive": [
                baseline.MIN_WINDOW_HOURS,
                baseline.MAX_WINDOW_HOURS,
            ],
            "maximum_adjacent_gap_minutes": baseline.MAX_GAP_HOURS * 60.0,
            "maximum_blocks_per_station": baseline.MAX_BLOCKS_PER_STATION,
            "selection": "earliest qualifying row-disjoint windows",
        },
        "containment_failures": int(
            (~compatible_k0["occupation_reference_covered"]).sum()
            + (~compatible_k0["deficit_reference_covered"]).sum()
            + (
                (~k1_detail["occupation_reference_covered_k1"]).sum()
                if len(k1_detail)
                else 0
            )
            + (
                (~k1_detail["deficit_reference_covered_k1"]).sum()
                if len(k1_detail)
                else 0
            )
        ),
        "k1_nesting_failures": int(
            (~k1_detail["k0_nested_in_k1"]).sum() if len(k1_detail) else 0
        ),
        "freeze_note": "09_论文写作输出/USGS_EXTERNAL_REPLICATION_FREEZE_NOTE_20260802.md",
        "freeze_note_sha256": _sha256(FREEZE_NOTE),
        "download_manifest": "04_真实数据/USGS_WDFN_external_2024/download_manifest.json",
        "download_manifest_sha256": _sha256(MANIFEST),
        "core_runtime": "03_代码与测试/current/robust_exact_current.py",
        "core_runtime_sha256": _sha256(CORE_RUNTIME),
        "analysis_script": "03_代码与测试/current/run_usgs_external_replication.py",
        "analysis_script_sha256": _sha256(Path(__file__).resolve()),
    }
    return (
        k0_detail,
        k1_detail,
        cadence_summary,
        block_cadence,
        station_blocks,
        data_audit,
        metadata,
    )


def write_outputs() -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (
        k0_detail,
        k1_detail,
        cadence_summary,
        block_cadence,
        station_blocks,
        data_audit,
        metadata,
    ) = run()
    k0_detail.drop(columns=["runtime_seconds"], errors="ignore").to_csv(
        OUTPUT / "usgs_external_k0_detail.csv", index=False
    )
    k1_detail.drop(columns=["runtime_seconds"], errors="ignore").to_csv(
        OUTPUT / "usgs_external_k1_detail.csv", index=False
    )
    cadence_summary.to_csv(
        OUTPUT / "usgs_external_cadence_summary.csv", index=False
    )
    block_cadence.to_csv(
        OUTPUT / "usgs_external_block_cadence_summary.csv", index=False
    )
    station_blocks.to_csv(
        OUTPUT / "usgs_external_station_block_summary.csv", index=False
    )
    data_audit.to_csv(OUTPUT / "usgs_external_data_audit.csv", index=False)
    (OUTPUT / "usgs_external_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# USGS external-network seven-day stress test",
        "",
        "All nine metadata-eligible Maryland USGS series were retained before "
        "endpoint computation. The Maryland-derived scenario labels were "
        "transferred unchanged. This is not external validation or population "
        "inference.",
        "",
        f"{metadata['series_with_qualifying_blocks']}/"
        f"{metadata['selected_series_count']} series supplied at least one "
        f"qualifying window; {metadata['compatible_block_count']}/"
        f"{metadata['qualifying_block_count']} station-blocks from "
        f"{metadata['compatible_station_count']} stations were compatible with "
        "the transferred L label.",
        "",
        "Main summaries below are medians across station-block phase medians. "
        "The phase cases are overlapping numerical scenarios, not independent "
        "replicates.",
        "",
        "| Cadence (h) | Compatible station-blocks | Stations | Within-block "
        "phase cases | Occupation width/T, median [IQR] | Deficit width/(HT), "
        "median [IQR] | Median k=1 increase: occupation; deficit | Failures |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in cadence_summary.itertuples(index=False):
        failures = int(
            row.containment_failures
            + row.deficit_containment_failures
            + row.k1_nesting_failures
        )
        lines.append(
            f"| {row.cadence_hours:g} | "
            f"{int(row.compatible_station_blocks)} | "
            f"{int(row.compatible_stations)} | "
            f"{int(row.within_block_phase_cases)} | "
            f"{row.median_block_occupation_width_fraction:.3f} "
            f"[{row.q25_block_occupation_width_fraction:.3f}, "
            f"{row.q75_block_occupation_width_fraction:.3f}] | "
            f"{row.median_block_deficit_width_fraction:.3f} "
            f"[{row.q25_block_deficit_width_fraction:.3f}, "
            f"{row.q75_block_deficit_width_fraction:.3f}] | "
            f"{row.median_block_occupation_k1_increase:.3f}; "
            f"{row.median_block_deficit_k1_increase:.3f} | {failures} |"
        )
    lines.extend(
        [
            "",
            "USGS source citation: U.S. Geological Survey (2026), USGS "
            "National Water Information System database, accessed August 2, "
            "2026, https://doi.org/10.5066/F7P55KJN.",
        ]
    )
    report = OUTPUT / "USGS_EXTERNAL_REPLICATION_REPORT.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    print(write_outputs())
