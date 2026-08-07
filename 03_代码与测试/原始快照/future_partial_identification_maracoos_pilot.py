# -*- coding: utf-8 -*-
"""Identity-aware masking pilot on eight independent MARACOOS/MDDNR station-years."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from future_partial_identification import (
    AffineSegment,
    InfeasibleObservationIntervals,
    identify_hypoxia_burden,
    integrate_occupation,
    integrate_oxygen_deficit,
)


PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / "data_external_future" / "MARACOOS_MDDNR_discovery_v3"
MANIFEST = DATA / "download_manifest.json"
OUTPUT = (
    PROJECT
    / "research_branch_outputs"
    / "future_agenda_20260716"
    / "partial_identification_maracoos_pilot_v3_collision_identity"
)

THRESHOLDS = (2.0, 5.0)
CADENCE_HOURS = (0.5, 1.0, 2.0, 4.0, 6.0)
BASE_CADENCE_HOURS = 0.25
MAX_SLOPE_GAP_HOURS = 0.30
WINDOW_POINTS = 30 * 24 * 4 + 1
WINDOW_HOURS = 30 * 24
MAX_BLOCKS_PER_STATION = 2
NOMINAL_CADENCE_SECONDS = 15 * 60
PARALLEL_COLLISION_SECONDS = NOMINAL_CADENCE_SECONDS / 2
MIN_PARALLEL_COLLISIONS = 10
MIN_MINORITY_COLLISION_SHARE = 0.05


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cross_group_collision_diagnostics(
    present_times: pd.Series,
    missing_times: pd.Series,
) -> tuple[int, int, float]:
    """Count minority rows with an opposite-group neighbor inside half a cadence."""

    # Force nanosecond resolution explicitly. Pandas may store timezone-aware
    # datetimes internally at microsecond resolution, so bare astype("int64")
    # is not guaranteed to have nanosecond units across runtimes.
    present_ns = np.sort(
        present_times.dt.tz_convert(None)
        .to_numpy(dtype="datetime64[ns]")
        .astype(np.int64)
    )
    missing_ns = np.sort(
        missing_times.dt.tz_convert(None)
        .to_numpy(dtype="datetime64[ns]")
        .astype(np.int64)
    )
    if present_ns.size == 0 or missing_ns.size == 0:
        return 0, int(min(present_ns.size, missing_ns.size)), 0.0
    if present_ns.size <= missing_ns.size:
        minority_ns, other_ns = present_ns, missing_ns
    else:
        minority_ns, other_ns = missing_ns, present_ns
    positions = np.searchsorted(other_ns, minority_ns, side="left")
    left_positions = np.clip(positions - 1, 0, other_ns.size - 1)
    right_positions = np.clip(positions, 0, other_ns.size - 1)
    nearest_ns = np.minimum(
        np.abs(minority_ns - other_ns[left_positions]),
        np.abs(minority_ns - other_ns[right_positions]),
    )
    collision_count = int(
        (nearest_ns < PARALLEL_COLLISION_SECONDS * 1_000_000_000).sum()
    )
    return collision_count, int(minority_ns.size), float(collision_count / minority_ns.size)


def load_station(item: dict[str, object]) -> tuple[pd.DataFrame, dict[str, object]]:
    path = DATA / str(item["csv_file"])
    actual_sha = sha256_file(path)
    if actual_sha != item["csv_sha256"]:
        raise ValueError(f"SHA-256 mismatch for {path.name}")
    raw = pd.read_csv(path, skiprows=[1], low_memory=False)
    required = {
        "time",
        "station_name",
        "mass_concentration_of_oxygen_in_sea_water",
    }
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"{path.name} missing columns {sorted(missing)}")
    raw_rows = len(raw)
    raw["time"] = pd.to_datetime(raw["time"], errors="coerce", utc=True)
    raw["do_mgL"] = pd.to_numeric(
        raw["mass_concentration_of_oxygen_in_sea_water"],
        errors="coerce",
    )
    if "sample_depth" in raw.columns:
        raw["sample_depth_numeric"] = pd.to_numeric(
            raw["sample_depth"],
            errors="coerce",
        )
    else:
        raw["sample_depth_numeric"] = np.nan

    # Identity correction frozen after a metadata-only audit exposed two
    # interleaved timestamp streams at Budds Landing. Other stations contain
    # long or periodic missing-depth segments rather than parallel records, so
    # they must not lose rows merely because the metadata flag changes.
    identity_rows = raw.loc[raw["time"].notna()].copy()
    identity_rows["sample_depth_present"] = identity_rows[
        "sample_depth_numeric"
    ].notna()
    present_rows = int(identity_rows["sample_depth_present"].sum())
    missing_rows = int((~identity_rows["sample_depth_present"]).sum())
    collision_count, minority_rows, minority_collision_share = (
        cross_group_collision_diagnostics(
            identity_rows.loc[identity_rows["sample_depth_present"], "time"],
            identity_rows.loc[~identity_rows["sample_depth_present"], "time"],
        )
    )
    parallel_stream_detected = (
        collision_count >= MIN_PARALLEL_COLLISIONS
        and minority_collision_share >= MIN_MINORITY_COLLISION_SHARE
    )
    select_present = present_rows >= missing_rows
    if parallel_stream_detected:
        selected_substream = (
            "sample_depth_present" if select_present else "sample_depth_missing"
        )
        selected_identity = identity_rows.loc[
            identity_rows["sample_depth_present"].eq(select_present)
        ].copy()
    else:
        selected_substream = "all_rows_no_parallel_collision"
        selected_identity = identity_rows.copy()
    parsed = selected_identity.loc[
        selected_identity["do_mgL"].notna()
        & selected_identity["do_mgL"].between(0.0, 50.0, inclusive="both")
    ].copy()
    parsed = parsed.sort_values("time", kind="stable")
    duplicate_rows = int(parsed.duplicated("time", keep=False).sum())
    parsed = parsed.drop_duplicates("time", keep="last").reset_index(drop=True)
    gaps = parsed["time"].diff().dt.total_seconds().div(3600.0)
    positive_gaps = gaps.loc[gaps.gt(0)]
    positive_gap_seconds = positive_gaps.mul(3600.0).round().astype(int)
    gap_counts = positive_gap_seconds.value_counts(dropna=True)
    modal_gap_seconds = int(gap_counts.index[0]) if not gap_counts.empty else None
    modal_gap_share = (
        float(gap_counts.iloc[0] / len(positive_gap_seconds))
        if len(positive_gap_seconds)
        else float("nan")
    )
    audit = {
        "dataset_id": item["dataset_id"],
        "candidate_year": item["candidate_year"],
        "split": "",
        "raw_rows": raw_rows,
        "time_parseable_sample_depth_present_rows": present_rows,
        "time_parseable_sample_depth_missing_rows": missing_rows,
        "cross_group_collision_count_lt_450s": collision_count,
        "minority_group_rows": minority_rows,
        "minority_group_collision_share": minority_collision_share,
        "parallel_stream_detected": parallel_stream_detected,
        "selected_substream": selected_substream,
        "selected_substream_rows_before_do_filter": len(selected_identity),
        "valid_unique_rows": len(parsed),
        "excluded_or_duplicate_rows": raw_rows - len(parsed),
        "duplicate_involved_rows": duplicate_rows,
        "modal_positive_gap_seconds": modal_gap_seconds,
        "modal_positive_gap_share": modal_gap_share,
        "first_time": str(parsed["time"].min()),
        "last_time": str(parsed["time"].max()),
        "do_min_mgL": float(parsed["do_mgL"].min()),
        "do_max_mgL": float(parsed["do_mgL"].max()),
        "median_positive_gap_hours": float(gaps.loc[gaps.gt(0)].median()),
        "share_positive_gaps_12_to_18min": float(
            positive_gaps.between(0.20, MAX_SLOPE_GAP_HOURS, inclusive="both").mean()
        ),
        "near_duplicate_positive_gaps_le_2min": int(positive_gaps.le(2 / 60).sum()),
        "share_positive_gaps_le_2min": float(positive_gaps.le(2 / 60).mean()),
        "share_gaps_le_18min": float(gaps.loc[gaps.gt(0)].le(MAX_SLOPE_GAP_HOURS).mean()),
        "csv_sha256": actual_sha,
    }
    return parsed, audit


def find_blocks(frame: pd.DataFrame) -> list[pd.DataFrame]:
    gaps = frame["time"].diff().dt.total_seconds().div(3600.0)
    run_id = (gaps.isna() | gaps.le(0) | gaps.gt(MAX_SLOPE_GAP_HOURS)).cumsum()
    blocks: list[pd.DataFrame] = []
    for _, run in frame.groupby(run_id, sort=True):
        run = run.reset_index(drop=True)
        for start in range(0, len(run) - WINDOW_POINTS + 1, WINDOW_POINTS - 1):
            candidate = run.iloc[start : start + WINDOW_POINTS].copy()
            elapsed = (
                candidate["time"] - candidate["time"].iloc[0]
            ).dt.total_seconds().to_numpy(dtype=float) / 3600.0
            adjacent = np.diff(elapsed)
            if (
                len(candidate) == WINDOW_POINTS
                and 719.5 <= elapsed[-1] <= 720.5
                and np.all(adjacent > 0)
                and np.all(adjacent <= MAX_SLOPE_GAP_HOURS)
            ):
                candidate["elapsed_hours"] = elapsed
                blocks.append(candidate)
                if len(blocks) >= MAX_BLOCKS_PER_STATION:
                    return blocks
    return blocks


def reference_burden(
    times: np.ndarray,
    values: np.ndarray,
    threshold: float,
) -> tuple[float, float, float]:
    normalized = times - times[0]
    segments: list[AffineSegment] = []
    slopes: list[float] = []
    for index, duration in enumerate(np.diff(normalized)):
        slope = float((values[index + 1] - values[index]) / duration)
        intercept = float(values[index] - slope * normalized[index])
        segments.append(
            AffineSegment(
                float(normalized[index]),
                float(normalized[index + 1]),
                slope,
                intercept,
            )
        )
        slopes.append(abs(slope))
    return (
        integrate_occupation(segments, threshold),
        integrate_oxygen_deficit(segments, threshold),
        max(slopes),
    )


def run() -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, object],
]:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest["status"] != "download_complete" or manifest["selected_count"] != 8:
        raise ValueError("independent-data manifest is not complete at eight stations")
    selected = sorted(manifest["selected"], key=lambda item: item["dataset_id"].casefold())

    station_frames: dict[str, pd.DataFrame] = {}
    audits: list[dict[str, object]] = []
    calibration_ids: list[str] = []
    validation_ids: list[str] = []
    for index, item in enumerate(selected):
        frame, audit = load_station(item)
        dataset_id = str(item["dataset_id"])
        split = "calibration" if index % 2 == 0 else "validation"
        audit["split"] = split
        station_frames[dataset_id] = frame
        audits.append(audit)
        (calibration_ids if split == "calibration" else validation_ids).append(dataset_id)

    calibration_slopes: list[np.ndarray] = []
    calibration_slope_rows: list[pd.DataFrame] = []
    for dataset_id in calibration_ids:
        frame = station_frames[dataset_id]
        gaps = frame["time"].diff().dt.total_seconds().div(3600.0)
        changes = frame["do_mgL"].diff().abs()
        eligible = gaps.gt(0) & gaps.le(MAX_SLOPE_GAP_HOURS)
        calibration_slopes.append((changes.loc[eligible] / gaps.loc[eligible]).to_numpy())
        slope_frame = pd.DataFrame(
            {
                "dataset_id": dataset_id,
                "previous_time": frame["time"].shift(1),
                "time": frame["time"],
                "gap_hours": gaps,
                "previous_do_mgL": frame["do_mgL"].shift(1),
                "do_mgL": frame["do_mgL"],
                "absolute_change_mgL": changes,
                "absolute_slope_mgL_per_hour": changes / gaps,
                "previous_sample_depth_m": frame["sample_depth_numeric"].shift(1),
                "sample_depth_m": frame["sample_depth_numeric"],
            }
        ).loc[eligible]
        slope_frame["sample_depth_changed"] = (
            slope_frame["previous_sample_depth_m"].notna()
            & slope_frame["sample_depth_m"].notna()
            & slope_frame["previous_sample_depth_m"].ne(slope_frame["sample_depth_m"])
        )
        calibration_slope_rows.append(slope_frame)
    slopes = np.concatenate(calibration_slopes)
    slopes = slopes[np.isfinite(slopes)]
    if slopes.size == 0:
        raise ValueError("no eligible calibration slopes")
    slope_parameters = {
        "calibration_q999": float(np.quantile(slopes, 0.999)),
        "calibration_max": float(np.max(slopes)),
    }
    slope_audit = pd.concat(calibration_slope_rows, ignore_index=True).sort_values(
        "absolute_slope_mgL_per_hour",
        ascending=False,
    )
    slope_station_summary = (
        slope_audit.groupby("dataset_id")["absolute_slope_mgL_per_hour"]
        .agg(
            n="size",
            q95=lambda values: values.quantile(0.95),
            q99=lambda values: values.quantile(0.99),
            q999=lambda values: values.quantile(0.999),
            maximum="max",
        )
        .reset_index()
    )

    validation_slope_rows: list[pd.DataFrame] = []
    for dataset_id in validation_ids:
        frame = station_frames[dataset_id]
        gaps = frame["time"].diff().dt.total_seconds().div(3600.0)
        changes = frame["do_mgL"].diff().abs()
        eligible = gaps.gt(0) & gaps.le(MAX_SLOPE_GAP_HOURS)
        validation_slope_rows.append(
            pd.DataFrame(
                {
                    "dataset_id": dataset_id,
                    "previous_time": frame["time"].shift(1),
                    "time": frame["time"],
                    "gap_hours": gaps,
                    "previous_do_mgL": frame["do_mgL"].shift(1),
                    "do_mgL": frame["do_mgL"],
                    "absolute_change_mgL": changes,
                    "absolute_slope_mgL_per_hour": changes / gaps,
                    "previous_sample_depth_m": frame["sample_depth_numeric"].shift(1),
                    "sample_depth_m": frame["sample_depth_numeric"],
                }
            ).loc[eligible]
        )
    validation_slope_audit = pd.concat(
        validation_slope_rows,
        ignore_index=True,
    ).sort_values("absolute_slope_mgL_per_hour", ascending=False)
    validation_slope_station_summary = (
        validation_slope_audit.groupby("dataset_id")["absolute_slope_mgL_per_hour"]
        .agg(
            n="size",
            q95=lambda values: values.quantile(0.95),
            q99=lambda values: values.quantile(0.99),
            q999=lambda values: values.quantile(0.999),
            maximum="max",
        )
        .reset_index()
    )

    blocks: list[dict[str, object]] = []
    for dataset_id in validation_ids:
        for block_index, block in enumerate(find_blocks(station_frames[dataset_id])):
            blocks.append(
                {
                    "dataset_id": dataset_id,
                    "block_index": block_index,
                    "start": str(block["time"].iloc[0]),
                    "end": str(block["time"].iloc[-1]),
                    "times": block["elapsed_hours"].to_numpy(dtype=float),
                    "values": block["do_mgL"].to_numpy(dtype=float),
                }
            )

    rows: list[dict[str, object]] = []
    for block in blocks:
        full_times = np.asarray(block["times"], dtype=float)
        full_values = np.asarray(block["values"], dtype=float)
        for cadence in CADENCE_HOURS:
            step = int(round(cadence / BASE_CADENCE_HOURS))
            for phase in range(step):
                sparse_indices = np.arange(phase, full_times.size, step, dtype=int)
                if sparse_indices.size < 2:
                    continue
                first_index = int(sparse_indices[0])
                last_index = int(sparse_indices[-1])
                reference_times = full_times[first_index : last_index + 1]
                reference_values = full_values[first_index : last_index + 1]
                sparse_times = full_times[sparse_indices]
                sparse_values = full_values[sparse_indices]
                for threshold in THRESHOLDS:
                    reference_occ, reference_deficit, reference_max_slope = reference_burden(
                        reference_times,
                        reference_values,
                        threshold,
                    )
                    for slope_kind, max_slope in slope_parameters.items():
                        base = {
                            "dataset_id": block["dataset_id"],
                            "block_index": block["block_index"],
                            "block_start": block["start"],
                            "block_end": block["end"],
                            "threshold_mgL": threshold,
                            "cadence_hours": cadence,
                            "phase_index": phase,
                            "slope_kind": slope_kind,
                            "max_slope_mgL_per_hour": max_slope,
                            "reference_max_slope_mgL_per_hour": reference_max_slope,
                            "reference_L_compatible": reference_max_slope <= max_slope + 1e-12,
                            "reference_occupation_hours": reference_occ,
                            "reference_deficit_mgL_hours": reference_deficit,
                        }
                        try:
                            result = identify_hypoxia_burden(
                                sparse_times,
                                sparse_values,
                                threshold=threshold,
                                max_slope=max_slope,
                            )
                        except InfeasibleObservationIntervals as exc:
                            rows.append({**base, "status": "infeasible", "error": str(exc)})
                            continue
                        occ_covered = (
                            result.occupation_lower - 1e-8
                            <= reference_occ
                            <= result.occupation_upper + 1e-8
                        )
                        deficit_covered = (
                            result.oxygen_deficit_lower - 1e-8
                            <= reference_deficit
                            <= result.oxygen_deficit_upper + 1e-8
                        )
                        rows.append(
                            {
                                **base,
                                "status": "feasible",
                                "error": "",
                                "horizon_hours": result.horizon,
                                "occupation_lower": result.occupation_lower,
                                "occupation_upper": result.occupation_upper,
                                "oxygen_deficit_lower": result.oxygen_deficit_lower,
                                "oxygen_deficit_upper": result.oxygen_deficit_upper,
                                "occupation_width_fraction": (
                                    result.occupation_upper - result.occupation_lower
                                )
                                / result.horizon,
                                "deficit_width_fraction_H_times_T": (
                                    result.oxygen_deficit_upper
                                    - result.oxygen_deficit_lower
                                )
                                / (threshold * result.horizon),
                                "occupation_reference_covered": occ_covered,
                                "deficit_reference_covered": deficit_covered,
                            }
                        )

    detail = pd.DataFrame(rows)
    audits_frame = pd.DataFrame(audits)
    if detail.empty:
        summary = pd.DataFrame()
        code_failures = 0
    else:
        group_keys = ["threshold_mgL", "cadence_hours", "slope_kind"]
        all_groups = detail.groupby(group_keys, dropna=False).agg(
            n_cases=("status", "size"),
            reference_L_compatible_share=("reference_L_compatible", "mean"),
        )
        feasible = detail.loc[detail["status"].eq("feasible")].copy()
        feasible_groups = feasible.groupby(group_keys, dropna=False).agg(
            n_feasible=("status", "size"),
            occupation_coverage=("occupation_reference_covered", "mean"),
            deficit_coverage=("deficit_reference_covered", "mean"),
            median_occupation_width_fraction=("occupation_width_fraction", "median"),
            mean_occupation_width_fraction=("occupation_width_fraction", "mean"),
            median_deficit_width_fraction=(
                "deficit_width_fraction_H_times_T",
                "median",
            ),
            mean_deficit_width_fraction=(
                "deficit_width_fraction_H_times_T",
                "mean",
            ),
        )
        compatible = feasible.loc[feasible["reference_L_compatible"]].copy()
        compatible_groups = compatible.groupby(group_keys, dropna=False).agg(
            n_L_compatible_feasible=("status", "size"),
            occupation_coverage_L_compatible=(
                "occupation_reference_covered",
                "mean",
            ),
            deficit_coverage_L_compatible=("deficit_reference_covered", "mean"),
            median_occupation_width_fraction_L_compatible=(
                "occupation_width_fraction",
                "median",
            ),
            median_deficit_width_fraction_L_compatible=(
                "deficit_width_fraction_H_times_T",
                "median",
            ),
        )
        summary = (
            all_groups.join(feasible_groups, how="left")
            .join(compatible_groups, how="left")
            .reset_index()
        )
        summary["feasible_share"] = summary["n_feasible"] / summary["n_cases"]
        # These columns are object-typed because infeasible rows carry missing
        # values.  Explicit False comparison avoids bitwise inversion turning
        # Python bool objects into -1/-2 integers.
        code_failures = int(
            compatible["occupation_reference_covered"].eq(False).sum()
            + compatible["deficit_reference_covered"].eq(False).sum()
        )

    metadata = {
        "status": (
            "blocked_no_complete_validation_blocks"
            if not blocks
            else "fail_code_containment"
            if code_failures
            else "pilot_complete_assumptions_not_validated"
        ),
        "source_manifest": str(MANIFEST),
        "source_manifest_sha256": sha256_file(MANIFEST),
        "calibration_dataset_ids": calibration_ids,
        "validation_dataset_ids": validation_ids,
        "calibration_eligible_slope_count": int(slopes.size),
        "slope_parameters_mgL_per_hour": slope_parameters,
        "slope_station_summary": slope_station_summary.to_dict(orient="records"),
        "validation_slope_station_summary": validation_slope_station_summary.to_dict(
            orient="records"
        ),
        "validation_eligible_slope_count": len(validation_slope_audit),
        "validation_share_slopes_le_calibration_q999": float(
            validation_slope_audit["absolute_slope_mgL_per_hour"]
            .le(slope_parameters["calibration_q999"])
            .mean()
        ),
        "validation_share_slopes_le_calibration_max": float(
            validation_slope_audit["absolute_slope_mgL_per_hour"]
            .le(slope_parameters["calibration_max"])
            .mean()
        ),
        "top_100_slope_depth_change_share": float(
            slope_audit.head(100)["sample_depth_changed"].mean()
        ),
        "validation_block_count": len(blocks),
        "validation_blocks": [
            {key: value for key, value in block.items() if key not in {"times", "values"}}
            for block in blocks
        ],
        "code_containment_failures_when_reference_L_compatible": code_failures,
        "reference_definition": "piecewise linear through complete 15-minute block",
        "reference_is_not_continuous_truth": True,
        "empirical_L_is_not_a_physical_guarantee": True,
        "stream_identity_rule": (
            "Split time-parseable rows by numeric sample_depth availability. "
            "Declare parallel streams only if at least 10 minority rows and at "
            "least 5% of the minority group have an opposite-group neighbor less "
            "than 450 seconds away. If parallel, retain the larger group (present "
            "wins ties); otherwise retain all rows. DO values are never inspected."
        ),
        "supersedes_initial_interleaved_stream_run": True,
    }
    return audits_frame, slope_audit, validation_slope_audit, detail, summary, metadata


def write_outputs() -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    audits, slope_audit, validation_slope_audit, detail, summary, metadata = run()
    audits.to_csv(OUTPUT / "station_data_audit.csv", index=False)
    slope_audit.head(1_000).to_csv(
        OUTPUT / "calibration_top_1000_slopes.csv",
        index=False,
    )
    validation_slope_audit.head(1_000).to_csv(
        OUTPUT / "validation_top_1000_slopes.csv",
        index=False,
    )
    detail.to_csv(OUTPUT / "masking_detail.csv", index=False)
    summary.to_csv(OUTPUT / "masking_summary.csv", index=False)
    (OUTPUT / "pilot_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    report = [
        "# MARACOOS/MDDNR collision-aware partial-identification pilot",
        "",
        f"Status: `{metadata['status']}`.",
        "",
        f"Calibration stations: {', '.join(metadata['calibration_dataset_ids'])}.",
        f"Validation stations: {', '.join(metadata['validation_dataset_ids'])}.",
        f"Complete 30-day validation blocks: {metadata['validation_block_count']}.",
        "Identity rule: split by numeric sample-depth availability; remove a subgroup "
        "only after repeated cross-group timestamp collisions below 450 seconds, "
        "without inspecting DO.",
        "",
        "Empirical sensitivity parameters (not physical guarantees):",
        "",
        f"- calibration q99.9 slope: {metadata['slope_parameters_mgL_per_hour']['calibration_q999']:.6g} mg/L/h",
        f"- calibration maximum slope: {metadata['slope_parameters_mgL_per_hour']['calibration_max']:.6g} mg/L/h",
        f"- sample-depth change share among top 100 slopes: {metadata['top_100_slope_depth_change_share']:.3f}",
        "",
        f"Code-containment failures when the operational reference was L-compatible: "
        f"{metadata['code_containment_failures_when_reference_L_compatible']}.",
        "",
    ]
    if not summary.empty:
        selected = summary.loc[summary["cadence_hours"].isin([0.5, 1.0, 2.0, 4.0])]
        report.extend(
            [
                "| H | Cadence h | L source | L-compatible share | Feasible share | Compatible n | Occupation width/T | Deficit width/(HT) |",
                "|---:|---:|---|---:|---:|---:|---:|---:|",
            ]
        )
        for row in selected.itertuples(index=False):
            report.append(
                f"| {row.threshold_mgL:g} | {row.cadence_hours:g} | {row.slope_kind} | "
                f"{row.reference_L_compatible_share:.3f} | {row.feasible_share:.3f} | "
                f"{row.n_L_compatible_feasible:.0f} | "
                f"{row.median_occupation_width_fraction_L_compatible:.3f} | "
                f"{row.median_deficit_width_fraction_L_compatible:.3f} |"
            )
    report.extend(
        [
            "",
            "Calibration slope distribution by station:",
            "",
            "| Dataset | n | q95 | q99 | q99.9 | max mg/L/h |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in metadata["slope_station_summary"]:
        report.append(
            f"| {row['dataset_id']} | {row['n']} | {row['q95']:.3f} | "
            f"{row['q99']:.3f} | {row['q999']:.3f} | {row['maximum']:.3f} |"
        )
    report.extend(
        [
            "",
            "Validation slope distribution by station (diagnostic only; not used to repair L):",
            "",
            "| Dataset | n | q95 | q99 | q99.9 | max mg/L/h |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in metadata["validation_slope_station_summary"]:
        report.append(
            f"| {row['dataset_id']} | {row['n']} | {row['q95']:.3f} | "
            f"{row['q99']:.3f} | {row['q999']:.3f} | {row['maximum']:.3f} |"
        )
    report.extend(
        [
            "",
            f"Validation slopes below calibration q99.9: "
            f"{metadata['validation_share_slopes_le_calibration_q999']:.5f}; "
            f"below calibration maximum: "
            f"{metadata['validation_share_slopes_le_calibration_max']:.5f}.",
        ]
    )
    report.extend(
        [
            "",
            "The 15-minute piecewise-linear record is an operational reference, not the "
            "unknown continuous truth. Validation data were not used to repair L.",
        ]
    )
    report_path = OUTPUT / "MARACOOS_PILOT_REPORT.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    return report_path


if __name__ == "__main__":
    print(write_outputs())
