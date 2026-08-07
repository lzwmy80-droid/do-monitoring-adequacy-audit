# -*- coding: utf-8 -*-
"""Build block-level and provenance assets for the EMAS revision O draft."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import run_maracoos_nonnegative_table2 as baseline


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
WRITING = PROJECT / "09_论文写作输出"
TABLES = WRITING / "tables"
MDDNR_OUTPUT = (
    PROJECT
    / "05_既有结果"
    / "partial_identification_maracoos_pilot_v3_collision_identity"
)
MDDNR_METADATA = MDDNR_OUTPUT / "pilot_metadata.json"
MDDNR_AUDIT_CSV = TABLES / "Table_S1_MDDNR_data_audit_revision_o.csv"
MDDNR_S1_MD = WRITING / "SUPPLEMENT_TABLE_S1_DATA_PROTOCOL_REVISION_O.md"
TABLE2_CSV = TABLES / "Table_2_revision_o_block_level.csv"
TABLE2_MD = TABLES / "Table_2_revision_o_block_level.md"
EXTERNAL = PROJECT / "05_既有结果" / "usgs_external_7day_replication"
TABLE3_CSV = TABLES / "Table_3_USGS_external_replication.csv"
TABLE3_MD = TABLES / "Table_3_USGS_external_replication.md"
USGS_S4_MD = WRITING / "SUPPLEMENT_TABLE_S4_USGS_DATA_AUDIT.md"
USGS_S5_MD = WRITING / "SUPPLEMENT_TABLE_S5_USGS_STATION_CADENCE.md"


def _fmt(value: float) -> str:
    return f"{float(value):.3f}"


def build_mddnr_audit() -> pd.DataFrame:
    manifest = json.loads(baseline.MANIFEST.read_text(encoding="utf-8"))
    metadata = json.loads(MDDNR_METADATA.read_text(encoding="utf-8"))
    calibration = set(metadata["calibration_dataset_ids"])
    evaluation = set(metadata["validation_dataset_ids"])
    selected = sorted(
        manifest["selected"], key=lambda item: str(item["dataset_id"]).casefold()
    )
    alternating_roles = [
        "slope-label partition" if index % 2 == 0 else "evaluation partition"
        for index in range(len(selected))
    ]
    expected_roles = [
        "slope-label partition"
        if str(item["dataset_id"]) in calibration
        else "evaluation partition"
        for item in selected
    ]
    if alternating_roles != expected_roles:
        raise ValueError("frozen split is not the documented alternating split")
    if calibration | evaluation != {str(item["dataset_id"]) for item in selected}:
        raise ValueError("frozen split does not cover every selected dataset")

    baseline._pilot.DATA = baseline.DATA
    baseline._pilot.MANIFEST = baseline.MANIFEST
    rows: list[dict[str, object]] = []
    for index, (item, role) in enumerate(zip(selected, expected_roles)):
        _frame, audit = baseline._pilot.load_station(item)
        dataset_id = str(item["dataset_id"])
        rows.append(
            {
                "sorted_index_zero_based": index,
                "role": role,
                "dataset_id": dataset_id,
                "year": int(item["candidate_year"]),
                "raw_download_rows": int(item["raw_data_rows"]),
                "manifest_numeric_DO_rows_before_identity_resolution": int(
                    item["valid_do_rows"]
                ),
                "identity_selected_rows_before_DO_filter": int(
                    audit["selected_substream_rows_before_do_filter"]
                ),
                "final_valid_unique_analysis_rows": int(audit["valid_unique_rows"]),
                "excluded_or_duplicate_rows_from_raw": int(
                    audit["excluded_or_duplicate_rows"]
                ),
                "parallel_stream_detected": bool(audit["parallel_stream_detected"]),
                "selected_substream": str(audit["selected_substream"]),
                "cross_group_collision_count_lt_450s": int(
                    audit["cross_group_collision_count_lt_450s"]
                ),
                "minority_group_collision_share": float(
                    audit["minority_group_collision_share"]
                ),
                "duplicate_involved_rows": int(audit["duplicate_involved_rows"]),
                "first_analysis_time_utc": str(audit["first_time"]),
                "last_analysis_time_utc": str(audit["last_time"]),
                "latitude": str(item["latitudes"][0]),
                "longitude": str(item["longitudes"][0]),
                "csv_sha256": str(item["csv_sha256"]),
                "erddap_dataset_url": (
                    "https://erddap.maracoos.org/erddap/tabledap/"
                    f"{dataset_id}.html"
                ),
                "exact_download_url": str(item["data_url"]),
            }
        )
    frame = pd.DataFrame(rows)
    TABLES.mkdir(parents=True, exist_ok=True)
    frame.to_csv(MDDNR_AUDIT_CSV, index=False)

    lines = [
        "# Table S1. Frozen MDDNR data inventory and complete preprocessing audit",
        "",
        "Data were downloaded on 16 July 2026 from the MARACOOS ERDDAP delivery "
        "of the Maryland Department of Natural Resources (MDDNR) Continuous "
        "Monitoring Program. The Data.gov catalog describes these station data "
        "as public and CC0-1.0. The frozen download-manifest SHA-256 is "
        f"`{baseline._sha256(baseline.MANIFEST)}`.",
        "",
        "Eligible dataset identifiers were sorted lexicographically. Even "
        "zero-based positions (0, 2, 4, 6) supplied slope-scenario labels; odd "
        "positions (1, 3, 5, 7) formed a deterministic evaluation partition. "
        "DO values did not enter station selection or partition assignment. "
        "Because Harris Creek Downstream and Upstream occur in different "
        "partitions, the evaluation partition is not claimed to be spatially "
        "independent or an external validation sample.",
        "",
        "| Role | Exact dataset ID | Year | Raw rows | Identity-selected rows | "
        "Final valid unique rows | Excluded from raw | Latitude | Longitude |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in frame.itertuples(index=False):
        lines.append(
            f"| {row.role} | `{row.dataset_id}` | {row.year} | "
            f"{row.raw_download_rows:,} | "
            f"{row.identity_selected_rows_before_DO_filter:,} | "
            f"{row.final_valid_unique_analysis_rows:,} | "
            f"{row.excluded_or_duplicate_rows_from_raw:,} | "
            f"{row.latitude} | {row.longitude} |"
        )
    lines.extend(
        [
            "",
            "## Exact preprocessing rules",
            "",
            "- Parse timestamps as UTC and DO as numeric. Rows with an "
            "unparseable timestamp, missing/non-numeric DO, or DO outside the "
            "inclusive plausibility screen `[0, 50] mg/L` do not enter analysis.",
            "- Diagnose a possible interleaved stream by splitting timestamped "
            "rows according to numeric sample-depth availability. A parallel "
            "stream is declared only when at least 10 minority-group rows have "
            "an opposite-group timestamp less than 450 s away and these "
            "collisions involve at least 5% of the minority group. If declared, "
            "retain the larger group (sample-depth-present wins a tie); otherwise "
            "retain all rows. This rule does not inspect DO values.",
            "- Sort stably by UTC time, diagnose duplicate-involved timestamps, "
            "and retain the last row at an exact duplicate timestamp. No input "
            "value is interpolated.",
            "- The manifest's numeric-DO count precedes identity resolution. "
            "For Budds Landing, 10,275 downloaded numeric-DO rows contain two "
            "interleaved streams; the frozen identity rule retains 8,833 rows "
            "for analysis and excludes 1,442 rows. The table therefore separates "
            "raw, identity-selected, and final analysis counts.",
            "- MDDNR annual metadata documents calibration/in-situ comparison, "
            "flag review, masking, and season-end QA/QC. The ±0.5 mg/L DO value "
            "is a comparison tolerance used in that QA process, not a pointwise "
            "error bound assumed by this analysis.",
            "",
            "## Per-dataset source and checksum register",
            "",
        ]
    )
    for row in frame.itertuples(index=False):
        lines.extend(
            [
                f"- `{row.dataset_id}` — CSV SHA-256 "
                f"`{row.csv_sha256}`; dataset page "
                f"<{row.erddap_dataset_url}>; exact frozen query "
                f"<{row.exact_download_url}>.",
            ]
        )
    lines.extend(
        [
            "",
            "Annual MDDNR process metadata: "
            "<https://eyesonthebay.dnr.maryland.gov/eyesonthebay/documents/"
            "metadata/MdDNR2017CMonProj.html>, "
            "<https://eyesonthebay.dnr.maryland.gov/eyesonthebay/documents/"
            "metadata/MdDNR2018CMonProj.html>, and "
            "<https://eyesonthebay.dnr.maryland.gov/eyesonthebay/documents/"
            "metadata/MdDNR2019CMonProj.html>.",
        ]
    )
    MDDNR_S1_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return frame


def build_table2() -> pd.DataFrame:
    detail = pd.read_csv(
        PROJECT
        / "05_既有结果"
        / "current_maracoos_7day_nonnegative_table2"
        / "table2_nonnegative_detail.csv"
    )
    detail = detail.loc[
        detail["status"].eq("feasible")
        & detail["reference_L_compatible"].eq(True)
    ].copy()
    block = (
        detail.groupby(
            ["dataset_id", "block_index", "cadence_hours"], as_index=False
        )
        .agg(
            compatible_phase_cases_n=("phase_index", "size"),
            median_occupation_width_over_T=(
                "occupation_width_fraction",
                "median",
            ),
            median_nonnegative_deficit_width_over_HT=(
                "deficit_width_fraction_nonnegative_H_times_T",
                "median",
            ),
        )
    )
    k1 = pd.read_csv(
        PROJECT
        / "05_既有结果"
        / "current_maracoos_7day_fixed_k1_sensitivity"
        / "fixed_k1_detail.csv"
    )
    k1_block = (
        k1.groupby(["dataset_id", "block_index", "cadence_hours"], as_index=False)
        .agg(
            median_occupation_k1_increase=(
                "occupation_width_fraction_increase",
                "median",
            ),
            median_deficit_k1_increase=(
                "deficit_width_fraction_increase",
                "median",
            ),
        )
    )
    block = block.merge(
        k1_block,
        on=["dataset_id", "block_index", "cadence_hours"],
        how="left",
        validate="one_to_one",
    )
    summary = (
        block.groupby("cadence_hours", as_index=False)
        .agg(
            compatible_station_blocks=("block_index", "size"),
            compatible_stations=("dataset_id", "nunique"),
            within_block_phase_cases=("compatible_phase_cases_n", "sum"),
            median_block_occupation_width=(
                "median_occupation_width_over_T",
                "median",
            ),
            q25_block_occupation_width=(
                "median_occupation_width_over_T",
                lambda values: float(values.quantile(0.25)),
            ),
            q75_block_occupation_width=(
                "median_occupation_width_over_T",
                lambda values: float(values.quantile(0.75)),
            ),
            median_block_deficit_width=(
                "median_nonnegative_deficit_width_over_HT",
                "median",
            ),
            q25_block_deficit_width=(
                "median_nonnegative_deficit_width_over_HT",
                lambda values: float(values.quantile(0.25)),
            ),
            q75_block_deficit_width=(
                "median_nonnegative_deficit_width_over_HT",
                lambda values: float(values.quantile(0.75)),
            ),
            median_block_occupation_k1_increase=(
                "median_occupation_k1_increase",
                "median",
            ),
            median_block_deficit_k1_increase=(
                "median_deficit_k1_increase",
                "median",
            ),
        )
    )
    summary.to_csv(TABLE2_CSV, index=False)
    lines = [
        "**Table 2. Station-block-level seven-day uncertainty widths in the "
        "deterministic MDDNR evaluation partition.** Main values are medians "
        "across eight compatible station-block phase medians; brackets give "
        "the interquartile range across station-blocks. Overlapping phases are "
        "within-block numerical scenarios, not independent replicates.",
        "",
        "| Cadence (h) | Compatible station-blocks (stations) | Phase scenarios | "
        "Occupation width/T, median [IQR] | Deficit width/(HT), median [IQR] | "
        "Median paired k=1 increase, occupation; deficit |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.cadence_hours:g} | "
            f"{int(row.compatible_station_blocks)} ({int(row.compatible_stations)}) | "
            f"{int(row.within_block_phase_cases)} | "
            f"{_fmt(row.median_block_occupation_width)} "
            f"[{_fmt(row.q25_block_occupation_width)}, "
            f"{_fmt(row.q75_block_occupation_width)}] | "
            f"{_fmt(row.median_block_deficit_width)} "
            f"[{_fmt(row.q25_block_deficit_width)}, "
            f"{_fmt(row.q75_block_deficit_width)}] | "
            f"{_fmt(row.median_block_occupation_k1_increase)}; "
            f"{_fmt(row.median_block_deficit_k1_increase)} |"
        )
    TABLE2_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary


def build_external_assets() -> pd.DataFrame | None:
    cadence_path = EXTERNAL / "usgs_external_cadence_summary.csv"
    if not cadence_path.exists():
        return None
    cadence = pd.read_csv(cadence_path)
    audit = pd.read_csv(EXTERNAL / "usgs_external_data_audit.csv")
    block_cadence = pd.read_csv(
        EXTERNAL / "usgs_external_block_cadence_summary.csv"
    )
    metadata = json.loads(
        (EXTERNAL / "usgs_external_metadata.json").read_text(encoding="utf-8")
    )

    main = cadence[
        [
            "cadence_hours",
            "compatible_station_blocks",
            "compatible_stations",
            "within_block_phase_cases",
            "median_block_occupation_width_fraction",
            "q25_block_occupation_width_fraction",
            "q75_block_occupation_width_fraction",
            "median_block_deficit_width_fraction",
            "q25_block_deficit_width_fraction",
            "q75_block_deficit_width_fraction",
            "median_block_occupation_k1_increase",
            "median_block_deficit_k1_increase",
        ]
    ].copy()
    main.to_csv(TABLE3_CSV, index=False)
    lines = [
        "**Table 3. Frozen-label USGS cross-network transfer stress test.** "
        "All nine metadata-eligible Maryland USGS dissolved-oxygen series were "
        "retained before endpoint computation. Values are medians across "
        "compatible station-block phase medians, with station-block IQRs. This "
        "is a deterministic transportability check, not external validation or "
        "population inference.",
        "",
        "| Cadence (h) | Compatible station-blocks (stations) | Phase scenarios | "
        "Occupation width/T, median [IQR] | Deficit width/(HT), median [IQR] | "
        "Median paired k=1 increase, occupation; deficit |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in main.itertuples(index=False):
        lines.append(
            f"| {row.cadence_hours:g} | "
            f"{int(row.compatible_station_blocks)} ({int(row.compatible_stations)}) | "
            f"{int(row.within_block_phase_cases)} | "
            f"{_fmt(row.median_block_occupation_width_fraction)} "
            f"[{_fmt(row.q25_block_occupation_width_fraction)}, "
            f"{_fmt(row.q75_block_occupation_width_fraction)}] | "
            f"{_fmt(row.median_block_deficit_width_fraction)} "
            f"[{_fmt(row.q25_block_deficit_width_fraction)}, "
            f"{_fmt(row.q75_block_deficit_width_fraction)}] | "
            f"{_fmt(row.median_block_occupation_k1_increase)}; "
            f"{_fmt(row.median_block_deficit_k1_increase)} |"
        )
    TABLE3_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    s4 = [
        "# Table S4. USGS external-network source and preprocessing audit",
        "",
        "The cohort was frozen from metadata before endpoint computation. "
        "All raw responses and adapted files are hash-locked. Only the legacy "
        "NWIS qualifier exactly `A` (approved for publication after processing "
        "and review) entered the frozen UTC interval. Native 5-min records were "
        "subselected at exact UTC quarter hours without averaging or "
        "interpolation. The NWIS instantaneous series statistic code is `00000`.",
        "",
        "| USGS site | Native gap (min) | Raw response rows | Frozen-interval rows | "
        "Approved-A rows | Quarter-hour rows | 7-day blocks | DO range (mg/L) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audit.itertuples(index=False):
        s4.append(
            f"| `{row.monitoring_location_id}` | "
            f"{row.native_modal_gap_minutes:g} | "
            f"{int(row.raw_response_rows):,} | "
            f"{int(row.rows_in_frozen_interval):,} | "
            f"{int(row.exact_A_rows_in_interval):,} | "
            f"{int(row.quarter_hour_reference_rows):,} | "
            f"{int(row.qualifying_7day_blocks)} | "
            f"{row.minimum_DO_mgL_in_reference:.2f}–"
            f"{row.maximum_DO_mgL_in_reference:.2f} |"
        )
    s4.extend(["", "## Exact URLs and checksums", ""])
    for row in audit.itertuples(index=False):
        s4.append(
            f"- `{row.monitoring_location_id}` — {row.site_name}; time-series "
            f"ID `{row.time_series_id}`; statistic `00000`: raw "
            f"SHA-256 `{row.raw_sha256}`; normalized SHA-256 "
            f"`{row.normalized_sha256}`; query <{row.source_query_url}>."
        )
    s4.extend(
        [
            "",
            "Recommended source citation: U.S. Geological Survey (2026), "
            "U.S. Geological Survey National Water Information System "
            "database, accessed 2 August 2026, "
            "<https://doi.org/10.5066/F7P55KJN>.",
        ]
    )
    USGS_S4_MD.write_text("\n".join(s4) + "\n", encoding="utf-8")

    station_cadence = (
        block_cadence.groupby(
            ["monitoring_location_id", "site_name", "cadence_hours"],
            as_index=False,
        )
        .agg(
            compatible_blocks=("block_index", "nunique"),
            phase_cases=("phase_cases", "sum"),
            median_occupation_width=(
                "median_occupation_width_fraction",
                "median",
            ),
            median_deficit_width=(
                "median_deficit_width_fraction_nonnegative",
                "median",
            ),
            median_occupation_k1_increase=(
                "median_occupation_width_fraction_increase",
                "median",
            ),
            median_deficit_k1_increase=(
                "median_deficit_width_fraction_increase",
                "median",
            ),
        )
    )
    s5 = [
        "# Table S5. USGS station-by-cadence descriptive results",
        "",
        "Each row summarizes compatible seven-day blocks within one selected "
        "USGS series and cadence. Values first take the median over phases "
        "within a block and then the median across that station's blocks. "
        "Neither phases nor blocks are treated as population replicates. Full "
        "station-block results are supplied in the machine-readable output.",
        "",
        "| USGS site | Cadence (h) | Blocks | Phase cases | Occupation width/T | "
        "Deficit width/(HT) | k=1 increase: occupation; deficit |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in station_cadence.itertuples(index=False):
        s5.append(
            f"| `{row.monitoring_location_id}` | {row.cadence_hours:g} | "
            f"{int(row.compatible_blocks)} | {int(row.phase_cases)} | "
            f"{_fmt(row.median_occupation_width)} | "
            f"{_fmt(row.median_deficit_width)} | "
            f"{_fmt(row.median_occupation_k1_increase)}; "
            f"{_fmt(row.median_deficit_k1_increase)} |"
        )
    s5.extend(
        [
            "",
            f"Cohort accounting: {metadata['series_with_qualifying_blocks']}/"
            f"{metadata['selected_series_count']} selected series supplied a "
            f"qualifying block; {metadata['compatible_block_count']}/"
            f"{metadata['qualifying_block_count']} blocks were compatible with "
            "the unchanged Maryland-derived L label.",
        ]
    )
    USGS_S5_MD.write_text("\n".join(s5) + "\n", encoding="utf-8")
    return main


def main() -> None:
    mddnr = build_mddnr_audit()
    table2 = build_table2()
    table3 = build_external_assets()
    print(
        json.dumps(
            {
                "mddnr_rows": len(mddnr),
                "table2_rows": len(table2),
                "external_ready": table3 is not None,
                "table3_rows": 0 if table3 is None else len(table3),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
