#!/usr/bin/env python3
"""Verify package integrity, source snapshots, and regenerated outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CURRENT = ROOT / "03_代码与测试" / "current"
SNAPSHOT = ROOT / "03_代码与测试" / "原始快照"
MDDNR = ROOT / "04_真实数据" / "MARACOOS_MDDNR_discovery_v3"
USGS = ROOT / "04_真实数据" / "USGS_WDFN_external_2024"
EXTERNAL = ROOT / "05_既有结果" / "usgs_external_7day_replication"
REFERENCE = ROOT / "reference_outputs"
VOLATILE_COLUMNS = {"runtime_seconds", "median_runtime_seconds"}

if str(SNAPSHOT) not in sys.path:
    sys.path.insert(0, str(SNAPSHOT))

import future_partial_identification_maracoos_pilot as _pilot  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_integrity() -> None:
    entries = json.loads((ROOT / "SHA256SUMS.json").read_text(encoding="utf-8"))
    failures = []
    for relative, expected in entries.items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
        elif sha256(path) != expected:
            failures.append(f"hash mismatch: {relative}")
    if failures:
        raise SystemExit("package integrity failed:\n" + "\n".join(failures))


def verify_mddnr_data() -> None:
    manifest_path = MDDNR / "download_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("catalog_license") != "CC0-1.0":
        raise SystemExit("unexpected MDDNR catalog licence")
    if manifest.get("selected_count") != 8:
        raise SystemExit("expected eight MDDNR selected stations")
    for item in manifest["selected"]:
        csv_path = MDDNR / item["csv_file"]
        info_path = MDDNR / item["info_file"]
        if sha256(csv_path) != item["csv_sha256"]:
            raise SystemExit(f"MDDNR data hash mismatch: {csv_path.name}")
        if sha256(info_path) != item["info_sha256"]:
            raise SystemExit(f"MDDNR metadata hash mismatch: {info_path.name}")


def verify_usgs_data() -> None:
    manifest = json.loads((USGS / "download_manifest.json").read_text(encoding="utf-8"))
    expected = (
        "USGS-01491000", "USGS-01579550", "USGS-01594441",
        "USGS-01638500", "USGS-01643580", "USGS-01646500",
        "USGS-01649190", "USGS-01649500", "USGS-01650800",
    )
    observed = tuple(str(item["monitoring_location_id"]) for item in manifest["series"])
    if observed != expected:
        raise SystemExit("USGS frozen cohort identity/order mismatch")
    if manifest.get("source_authority") != "U.S. Geological Survey":
        raise SystemExit("unexpected USGS authority")
    if manifest.get("parameter_code") != "00300" or manifest.get("unit") != "mg/l":
        raise SystemExit("unexpected USGS parameter identity")
    if "10.5066/F7P55KJN" not in manifest.get("recommended_dataset_citation", ""):
        raise SystemExit("USGS dataset citation DOI missing")
    for item in manifest["series"]:
        raw = USGS / item["raw_file"]
        normalized = USGS / item["normalized_file"]
        if sha256(raw) != item["raw_sha256"]:
            raise SystemExit(f"USGS raw hash mismatch: {raw.name}")
        if sha256(normalized) != item["normalized_sha256"]:
            raise SystemExit(f"USGS normalized hash mismatch: {normalized.name}")
        approved = item.get("qualifier_definitions", {}).get("A", "")
        if "Approved for publication" not in approved:
            raise SystemExit("USGS approved qualifier definition missing")


def verify_calibration_slope() -> None:
    manifest = json.loads((MDDNR / "download_manifest.json").read_text(encoding="utf-8"))
    metadata_path = (
        ROOT / "05_既有结果" /
        "partial_identification_maracoos_pilot_v3_collision_identity" /
        "pilot_metadata.json"
    )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    selected = {str(item["dataset_id"]): item for item in manifest["selected"]}
    _pilot.DATA = MDDNR
    slope_arrays = []
    for dataset_id in metadata["calibration_dataset_ids"]:
        frame, _audit = _pilot.load_station(selected[dataset_id])
        gaps = frame["time"].diff().dt.total_seconds().div(3600.0)
        changes = frame["do_mgL"].diff().abs()
        eligible = gaps.gt(0) & gaps.le(_pilot.MAX_SLOPE_GAP_HOURS)
        slope_arrays.append((changes.loc[eligible] / gaps.loc[eligible]).to_numpy())
    slopes = np.concatenate(slope_arrays)
    slopes = slopes[np.isfinite(slopes)]
    q999 = float(np.quantile(slopes, 0.999))
    maximum = float(np.max(slopes))
    expected = metadata["slope_parameters_mgL_per_hour"]
    if slopes.size != 28968:
        raise SystemExit(f"unexpected calibration slope count: {slopes.size}")
    if not math.isclose(q999, float(expected["calibration_q999"]), rel_tol=0.0, abs_tol=1e-14):
        raise SystemExit("calibration q99.9 slope mismatch")
    if not math.isclose(maximum, float(expected["calibration_max"]), rel_tol=0.0, abs_tol=1e-14):
        raise SystemExit("calibration maximum slope mismatch")


def parse_number(value: str) -> float | None:
    if value == "":
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def compare_csv(actual: Path, expected: Path) -> None:
    with actual.open("r", encoding="utf-8", newline="") as handle:
        actual_rows = list(csv.DictReader(handle))
    with expected.open("r", encoding="utf-8", newline="") as handle:
        expected_rows = list(csv.DictReader(handle))
    if len(actual_rows) != len(expected_rows):
        raise SystemExit(f"row-count mismatch for {actual.name}: {len(actual_rows)} != {len(expected_rows)}")
    for index, (left, right) in enumerate(zip(actual_rows, expected_rows)):
        columns = (set(left) | set(right)) - VOLATILE_COLUMNS
        for column in sorted(columns):
            a = left.get(column, "")
            b = right.get(column, "")
            if a == b:
                continue
            an = parse_number(a)
            bn = parse_number(b)
            if an is not None and bn is not None and math.isclose(an, bn, rel_tol=1e-11, abs_tol=1e-11):
                continue
            raise SystemExit(
                f"value mismatch {actual.name} row {index + 2} "
                f"column {column}: {a!r} != {b!r}"
            )


def verify_external_receipt() -> None:
    receipt = json.loads(
        (EXTERNAL / "USGS_EXTERNAL_REPLICATION_VERIFICATION.json").read_text(encoding="utf-8")
    )
    expected_counts = {
        "selected_series_count": 9,
        "series_with_qualifying_blocks": 9,
        "qualifying_block_count": 36,
        "compatible_block_count": 36,
        "all_phase_case_count": 1080,
        "compatible_k0_case_count": 1080,
        "fixed_k1_case_count": 1080,
        "containment_failures": 0,
        "nesting_failures": 0,
    }
    if receipt.get("status") != "passed":
        raise SystemExit("USGS verification receipt is not passed")
    for key, expected in expected_counts.items():
        if receipt.get(key) != expected:
            raise SystemExit(f"USGS receipt mismatch: {key}")
    for name, expected in receipt["output_hashes"].items():
        if sha256(EXTERNAL / name) != expected:
            raise SystemExit(f"USGS frozen output hash mismatch: {name}")


def verify_frozen_release_outputs() -> None:
    csv_paths = (
        "05_既有结果/usgs_external_7day_replication/usgs_external_k0_detail.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_k1_detail.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_cadence_summary.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_block_cadence_summary.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_station_block_summary.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_data_audit.csv",
        "09_论文写作输出/tables/Table_S1_MDDNR_data_audit_revision_o.csv",
        "09_论文写作输出/tables/Table_2_revision_o_block_level.csv",
        "09_论文写作输出/tables/Table_3_USGS_external_replication.csv",
    )
    for relative in csv_paths:
        compare_csv(ROOT / relative, REFERENCE / relative)
    text_paths = (
        "05_既有结果/usgs_external_7day_replication/usgs_external_metadata.json",
        "05_既有结果/usgs_external_7day_replication/USGS_EXTERNAL_REPLICATION_REPORT.md",
        "09_论文写作输出/tables/Table_2_revision_o_block_level.md",
        "09_论文写作输出/tables/Table_3_USGS_external_replication.md",
        "09_论文写作输出/SUPPLEMENT_TABLE_S1_DATA_PROTOCOL_REVISION_O.md",
        "09_论文写作输出/SUPPLEMENT_TABLE_S4_USGS_DATA_AUDIT.md",
        "09_论文写作输出/SUPPLEMENT_TABLE_S5_USGS_STATION_CADENCE.md",
    )
    for relative in text_paths:
        left = (ROOT / relative).read_text(encoding="utf-8")
        right = (REFERENCE / relative).read_text(encoding="utf-8")
        if left != right:
            raise SystemExit(f"frozen release output mismatch: {relative}")


def verify_regenerated_outputs() -> None:
    csv_paths = (
        "05_既有结果/current_maracoos_7day_nonnegative_table2/table2_nonnegative_detail.csv",
        "05_既有结果/current_maracoos_7day_nonnegative_table2/table2_nonnegative_summary.csv",
        "05_既有结果/current_maracoos_24h_nonnegative_trilemma/trilemma_nonnegative_detail.csv",
        "05_既有结果/current_maracoos_7day_fixed_k1_sensitivity/fixed_k1_detail.csv",
        "05_既有结果/current_maracoos_7day_fixed_k1_sensitivity/fixed_k1_summary.csv",
        "09_论文写作输出/tables/Table_1_model_assumptions.csv",
        "09_论文写作输出/tables/Table_2_7day_nonnegative.csv",
        "09_论文写作输出/tables/Table_S2_fixed_k1_sensitivity.csv",
        "09_论文写作输出/tables/Table_S3_block_level_widths.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_k0_detail.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_k1_detail.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_cadence_summary.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_block_cadence_summary.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_station_block_summary.csv",
        "05_既有结果/usgs_external_7day_replication/usgs_external_data_audit.csv",
        "09_论文写作输出/tables/Table_S1_MDDNR_data_audit_revision_o.csv",
        "09_论文写作输出/tables/Table_2_revision_o_block_level.csv",
        "09_论文写作输出/tables/Table_3_USGS_external_replication.csv",
    )
    for relative in csv_paths:
        compare_csv(ROOT / relative, REFERENCE / relative)

    binary_paths = (
        "09_论文写作输出/figures/Figure_1_endpoint_identification_concept.svg",
        "09_论文写作输出/figures/Figure_2_L_k_width_nonnegative.svg",
        "09_论文写作输出/figures/Figure_3_common_anchor_geometry.svg",
    )
    for relative in binary_paths:
        if (ROOT / relative).read_bytes() != (REFERENCE / relative).read_bytes():
            raise SystemExit(f"canonical binary output mismatch: {relative}")

    text_paths = (
        "05_既有结果/usgs_external_7day_replication/usgs_external_metadata.json",
        "05_既有结果/usgs_external_7day_replication/USGS_EXTERNAL_REPLICATION_REPORT.md",
        "09_论文写作输出/tables/Table_2_revision_o_block_level.md",
        "09_论文写作输出/tables/Table_3_USGS_external_replication.md",
        "09_论文写作输出/SUPPLEMENT_TABLE_S1_DATA_PROTOCOL_REVISION_O.md",
        "09_论文写作输出/SUPPLEMENT_TABLE_S4_USGS_DATA_AUDIT.md",
        "09_论文写作输出/SUPPLEMENT_TABLE_S5_USGS_STATION_CADENCE.md",
    )
    for relative in text_paths:
        left = (ROOT / relative).read_text(encoding="utf-8")
        right = (REFERENCE / relative).read_text(encoding="utf-8")
        if left != right:
            raise SystemExit(f"canonical text output mismatch: {relative}")

    detail = ROOT / csv_paths[0]
    with detail.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    feasible = [row for row in rows if row["status"] == "feasible"]
    compatible = [row for row in feasible if row["reference_L_compatible"] == "True"]
    if (len(rows), len(feasible), len(compatible)) != (480, 468, 240):
        raise SystemExit("unexpected MDDNR Table 2 cohort counts")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-only", action="store_true")
    args = parser.parse_args()
    verify_integrity()
    verify_mddnr_data()
    verify_usgs_data()
    verify_calibration_slope()
    verify_external_receipt()
    verify_frozen_release_outputs()
    expected_runtime = "5e06c11ea2bff3e1fd2713017736cb8343fabd4d9b4b332af234025a45891b2a"
    if sha256(CURRENT / "robust_exact_current.py") != expected_runtime:
        raise SystemExit("scientific runtime hash mismatch")
    if not args.package_only:
        verify_regenerated_outputs()
    print("verification: PASS")


if __name__ == "__main__":
    main()
