# -*- coding: utf-8 -*-
"""Verify the frozen USGS external-network replication artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
DATA = PROJECT / "04_真实数据" / "USGS_WDFN_external_2024"
OUTPUT = PROJECT / "05_既有结果" / "usgs_external_7day_replication"
FREEZE_NOTE = (
    PROJECT
    / "09_论文写作输出"
    / "USGS_EXTERNAL_REPLICATION_FREEZE_NOTE_20260802.md"
)
CORE = HERE / "robust_exact_current.py"
ANALYSIS = HERE / "run_usgs_external_replication.py"
EXPECTED_SITES = (
    "USGS-01491000",
    "USGS-01579550",
    "USGS-01594441",
    "USGS-01638500",
    "USGS-01643580",
    "USGS-01646500",
    "USGS-01649190",
    "USGS-01649500",
    "USGS-01650800",
)
EXPECTED_PHASES = {0.5: 2, 1.0: 4, 2.0: 8, 4.0: 16}
TOLERANCE = 1e-8


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _case_keys(frame: pd.DataFrame) -> set[tuple[str, int, float, int]]:
    return {
        (
            str(row.monitoring_location_id),
            int(row.block_index),
            float(row.cadence_hours),
            int(row.phase_index),
        )
        for row in frame.itertuples(index=False)
    }


def verify() -> Path:
    manifest_path = DATA / "download_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    sites = tuple(
        str(item["monitoring_location_id"]) for item in manifest["series"]
    )
    if sites != EXPECTED_SITES:
        raise AssertionError("selected USGS site denominator changed")
    for item in manifest["series"]:
        raw = DATA / str(item["raw_file"])
        normalized = DATA / str(item["normalized_file"])
        if _sha256(raw) != str(item["raw_sha256"]):
            raise AssertionError(f"raw hash failed for {item['monitoring_location_id']}")
        if _sha256(normalized) != str(item["normalized_sha256"]):
            raise AssertionError(
                f"normalized hash failed for {item['monitoring_location_id']}"
            )

    metadata_path = OUTPUT / "usgs_external_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    hash_checks = {
        "freeze_note_sha256": _sha256(FREEZE_NOTE),
        "download_manifest_sha256": _sha256(manifest_path),
        "core_runtime_sha256": _sha256(CORE),
        "analysis_script_sha256": _sha256(ANALYSIS),
    }
    for key, observed in hash_checks.items():
        if metadata[key] != observed:
            raise AssertionError(f"metadata hash mismatch: {key}")
    if metadata["selected_series_count"] != len(EXPECTED_SITES):
        raise AssertionError("selected-series denominator is not nine")

    k0 = pd.read_csv(OUTPUT / "usgs_external_k0_detail.csv")
    k1 = pd.read_csv(OUTPUT / "usgs_external_k1_detail.csv")
    blocks = pd.read_csv(OUTPUT / "usgs_external_station_block_summary.csv")
    cadence = pd.read_csv(OUTPUT / "usgs_external_cadence_summary.csv")
    audit = pd.read_csv(OUTPUT / "usgs_external_data_audit.csv")
    if tuple(audit["monitoring_location_id"]) != EXPECTED_SITES:
        raise AssertionError("data-audit site order or denominator changed")
    if audit["duplicate_reference_timestamps"].sum() != 0:
        raise AssertionError("duplicate reference timestamps were not stopped")
    if audit["qualifying_7day_blocks"].sum() != len(blocks):
        raise AssertionError("block attrition accounting is inconsistent")

    expected_case_count = 0
    for row in blocks.itertuples(index=False):
        for phase_count in EXPECTED_PHASES.values():
            expected_case_count += phase_count
    if len(k0) != expected_case_count:
        raise AssertionError("all-phase case denominator is incomplete")
    phase_counts = (
        k0.groupby(
            ["monitoring_location_id", "block_index", "cadence_hours"]
        )["phase_index"]
        .nunique()
        .reset_index(name="observed")
    )
    phase_counts["expected"] = phase_counts["cadence_hours"].map(EXPECTED_PHASES)
    if not phase_counts["observed"].eq(phase_counts["expected"]).all():
        raise AssertionError("a block-cadence group lacks frozen phases")

    compatible_k0 = k0.loc[
        k0["status"].eq("feasible") & k0["reference_L_compatible"].eq(True)
    ].copy()
    if not compatible_k0["occupation_reference_covered"].all():
        raise AssertionError("k=0 occupation containment failure")
    if not compatible_k0["deficit_reference_covered"].all():
        raise AssertionError("k=0 deficit containment failure")
    if not compatible_k0["occupation_width_fraction"].between(
        -TOLERANCE, 1.0 + TOLERANCE
    ).all():
        raise AssertionError("k=0 occupation width outside physical range")
    if not compatible_k0[
        "deficit_width_fraction_nonnegative_H_times_T"
    ].between(-TOLERANCE, 1.0 + TOLERANCE).all():
        raise AssertionError("k=0 deficit width outside physical range")
    if not (
        compatible_k0["deficit_upper_nonnegative_mgL_hours"]
        <= compatible_k0["deficit_upper_unconstrained_mgL_hours"] + TOLERANCE
    ).all():
        raise AssertionError("state floor widened a deficit upper endpoint")

    if _case_keys(k1) != _case_keys(compatible_k0):
        raise AssertionError("k=1 cases differ from compatible k=0 cases")
    if not k1["k0_nested_in_k1"].all():
        raise AssertionError("k=0 is not nested within k=1")
    if not k1["occupation_reference_covered_k1"].all():
        raise AssertionError("k=1 occupation containment failure")
    if not k1["deficit_reference_covered_k1"].all():
        raise AssertionError("k=1 deficit containment failure")
    if not k1["common_anchor_condition_M_gt_2k"].all():
        raise AssertionError("M > 2k failed")
    count_columns = [
        column for column in k1.columns if column.endswith("_deleted_count")
    ]
    if not count_columns or (k1[count_columns] > 1).any().any():
        raise AssertionError("an endpoint witness exceeded k=1")
    if not k1["occupation_width_fraction_k1"].between(
        -TOLERANCE, 1.0 + TOLERANCE
    ).all():
        raise AssertionError("k=1 occupation width outside physical range")
    if not k1["deficit_width_fraction_k1"].between(
        -TOLERANCE, 1.0 + TOLERANCE
    ).all():
        raise AssertionError("k=1 deficit width outside physical range")

    if set(cadence["cadence_hours"]) != set(EXPECTED_PHASES):
        raise AssertionError("cadence summary is incomplete")
    failure_columns = [
        "containment_failures",
        "deficit_containment_failures",
        "k1_nesting_failures",
    ]
    if cadence[failure_columns].sum().sum() != 0:
        raise AssertionError("summary reports a deterministic failure")
    compatible_blocks = int(blocks["reference_L_compatible"].sum())
    if not cadence["compatible_station_blocks"].eq(compatible_blocks).all():
        raise AssertionError("cadence-level compatible block denominator changed")

    receipt = {
        "status": "passed",
        "selected_series_count": len(EXPECTED_SITES),
        "series_with_qualifying_blocks": int(
            audit["qualifying_7day_blocks"].gt(0).sum()
        ),
        "qualifying_block_count": len(blocks),
        "compatible_block_count": compatible_blocks,
        "all_phase_case_count": len(k0),
        "compatible_k0_case_count": len(compatible_k0),
        "fixed_k1_case_count": len(k1),
        "phase_geometry_groups": len(phase_counts),
        "containment_failures": 0,
        "nesting_failures": 0,
        "source_hashes_checked": 2 * len(EXPECTED_SITES),
        "frozen_file_hashes": hash_checks,
        "output_hashes": {
            path.name: _sha256(path)
            for path in sorted(OUTPUT.iterdir())
            if path.is_file()
            and path.name != "USGS_EXTERNAL_REPLICATION_VERIFICATION.json"
        },
    }
    receipt_path = OUTPUT / "USGS_EXTERNAL_REPLICATION_VERIFICATION.json"
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return receipt_path


if __name__ == "__main__":
    print(verify())
