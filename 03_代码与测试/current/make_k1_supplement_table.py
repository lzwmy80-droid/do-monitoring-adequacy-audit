# -*- coding: utf-8 -*-
"""Freeze Supplementary Table S2 from the fixed-k=1 sensitivity outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SOURCE = (
    PROJECT / "05_既有结果" / "current_maracoos_7day_fixed_k1_sensitivity"
)
OUTPUT = PROJECT / "09_论文写作输出"
TABLES = OUTPUT / "tables"
SUMMARY = SOURCE / "fixed_k1_summary.csv"
DETAIL = SOURCE / "fixed_k1_detail.csv"
METADATA = SOURCE / "fixed_k1_metadata.json"
ANALYSIS_SCRIPT = HERE / "run_maracoos_fixed_k1_sensitivity.py"
FREEZE_NOTE = OUTPUT / "K1_SENSITIVITY_FREEZE_NOTE_20260727.md"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _range_text(low: float, high: float, digits: int) -> str:
    return f"{low:.{digits}f}–{high:.{digits}f}"


def build() -> Path:
    TABLES.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(SUMMARY)
    metadata = json.loads(METADATA.read_text(encoding="utf-8"))

    if int(summary["cases"].sum()) != 240:
        raise ValueError("Supplementary Table S2 expects exactly 240 cases")
    if metadata["containment_failures"] != 0 or metadata["nesting_failures"] != 0:
        raise ValueError("fixed-k=1 diagnostics contain failures")
    if not all(
        metadata[name]
        for name in (
            "all_minimum_replacements_equal_zero",
            "all_cases_satisfy_M_gt_2k",
            "all_physical_width_checks_passed",
        )
    ):
        raise ValueError("a fixed-k=1 audit assertion failed")

    table = pd.DataFrame(
        {
            "cadence_hours": summary["cadence_hours"],
            "cases": summary["cases"].astype(int),
            "sparse_reports_min": summary["min_sparse_reports"].astype(int),
            "sparse_reports_max": summary["max_sparse_reports"].astype(int),
            "permitted_fraction_percent_min": (
                100 * summary["min_permitted_replacement_fraction"]
            ),
            "permitted_fraction_percent_max": (
                100 * summary["max_permitted_replacement_fraction"]
            ),
            "median_occupation_width_k0": (
                summary["median_occupation_width_fraction_k0"]
            ),
            "median_occupation_width_k1": (
                summary["median_occupation_width_fraction_k1"]
            ),
            "median_paired_occupation_increase": (
                summary["median_occupation_width_fraction_increase"]
            ),
            "median_deficit_width_k0": (
                summary["median_deficit_width_fraction_k0"]
            ),
            "median_deficit_width_k1": (
                summary["median_deficit_width_fraction_k1"]
            ),
            "median_paired_deficit_increase": (
                summary["median_deficit_width_fraction_increase"]
            ),
            "cases_any_endpoint_uses_deletion": (
                summary["cases_any_endpoint_uses_deletion"].astype(int)
            ),
            "diagnostic_failures": (
                summary["occupation_containment_failures_k1"]
                + summary["deficit_containment_failures_k1"]
                + summary["nesting_failures"]
            ).astype(int),
        }
    )

    csv_path = TABLES / "Table_S2_fixed_k1_sensitivity.csv"
    md_path = TABLES / "Table_S2_fixed_k1_sensitivity.md"
    tex_path = TABLES / "Table_S2_fixed_k1_sensitivity.tex"
    caption_path = OUTPUT / "SUPPLEMENT_TABLE_S2_CAPTION.md"
    table.to_csv(csv_path, index=False)

    md_lines = [
        r"| Cadence (h) | Cases | Reports \(M\) | Permitted \(100/M\) (%) | "
        r"Occ. width \(k=0\) | Occ. width \(k=1\) | Paired Δ | "
        r"Deficit width \(k=0\) | Deficit width \(k=1\) | Paired Δ | "
        "Witness uses permitted deletion | Containment/nesting failures |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in table.itertuples(index=False):
        md_lines.append(
            f"| {row.cadence_hours:g} | {row.cases} | "
            f"{row.sparse_reports_min}–{row.sparse_reports_max} | "
            f"{_range_text(row.permitted_fraction_percent_min, row.permitted_fraction_percent_max, 3)} | "
            f"{row.median_occupation_width_k0:.3f} | "
            f"{row.median_occupation_width_k1:.3f} | "
            f"{row.median_paired_occupation_increase:.3f} | "
            f"{row.median_deficit_width_k0:.3f} | "
            f"{row.median_deficit_width_k1:.3f} | "
            f"{row.median_paired_deficit_increase:.3f} | "
            f"{row.cases_any_endpoint_uses_deletion}/{row.cases} | "
            f"{row.diagnostic_failures} |"
        )
    md_lines.extend(
        [
            "",
            "Note: failures are the combined deterministic 15-min "
            "reference-containment and k=0-within-k=1 nesting failures. "
            "The pooled cadence-phase cases are not independent windows. "
            "A witness deletion is an endpoint certificate, not a verified "
            "sensor error.",
        ]
    )
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    tex_lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Fixed-count sensitivity to at most one arbitrary replacement.}",
        r"\label{tab:fixed_k1_sensitivity}",
        r"\resizebox{\textwidth}{!}{%",
        r"\begin{tabular}{rrrrrrrrrrrr}",
        r"\toprule",
        r"Cadence & Cases & $M$ & $100/M$ (\%) & "
        r"$W_O(0)$ & $W_O(1)$ & paired $\Delta_O$ & "
        r"$W_D(0)$ & $W_D(1)$ & paired $\Delta_D$ & "
        r"witness uses permitted deletion & check failures \\",
        r"\midrule",
    ]
    for row in table.itertuples(index=False):
        tex_lines.append(
            f"{row.cadence_hours:g} & {row.cases} & "
            f"{row.sparse_reports_min}--{row.sparse_reports_max} & "
            f"{row.permitted_fraction_percent_min:.3f}--"
            f"{row.permitted_fraction_percent_max:.3f} & "
            f"{row.median_occupation_width_k0:.3f} & "
            f"{row.median_occupation_width_k1:.3f} & "
            f"{row.median_paired_occupation_increase:.3f} & "
            f"{row.median_deficit_width_k0:.3f} & "
            f"{row.median_deficit_width_k1:.3f} & "
            f"{row.median_paired_deficit_increase:.3f} & "
            f"{row.cases_any_endpoint_uses_deletion}/{row.cases} & "
            f"{row.diagnostic_failures} \\\\"
        )
    tex_lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}%",
            r"}",
            r"\begin{minipage}{0.98\textwidth}",
            r"\footnotesize Notes: Check failures combine deterministic "
            r"15-min reference containment and $k=0$-within-$k=1$ nesting. "
            r"Pooled cadence--phase cases are not independent windows. "
            r"A witness deletion is an endpoint certificate, not a verified "
            r"sensor error.",
            r"\end{minipage}",
            r"\end{table*}",
        ]
    )
    tex_path.write_text("\n".join(tex_lines) + "\n", encoding="utf-8")

    caption_path.write_text(
        "\n".join(
            [
                "# Supplementary Table S2 caption",
                "",
                "**Fixed-count sensitivity to at most one arbitrary "
                "replacement.** All and only the 240 seven-day cases that "
                r"were feasible and reference-\(L\)-compatible in the frozen "
                r"\(k=0\) analysis were re-evaluated with the externally fixed "
                r"budget \(k=1\), \(H=5\) mg/L, the calibration-only q99.9 "
                r"slope label, and \(x(t)\ge0\) mg/L. Occupation and deficit "
                r"widths are normalized by \(T\) and \(HT\), respectively. "
                r"Paired Δ is the median within-case increase from \(k=0\) to "
                r"\(k=1\). “Witness uses permitted deletion” counts cases in "
                "which at least one endpoint certificate uses that option; "
                "it is not an error count. Failures combine deterministic "
                "15-min "
                r"reference-containment and \(k=0\)-within-\(k=1\) nesting "
                "checks. The fixed count implies a cadence-dependent "
                "permitted fraction, and pooled phase cases are not "
                "independent windows.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    assets = (csv_path, md_path, tex_path, caption_path)
    sources = (
        SUMMARY,
        DETAIL,
        METADATA,
        ANALYSIS_SCRIPT,
        FREEZE_NOTE,
        Path(__file__).resolve(),
    )
    manifest = OUTPUT / "K1_SENSITIVITY_ASSET_MANIFEST_20260727.json"
    payload = {
        "status": "fixed_k1_supplementary_table_frozen",
        "scope": (
            "paired sensitivity on the pre-frozen 240 compatible cases; "
            "fixed count k=1, not an estimated error count or rate"
        ),
        "sources": {str(path): _sha256(path) for path in sources},
        "assets": {str(path): _sha256(path) for path in assets},
    }
    manifest.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


if __name__ == "__main__":
    print(build())
