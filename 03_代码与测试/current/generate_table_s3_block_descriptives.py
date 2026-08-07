#!/usr/bin/env python3
"""Generate Supplementary Table S3 from the frozen seven-day detail file.

The table changes the reporting unit from overlapping cadence–phase cases to
station-block × cadence summaries.  It is descriptive only: no independence
assumption, confidence interval, or population inference is introduced.
"""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "05_既有结果"
    / "current_maracoos_7day_nonnegative_table2"
    / "table2_nonnegative_detail.csv"
)
OUT_DIR = ROOT / "09_论文写作输出" / "tables"
CSV_OUT = OUT_DIR / "Table_S3_block_level_widths.csv"
MD_OUT = OUT_DIR / "Table_S3_block_level_widths.md"
TEX_OUT = OUT_DIR / "Table_S3_block_level_widths.tex"
CAPTION_OUT = ROOT / "09_论文写作输出" / "SUPPLEMENT_TABLE_S3_CAPTION.md"
MANIFEST_OUT = (
    ROOT / "09_论文写作输出" / "TABLE_S3_ASSET_MANIFEST_20260727.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(character, character) for character in text)


def main() -> None:
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle))

    compatible = [
        row
        for row in source_rows
        if row["status"] == "feasible"
        and row["reference_L_compatible"].strip().lower() == "true"
    ]
    if len(compatible) != 240:
        raise RuntimeError(f"Expected 240 compatible cases, found {len(compatible)}")

    grouped: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in compatible:
        key = (
            row["dataset_id"],
            row["block_index"],
            row["block_start"],
            row["block_end"],
            row["cadence_hours"],
        )
        grouped[key].append(row)

    station_blocks = {(key[0], key[1]) for key in grouped}
    stations = {key[0] for key in grouped}
    if len(station_blocks) != 8 or len(stations) != 2 or len(grouped) != 32:
        raise RuntimeError(
            "Frozen cohort mismatch: "
            f"stations={len(stations)}, station_blocks={len(station_blocks)}, "
            f"block-cadence groups={len(grouped)}"
        )

    rows: list[dict[str, str]] = []
    for key in sorted(
        grouped,
        key=lambda value: (
            value[0],
            int(value[1]),
            float(value[4]),
        ),
    ):
        dataset_id, block_index, block_start, block_end, cadence = key
        cases = grouped[key]
        rows.append(
            {
                "station": dataset_id.removeprefix("mddnr_").replace("_", " "),
                "block_index": block_index,
                "block_start_utc": block_start,
                "block_end_utc": block_end,
                "cadence_hours": f"{float(cadence):g}",
                "compatible_phase_cases_n": str(len(cases)),
                "median_occupation_width_over_T": (
                    f"{median([float(row['occupation_width_fraction']) for row in cases]):.3f}"
                ),
                "median_nonnegative_deficit_width_over_HT": (
                    f"{median([float(row['deficit_width_fraction_nonnegative_H_times_T']) for row in cases]):.3f}"
                ),
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with CSV_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    display_headers = [
        "Station",
        "Block",
        "Start (UTC)",
        "End (UTC)",
        "Cadence (h)",
        "Compatible phase cases, n",
        "Median occupation width / T",
        "Median nonnegative deficit width / (HT)",
    ]
    display_rows = [
        [
            row["station"],
            row["block_index"],
            row["block_start_utc"],
            row["block_end_utc"],
            row["cadence_hours"],
            row["compatible_phase_cases_n"],
            row["median_occupation_width_over_T"],
            row["median_nonnegative_deficit_width_over_HT"],
        ]
        for row in rows
    ]

    markdown = [
        "# Supplementary Table S3",
        "",
        "**Station-block-level descriptive widths for the seven-day sampling "
        "baseline.** Each row is one compatible station-block × cadence group. "
        "Values are medians over overlapping phase offsets within that group; "
        "phase cases are not independent windows. No inferential interval or "
        "population weighting is applied.",
        "",
        "| " + " | ".join(display_headers) + " |",
        "|" + "|".join(["---"] * len(display_headers)) + "|",
    ]
    markdown.extend("| " + " | ".join(row) + " |" for row in display_rows)
    markdown.append("")
    MD_OUT.write_text("\n".join(markdown), encoding="utf-8")

    column_spec = "llp{2.5cm}p{2.5cm}rrrr"
    tex_lines = [
        r"\begin{longtable}{" + column_spec + "}",
        r"\caption{Station-block-level descriptive widths for the seven-day "
        r"sampling baseline. Values are medians over overlapping phase offsets; "
        r"phase cases are not independent windows.}\\",
        r"\toprule",
        " & ".join(latex_escape(value) for value in display_headers) + r" \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        " & ".join(latex_escape(value) for value in display_headers) + r" \\",
        r"\midrule",
        r"\endhead",
    ]
    tex_lines.extend(
        " & ".join(latex_escape(value) for value in row) + r" \\" for row in display_rows
    )
    tex_lines.extend([r"\bottomrule", r"\end{longtable}", ""])
    TEX_OUT.write_text("\n".join(tex_lines), encoding="utf-8")

    caption = r"""# Supplementary Table S3 caption

**Station-block-level descriptive widths for the seven-day** \(k=0\) **sampling
baseline.** Each row is one reference-\(L\)-compatible station-block × cadence
group. Occupation and nonnegative cumulative-deficit widths are normalized by
\(T\) and \(HT\), respectively, and summarized by the median over overlapping
phase offsets within the group. Phase cases are not independent windows; the
table is descriptive and supplies neither confidence intervals nor
population-weighted inference.
"""
    CAPTION_OUT.write_text(caption, encoding="utf-8")

    manifest = {
        "asset": "Supplementary Table S3",
        "frozen_input": {
            "path": str(SOURCE.relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(SOURCE),
            "compatible_cases": len(compatible),
            "stations": len(stations),
            "station_blocks": len(station_blocks),
            "station_block_cadence_rows": len(rows),
        },
        "generator": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "outputs": [
            {
                "path": str(path.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(path),
            }
            for path in (CSV_OUT, MD_OUT, TEX_OUT, CAPTION_OUT)
        ],
    }
    MANIFEST_OUT.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "compatible_cases": len(compatible),
                "stations": len(stations),
                "station_blocks": len(station_blocks),
                "rows": len(rows),
                "csv": str(CSV_OUT),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
