# -*- coding: utf-8 -*-
"""Build the frozen Figure 2 and Table 2 manuscript assets."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from pathlib import Path

import cairosvg


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
TRILEMMA = (
    PROJECT
    / "05_既有结果"
    / "current_maracoos_24h_nonnegative_trilemma"
    / "trilemma_nonnegative_detail.csv"
)
TABLE2 = (
    PROJECT
    / "05_既有结果"
    / "current_maracoos_7day_nonnegative_table2"
    / "table2_nonnegative_summary.csv"
)
TABLE2_DETAIL = TABLE2.with_name("table2_nonnegative_detail.csv")
TRILEMMA_METADATA = TRILEMMA.with_name("trilemma_nonnegative_metadata.json")
TABLE2_METADATA = TABLE2.with_name("table2_nonnegative_metadata.json")
RUNTIME = HERE / "robust_exact_current.py"
OUTPUT = PROJECT / "09_论文写作输出"
FIGURES = OUTPUT / "figures"
TABLES = OUTPUT / "tables"

WIDTH = 1780
HEIGHT = 560
PANEL_WIDTH = 470
PANEL_HEIGHT = 365
PANEL_TOP = 105
PANEL_LEFTS = (75, 625, 1175)
X_MIN = math.log2(4.0)
X_MAX = math.log2(64.0)
L_TICKS = (4.0, 8.0, 13.20396, 16.0, 25.52, 32.0, 64.0)

STATIONS = (
    "mddnr_Bishopville_Prong",
    "mddnr_Camp_Tockwogh",
    "mddnr_Greys_Creek",
    "mddnr_Harris_Creek_Upstream",
)
STATION_LABELS = {
    "mddnr_Bishopville_Prong": "Bishopville Prong",
    "mddnr_Camp_Tockwogh": "Camp Tockwogh",
    "mddnr_Greys_Creek": "Greys Creek",
    "mddnr_Harris_Creek_Upstream": "Harris Creek Upstream",
}
COLORS = {
    "mddnr_Bishopville_Prong": "#0072B2",
    "mddnr_Camp_Tockwogh": "#E69F00",
    "mddnr_Greys_Creek": "#009E73",
    "mddnr_Harris_Creek_Upstream": "#CC79A7",
}
MARKERS = {
    "mddnr_Bishopville_Prong": "circle",
    "mddnr_Camp_Tockwogh": "square",
    "mddnr_Greys_Creek": "triangle",
    "mddnr_Harris_Creek_Upstream": "diamond",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _x(value: float, left: float) -> float:
    fraction = (math.log2(value) - X_MIN) / (X_MAX - X_MIN)
    return left + fraction * PANEL_WIDTH


def _y(value: float, maximum: float) -> float:
    return PANEL_TOP + PANEL_HEIGHT * (1.0 - value / maximum)


def _text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 20,
    anchor: str = "middle",
    weight: str = "normal",
    fill: str = "#222222",
    rotate: float | None = None,
) -> str:
    transform = (
        f' transform="rotate({rotate:g} {x:.2f} {y:.2f})"'
        if rotate is not None
        else ""
    )
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="Arial,Helvetica,sans-serif" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}"{transform}>'
        f"{html.escape(value)}</text>"
    )


def _marker(
    x: float,
    y: float,
    kind: str,
    color: str,
    *,
    size: float = 7.0,
    opacity: float = 1.0,
) -> str:
    common = (
        f'fill="{color}" stroke="#ffffff" stroke-width="1.5" '
        f'opacity="{opacity:.3f}"'
    )
    if kind == "circle":
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" {common}/>'
    if kind == "square":
        return (
            f'<rect x="{x-size:.2f}" y="{y-size:.2f}" '
            f'width="{2*size:.2f}" height="{2*size:.2f}" {common}/>'
        )
    if kind == "triangle":
        points = (
            f"{x:.2f},{y-size-1:.2f} "
            f"{x-size-1:.2f},{y+size:.2f} "
            f"{x+size+1:.2f},{y+size:.2f}"
        )
        return f'<polygon points="{points}" {common}/>'
    points = (
        f"{x:.2f},{y-size-1:.2f} "
        f"{x+size+1:.2f},{y:.2f} "
        f"{x:.2f},{y+size+1:.2f} "
        f"{x-size-1:.2f},{y:.2f}"
    )
    return f'<polygon points="{points}" {common}/>'


def _panel_axes(
    left: float,
    letter: str,
    y_label: str,
    y_max: float,
    y_ticks: tuple[float, ...],
) -> list[str]:
    elements = [
        f'<rect x="{left:.2f}" y="{PANEL_TOP:.2f}" '
        f'width="{PANEL_WIDTH}" height="{PANEL_HEIGHT}" '
        'fill="#ffffff" stroke="#333333" stroke-width="1.5"/>',
        _text(left - 52, PANEL_TOP - 24, letter, size=27, anchor="start", weight="bold"),
    ]
    for tick in y_ticks:
        position = _y(tick, y_max)
        elements.append(
            f'<line x1="{left:.2f}" y1="{position:.2f}" '
            f'x2="{left+PANEL_WIDTH:.2f}" y2="{position:.2f}" '
            'stroke="#d9d9d9" stroke-width="1"/>'
        )
        label = f"{tick:.2f}".rstrip("0").rstrip(".")
        elements.append(
            _text(left - 11, position + 7, label, size=17, anchor="end")
        )
    for tick in L_TICKS:
        position = _x(tick, left)
        calibration = math.isclose(tick, 13.20396) or math.isclose(tick, 25.52)
        elements.append(
            f'<line x1="{position:.2f}" y1="{PANEL_TOP:.2f}" '
            f'x2="{position:.2f}" y2="{PANEL_TOP+PANEL_HEIGHT:.2f}" '
            f'stroke="{"#b0b0b0" if calibration else "#eeeeee"}" '
            f'stroke-width="{"1.5" if calibration else "1"}" '
            f'stroke-dasharray="{"5,4" if calibration else "none"}"/>'
        )
        label = (
            "13.2"
            if math.isclose(tick, 13.20396)
            else "25.5"
            if math.isclose(tick, 25.52)
            else f"{tick:g}"
        )
        elements.append(
            _text(position, PANEL_TOP + PANEL_HEIGHT + 26, label, size=16)
        )
    elements.extend(
        [
            _text(
                left + PANEL_WIDTH / 2,
                PANEL_TOP + PANEL_HEIGHT + 58,
                "Lipschitz bound L (mg/L/h; log2 scale)",
                size=18,
            ),
            _text(
                left - 55,
                PANEL_TOP + PANEL_HEIGHT / 2,
                y_label,
                size=18,
                rotate=-90,
            ),
        ]
    )
    return elements


def _line(
    points: list[tuple[float, float]],
    color: str,
    *,
    dashed: bool = False,
    opacity: float = 1.0,
    width: float = 3.0,
) -> str:
    path = " ".join(
        ("M" if index == 0 else "L") + f" {x:.2f} {y:.2f}"
        for index, (x, y) in enumerate(points)
    )
    dash = ' stroke-dasharray="9,7"' if dashed else ""
    return (
        f'<path d="{path}" fill="none" stroke="{color}" '
        f'stroke-width="{width:.2f}" opacity="{opacity:.3f}"{dash}/>'
    )


def build_figure() -> Path:
    rows = _read_rows(TRILEMMA)
    if len(rows) != 28:
        raise ValueError(f"expected 28 frozen trilemma rows, found {len(rows)}")
    grouped: dict[str, list[dict[str, float | str]]] = {
        station: [] for station in STATIONS
    }
    for row in rows:
        station = row["dataset_id"]
        if station not in grouped:
            raise ValueError(f"unexpected validation station: {station}")
        grouped[station].append(
            {
                "L": float(row["L_mgL_per_hour"]),
                "k": float(row["k_min_fraction"]),
                "occupation": float(row["occupation_width_fraction"]),
                "deficit_unconstrained": float(
                    row[
                        "deficit_width_fraction_unconstrained_H_times_T"
                    ]
                ),
                "deficit_nonnegative": float(
                    row[
                        "deficit_width_fraction_nonnegative_H_times_T"
                    ]
                ),
            }
        )
    for station_rows in grouped.values():
        station_rows.sort(key=lambda row: float(row["L"]))
        if len(station_rows) != 7:
            raise ValueError("each station must have seven frozen L values")

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" '
        f'height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    legend_x = 90
    for station in STATIONS:
        color = COLORS[station]
        elements.append(
            _marker(
                legend_x,
                31,
                MARKERS[station],
                color,
                size=6.5,
            )
        )
        elements.append(
            _text(
                legend_x + 14,
                38,
                STATION_LABELS[station],
                size=18,
                anchor="start",
            )
        )
        legend_x += {
            "mddnr_Bishopville_Prong": 245,
            "mddnr_Camp_Tockwogh": 225,
            "mddnr_Greys_Creek": 185,
            "mddnr_Harris_Creek_Upstream": 285,
        }[station]
    style_x = 1165
    elements.extend(
        [
            f'<line x1="{style_x}" y1="27" x2="{style_x+42}" y2="27" '
            'stroke="#444444" stroke-width="3"/>',
            _text(style_x + 50, 34, "DO ≥ 0", size=17, anchor="start"),
            f'<line x1="{style_x+155}" y1="27" x2="{style_x+197}" y2="27" '
            'stroke="#777777" stroke-width="3" stroke-dasharray="9,7" '
            'opacity="0.65"/>',
            _text(
                style_x + 205,
                34,
                "unbounded-state audit",
                size=17,
                anchor="start",
            ),
        ]
    )

    elements.extend(
        _panel_axes(
            PANEL_LEFTS[0],
            "A",
            "Minimum replacement fraction, k_min / M",
            0.22,
            (0.0, 0.05, 0.10, 0.15, 0.20),
        )
    )
    elements.extend(
        _panel_axes(
            PANEL_LEFTS[1],
            "B",
            "Occupation endpoint width / T",
            0.90,
            (0.0, 0.2, 0.4, 0.6, 0.8),
        )
    )
    elements.extend(
        _panel_axes(
            PANEL_LEFTS[2],
            "C",
            "Deficit endpoint width / (H T)",
            1.35,
            (0.0, 0.25, 0.50, 0.75, 1.00, 1.25),
        )
    )
    physical_y = _y(1.0, 1.35)
    elements.extend(
        [
            f'<line x1="{PANEL_LEFTS[2]:.2f}" y1="{physical_y:.2f}" '
            f'x2="{PANEL_LEFTS[2]+PANEL_WIDTH:.2f}" y2="{physical_y:.2f}" '
            'stroke="#D55E00" stroke-width="2" stroke-dasharray="4,4"/>',
            _text(
                PANEL_LEFTS[2] + PANEL_WIDTH - 6,
                physical_y - 8,
                "physical maximum",
                size=15,
                anchor="end",
                fill="#B34700",
            ),
        ]
    )
    for label, value, offset, anchor in (
        ("cal. q99.9", 13.20396, -7.0, "end"),
        ("cal. max", 25.52, 7.0, "start"),
    ):
        elements.append(
            _text(
                _x(value, PANEL_LEFTS[0]) + offset,
                PANEL_TOP + 20,
                label,
                size=13,
                anchor=anchor,
                fill="#666666",
            )
        )

    for station in STATIONS:
        station_rows = grouped[station]
        color = COLORS[station]
        marker = MARKERS[station]
        series = (
            ("k", PANEL_LEFTS[0], 0.22),
            ("occupation", PANEL_LEFTS[1], 0.90),
            ("deficit_nonnegative", PANEL_LEFTS[2], 1.35),
        )
        for key, left, maximum in series:
            points = [
                (_x(float(row["L"]), left), _y(float(row[key]), maximum))
                for row in station_rows
            ]
            elements.append(_line(points, color))
            elements.extend(
                _marker(x, y, marker, color) for x, y in points
            )
        audit_points = [
            (
                _x(float(row["L"]), PANEL_LEFTS[2]),
                _y(float(row["deficit_unconstrained"]), 1.35),
            )
            for row in station_rows
        ]
        elements.append(
            _line(
                audit_points,
                color,
                dashed=True,
                opacity=0.43,
                width=2.4,
            )
        )

    elements.append("</svg>")
    FIGURES.mkdir(parents=True, exist_ok=True)
    svg_path = FIGURES / "Figure_2_L_k_width_nonnegative.svg"
    svg_path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    cairosvg.svg2png(
        bytestring=svg_path.read_bytes(),
        write_to=str(FIGURES / "Figure_2_L_k_width_nonnegative.png"),
        output_width=2 * WIDTH,
        output_height=2 * HEIGHT,
    )
    cairosvg.svg2pdf(
        bytestring=svg_path.read_bytes(),
        write_to=str(FIGURES / "Figure_2_L_k_width_nonnegative.pdf"),
    )
    return svg_path


def build_table() -> tuple[Path, Path, Path]:
    rows = _read_rows(TABLE2)
    if len(rows) != 4:
        raise ValueError(f"expected four frozen Table 2 rows, found {len(rows)}")
    rows.sort(key=lambda row: float(row["cadence_hours"]))
    detail_rows = _read_rows(TABLE2_DETAIL)
    if len(detail_rows) != 480:
        raise ValueError(
            f"expected 480 frozen Table 2 detail rows, found {len(detail_rows)}"
        )
    counts: dict[float, dict[str, int]] = {}
    for row in detail_rows:
        cadence = float(row["cadence_hours"])
        cadence_counts = counts.setdefault(
            cadence,
            {"total": 0, "feasible": 0, "compatible": 0},
        )
        cadence_counts["total"] += 1
        feasible = row["status"] == "feasible"
        if feasible:
            cadence_counts["feasible"] += 1
        compatible = row["reference_L_compatible"].strip().lower() == "true"
        if feasible and compatible:
            cadence_counts["compatible"] += 1
    for row in rows:
        cadence = float(row["cadence_hours"])
        expected = int(float(row["compatible_feasible_cases"]))
        if counts[cadence]["compatible"] != expected:
            raise ValueError(
                f"compatible case mismatch at cadence {cadence:g}: "
                f"{counts[cadence]['compatible']} != {expected}"
            )
    TABLES.mkdir(parents=True, exist_ok=True)
    csv_path = TABLES / "Table_2_7day_nonnegative.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Cadence (h)",
                "Total pooled cases",
                "Algorithm-feasible cases",
                "Reference-compatible cases",
                "Median occupation width / T",
                "Median deficit width / (H T)",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    f"{float(row['cadence_hours']):g}",
                    counts[float(row["cadence_hours"])]["total"],
                    counts[float(row["cadence_hours"])]["feasible"],
                    counts[float(row["cadence_hours"])]["compatible"],
                    f"{float(row['median_occupation_width_fraction']):.3f}",
                    f"{float(row['median_deficit_width_fraction_nonnegative']):.3f}",
                ]
            )

    md_path = TABLES / "Table_2_7day_nonnegative.md"
    md_lines = [
        "| Cadence (h) | Total pooled cases | Algorithm-feasible cases | "
        "Reference-compatible cases | "
        "Median occupation width / T | Median deficit width / (H T) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        cadence = float(row["cadence_hours"])
        md_lines.append(
            f"| {cadence:g} | "
            f"{counts[cadence]['total']} | "
            f"{counts[cadence]['feasible']} | "
            f"{counts[cadence]['compatible']} | "
            f"{float(row['median_occupation_width_fraction']):.3f} | "
            f"{float(row['median_deficit_width_fraction_nonnegative']):.3f} |"
        )
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    tex_path = TABLES / "Table_2_7day_nonnegative.tex"
    tex_lines = [
        r"\begin{tabular}{rrrrrr}",
        r"\toprule",
        r"Cadence (h) & Total & Feasible & Compatible & "
        r"Occupation width$/T$ & Deficit width$/(HT)$ \\",
        r"\midrule",
    ]
    for row in rows:
        cadence = float(row["cadence_hours"])
        tex_lines.append(
            f"{cadence:g} & "
            f"{counts[cadence]['total']} & "
            f"{counts[cadence]['feasible']} & "
            f"{counts[cadence]['compatible']} & "
            f"{float(row['median_occupation_width_fraction']):.3f} & "
            f"{float(row['median_deficit_width_fraction_nonnegative']):.3f} "
            r"\\"
        )
    tex_lines.extend([r"\bottomrule", r"\end{tabular}"])
    tex_path.write_text("\n".join(tex_lines) + "\n", encoding="utf-8")
    return csv_path, md_path, tex_path


def write_caption_and_manifest(
    figure_path: Path,
    table_paths: tuple[Path, Path, Path],
) -> Path:
    captions = OUTPUT / "FROZEN_FIGURE_TABLE_CAPTIONS.md"
    captions.write_text(
        "\n".join(
            [
                "# Frozen captions",
                "",
                "## Figure 2",
                "",
                "**Conditional observability trade-off under sparse record-value "
                "corruption.** Each curve uses the earliest complete 24-h "
                "15-min block at one frozen holdout station and threshold "
                "$H=5$ mg L$^{-1}$. For each candidate Lipschitz bound $L$, "
                "$k$ is set to the minimum discrete replacement budget needed "
                "for compatibility; this is a data-dependent feasibility "
                "witness, not a known error count. "
                "(A) Required replacement fraction; (B) sharp occupation "
                "endpoint width; (C) sharp cumulative-deficit endpoint width. "
                "For the nonnegative reports and $H>0$ used here, panels A and "
                "B are shared by the constrained and unbounded-state "
                "formulations. In panel C, solid curves impose "
                "$x(t)\\geq0$ mg L$^{-1}$, whereas faint dashed curves show "
                "the unbounded-state audit. The empirical surface is a "
                "conditional stress test and does not validate $L$ or label "
                "deleted reports as errors.",
                "",
                "## Table 2",
                "",
                "**Seven-day $k=0$ sampling baseline under the physical "
                "dissolved-oxygen state constraint.** Width entries are medians "
                "over all sampling "
                "phases in the reference-compatible feasible subset; the three "
                "case-count columns expose the preceding filters. Frozen "
                "holdout blocks contribute to the width summary only when their "
                "15-min operational "
                "reference is compatible with the calibration-only q99.9 "
                "slope bound ($L=13.204$ mg L$^{-1}$ h$^{-1}$), at "
                "$H=5$ mg L$^{-1}$ and $k=0$. The latent path is constrained "
                "by $x(t)\\geq0$ mg L$^{-1}$. The compatible summary comprises "
                "8/16 station-window blocks, 2/4 holdout stations, and "
                "240/480 pooled cadence–phase cases. The operational reference "
                "is not continuous truth, and the empirical slope is not asserted "
                "to be a physical bound.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    assets = [
        figure_path,
        FIGURES / "Figure_2_L_k_width_nonnegative.png",
        FIGURES / "Figure_2_L_k_width_nonnegative.pdf",
        *table_paths,
        captions,
    ]
    manifest = {
        "status": "frozen_assets_built",
        "source_files": {
            str(TRILEMMA): _sha256(TRILEMMA),
            str(TRILEMMA_METADATA): _sha256(TRILEMMA_METADATA),
            str(TABLE2): _sha256(TABLE2),
            str(TABLE2_DETAIL): _sha256(TABLE2_DETAIL),
            str(TABLE2_METADATA): _sha256(TABLE2_METADATA),
            str(RUNTIME): _sha256(RUNTIME),
            str(HERE / "make_frozen_paper_assets.py"): _sha256(
                HERE / "make_frozen_paper_assets.py"
            ),
        },
        "assets": {str(path): _sha256(path) for path in assets},
        "frozen_parameters": {
            "threshold_mgL": 5.0,
            "state_lower_bound_mgL": 0.0,
            "table2_L_mgL_per_hour": 13.203959999999936,
            "table2_max_replacements": 0,
            "table2_cadence_hours": [0.5, 1.0, 2.0, 4.0],
            "figure2_parameters": {
                "L_grid_mgL_per_hour": list(L_TICKS),
                "replacement_budget_rule": (
                    "station-specific minimum discrete replacement budget "
                    "k_min(L) on the frozen 24-h block"
                ),
                "window_rule": (
                    "earliest complete 24-h 15-min block at each frozen "
                    "validation station"
                ),
                "validation_station_ids": list(STATIONS),
            },
        },
    }
    manifest_path = OUTPUT / "FROZEN_ASSET_MANIFEST_20260727.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest_path


def main() -> None:
    figure = build_figure()
    tables = build_table()
    manifest = write_caption_and_manifest(figure, tables)
    print(figure)
    print(manifest)


if __name__ == "__main__":
    main()
