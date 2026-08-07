# -*- coding: utf-8 -*-
"""Build publication-oriented conceptual Figures 1 and 3."""

from __future__ import annotations

import hashlib
import html
import json
import math
from pathlib import Path

import cairosvg


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OUTPUT = PROJECT / "09_论文写作输出"
FIGURES = OUTPUT / "figures"

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
PINK = "#CC79A7"
RED = "#D55E00"
GRAY = "#666666"
LIGHT_GRAY = "#D9D9D9"
INK = "#202020"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 19,
    anchor: str = "middle",
    weight: str = "normal",
    fill: str = INK,
    rotate: float | None = None,
    italic: bool = False,
) -> str:
    transform = (
        f' transform="rotate({rotate:g} {x:.2f} {y:.2f})"'
        if rotate is not None
        else ""
    )
    style = ' font-style="italic"' if italic else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="Arial,Helvetica,sans-serif" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}"{style}{transform}>'
        f"{html.escape(value)}</text>"
    )


def _line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    stroke: str = INK,
    width: float = 2.0,
    dash: str | None = None,
    opacity: float = 1.0,
) -> str:
    dash_part = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{stroke}" stroke-width="{width:.2f}" '
        f'opacity="{opacity:.3f}"{dash_part}/>'
    )


def _path(
    points: list[tuple[float, float]],
    *,
    stroke: str = INK,
    width: float = 2.5,
    fill: str = "none",
    opacity: float = 1.0,
    dash: str | None = None,
) -> str:
    commands = " ".join(
        ("M" if index == 0 else "L") + f" {x:.2f} {y:.2f}"
        for index, (x, y) in enumerate(points)
    )
    dash_part = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<path d="{commands}" stroke="{stroke}" stroke-width="{width:.2f}" '
        f'fill="{fill}" opacity="{opacity:.3f}"{dash_part}/>'
    )


def _circle(
    x: float,
    y: float,
    *,
    radius: float = 7.0,
    fill: str = BLUE,
    stroke: str = "#ffffff",
    width: float = 1.5,
) -> str:
    return (
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width:.2f}"/>'
    )


def _cross(x: float, y: float, *, color: str = RED, size: float = 8.0) -> str:
    return "\n".join(
        [
            _line(x - size, y - size, x + size, y + size, stroke=color, width=3),
            _line(x - size, y + size, x + size, y - size, stroke=color, width=3),
        ]
    )


def _arrow(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    color: str = GRAY,
    width: float = 2.0,
) -> str:
    angle = math.atan2(y2 - y1, x2 - x1)
    head = 9.0
    wing = 0.55
    points = [
        (x2, y2),
        (
            x2 - head * math.cos(angle - wing),
            y2 - head * math.sin(angle - wing),
        ),
        (
            x2 - head * math.cos(angle + wing),
            y2 - head * math.sin(angle + wing),
        ),
    ]
    polygon = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return "\n".join(
        [
            _line(x1, y1, x2, y2, stroke=color, width=width),
            f'<polygon points="{polygon}" fill="{color}"/>',
        ]
    )


def _axes(
    left: float,
    top: float,
    width: float,
    height: float,
    *,
    x_label: str,
    y_label: str,
) -> list[str]:
    return [
        f'<rect x="{left:.2f}" y="{top:.2f}" width="{width:.2f}" '
        f'height="{height:.2f}" fill="#ffffff" stroke="#444444" '
        'stroke-width="1.5"/>',
        _text(left + width / 2, top + height + 42, x_label, size=17),
        _text(
            left - 46,
            top + height / 2,
            y_label,
            size=17,
            rotate=-90,
        ),
    ]


def _map_x(t: float, left: float, width: float, t_min: float, t_max: float) -> float:
    return left + width * (t - t_min) / (t_max - t_min)


def _map_y(v: float, top: float, height: float, v_min: float, v_max: float) -> float:
    return top + height * (v_max - v) / (v_max - v_min)


def _write_all(svg_path: Path, svg: str) -> list[Path]:
    svg_path.write_text(svg, encoding="utf-8")
    png_path = svg_path.with_suffix(".png")
    pdf_path = svg_path.with_suffix(".pdf")
    root = svg.split('viewBox="0 0 ', 1)[1].split('"', 1)[0].split()
    width, height = int(root[0]), int(root[1])
    cairosvg.svg2png(
        bytestring=svg.encode("utf-8"),
        write_to=str(png_path),
        output_width=2 * width,
        output_height=2 * height,
    )
    cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=str(pdf_path))
    return [svg_path, png_path, pdf_path]


def build_figure_1() -> list[Path]:
    width, height = 1800, 660
    panel_top, panel_height, panel_width = 185, 330, 460
    panel_lefts = (80, 670, 1260)
    times = [float(value) for value in range(7)]
    reports = [6.0, 5.2, 4.5, 7.2, 4.2, 5.0, 5.8]
    retained = [0, 1, 2, 4, 5, 6]
    lipschitz = 1.0
    threshold = 5.0
    value_min, value_max = 3.2, 7.8

    grid_times = [6.0 * index / 300 for index in range(301)]
    lower = [
        max(reports[index] - lipschitz * abs(t - times[index]) for index in retained)
        for t in grid_times
    ]
    upper = [
        min(reports[index] + lipschitz * abs(t - times[index]) for index in retained)
        for t in grid_times
    ]

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        _text(
            width / 2,
            42,
            "From uncertain records to four sharp burden endpoints",
            size=27,
            weight="bold",
        ),
    ]

    for index, (left, letter, title) in enumerate(
        zip(
            panel_lefts,
            ("A", "B", "C"),
            (
                "Published reports and an unknown corrupted subset",
                "Feasible Lipschitz band for one retained set",
                "Endpoint-specific extremal paths",
            ),
        )
    ):
        elements.append(_text(left - 36, 120, letter, size=30, weight="bold"))
        elements.append(_text(left + panel_width / 2, 120, title, size=18))
        if index < 2:
            elements.extend(
                _axes(
                    left,
                    panel_top,
                    panel_width,
                    panel_height,
                    x_label="Time",
                    y_label="Dissolved oxygen",
                )
            )

    left = panel_lefts[0]
    threshold_y = _map_y(
        threshold, panel_top, panel_height, value_min, value_max
    )
    elements.append(
        _line(
            left,
            threshold_y,
            left + panel_width,
            threshold_y,
            stroke=ORANGE,
            width=2,
            dash="7,5",
        )
    )
    elements.append(
        _text(
            left + panel_width - 6,
            threshold_y - 9,
            "threshold H",
            size=15,
            anchor="end",
            fill="#A86200",
        )
    )
    report_points = [
        (
            _map_x(t, left, panel_width, 0, 6),
            _map_y(v, panel_top, panel_height, value_min, value_max),
        )
        for t, v in zip(times, reports)
    ]
    elements.append(
        _path(report_points, stroke="#A0A0A0", width=2, dash="5,5")
    )
    for index, (x, y) in enumerate(report_points):
        if index == 3:
            elements.append(_cross(x, y))
            elements.append(
                _text(
                    x + 12,
                    y - 17,
                    "candidate corruption",
                    size=14,
                    anchor="start",
                    fill=RED,
                )
            )
        else:
            elements.append(_circle(x, y))
    left = panel_lefts[1]
    threshold_y = _map_y(
        threshold, panel_top, panel_height, value_min, value_max
    )
    lower_points = [
        (
            _map_x(t, left, panel_width, 0, 6),
            _map_y(v, panel_top, panel_height, value_min, value_max),
        )
        for t, v in zip(grid_times, lower)
    ]
    upper_points = [
        (
            _map_x(t, left, panel_width, 0, 6),
            _map_y(v, panel_top, panel_height, value_min, value_max),
        )
        for t, v in zip(grid_times, upper)
    ]
    band_points = upper_points + list(reversed(lower_points))
    band = " ".join(f"{x:.2f},{y:.2f}" for x, y in band_points)
    elements.append(
        f'<polygon points="{band}" fill="{BLUE}" opacity="0.14" stroke="none"/>'
    )
    elements.append(_path(upper_points, stroke=GREEN, width=3))
    elements.append(_path(lower_points, stroke=PINK, width=3))
    elements.append(
        _line(
            left,
            threshold_y,
            left + panel_width,
            threshold_y,
            stroke=ORANGE,
            width=2,
            dash="7,5",
        )
    )
    for index in retained:
        x = _map_x(times[index], left, panel_width, 0, 6)
        y = _map_y(reports[index], panel_top, panel_height, value_min, value_max)
        elements.append(_circle(x, y, radius=5.5))
    elements.append(
        _text(
            left + 12,
            panel_top + 27,
            "upper envelope u_S",
            size=15,
            anchor="start",
            fill=GREEN,
        )
    )
    elements.append(
        _text(
            left + 12,
            panel_top + 50,
            "lower envelope ℓ_S",
            size=15,
            anchor="start",
            fill=PINK,
        )
    )
    elements.append(
        _text(
            left + panel_width - 8,
            panel_top + panel_height - 13,
            "all feasible paths lie in the band",
            size=14,
            anchor="end",
            fill=GRAY,
        )
    )

    left = panel_lefts[2]
    mini_height = 130
    mini_gap = 54
    mini_tops = (panel_top, panel_top + mini_height + mini_gap)
    for mini_top in mini_tops:
        elements.append(
            f'<rect x="{left:.2f}" y="{mini_top:.2f}" width="{panel_width:.2f}" '
            f'height="{mini_height:.2f}" fill="#ffffff" stroke="#555555" '
            'stroke-width="1.2"/>'
        )
        h_y = _map_y(threshold, mini_top, mini_height, value_min, value_max)
        elements.append(
            _line(
                left,
                h_y,
                left + panel_width,
                h_y,
                stroke=ORANGE,
                width=1.8,
                dash="7,5",
            )
        )
    upper_mini = [
        (
            _map_x(t, left, panel_width, 0, 6),
            _map_y(v, mini_tops[0], mini_height, value_min, value_max),
        )
        for t, v in zip(grid_times, upper)
    ]
    lower_mini = [
        (
            _map_x(t, left, panel_width, 0, 6),
            _map_y(v, mini_tops[1], mini_height, value_min, value_max),
        )
        for t, v in zip(grid_times, lower)
    ]
    elements.append(_path(upper_mini, stroke=GREEN, width=3))
    elements.append(_path(lower_mini, stroke=PINK, width=3))
    elements.append(
        _text(
            left + 12,
            mini_tops[0] + 25,
            "u_S gives both lower endpoints",
            size=16,
            anchor="start",
            fill=GREEN,
        )
    )
    elements.append(
        _text(
            left + 12,
            mini_tops[1] + 25,
            "max(ℓ_S, B) gives both upper endpoints",
            size=16,
            anchor="start",
            fill=PINK,
        )
    )
    elements.append(
        _text(
            left + panel_width / 2,
            panel_top + panel_height + 42,
            "Separate DPs optimize over every admissible retained set",
            size=17,
        )
    )

    elements.append(_arrow(555, 350, 645, 350, color=GRAY))
    elements.append(_arrow(1145, 350, 1235, 350, color=GRAY))
    elements.append(
        _text(
            width / 2,
            625,
            "Deleted indices are endpoint certificates, not verified bad sensor records.",
            size=17,
            fill=GRAY,
        )
    )
    elements.append("</svg>")
    return _write_all(
        FIGURES / "Figure_1_endpoint_identification_concept.svg",
        "\n".join(elements) + "\n",
    )


def build_figure_3() -> list[Path]:
    width, height = 1800, 660
    panel_lefts = (80, 660, 1240)
    panel_width = 480
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        _text(
            width / 2,
            42,
            "Common-anchor geometry and the report-uniform boundary",
            size=27,
            weight="bold",
        ),
    ]

    for left, letter, title in zip(
        panel_lefts,
        ("A", "B", "C"),
        (
            "Two explanations share anchors when M > 2k",
            "Deficit difference is controlled by missing-run geometry",
            "Without a state bound, M ≤ 2k admits a counterexample",
        ),
    ):
        elements.append(_text(left - 35, 120, letter, size=30, weight="bold"))
        elements.append(_text(left + panel_width / 2, 120, title, size=18))

    left = panel_lefts[0]
    row_ys = (250, 390)
    labels = ("retained by explanation x", "retained by explanation z")
    set_x = {0, 2, 4, 6}
    set_z = {1, 2, 5, 6}
    common = set_x & set_z
    for row_y, label, retained_set, color in zip(
        row_ys,
        labels,
        (set_x, set_z),
        (BLUE, PINK),
    ):
        elements.append(_line(left + 20, row_y, left + 450, row_y, stroke="#AAAAAA"))
        elements.append(
            _text(left + 20, row_y - 34, label, size=16, anchor="start", fill=color)
        )
        for index in range(7):
            x = left + 25 + index * 70
            if index in common:
                elements.append(
                    _circle(x, row_y, radius=8, fill=INK, stroke="#ffffff")
                )
            elif index in retained_set:
                elements.append(
                    _circle(x, row_y, radius=7, fill=color, stroke="#ffffff")
                )
            else:
                elements.append(
                    _circle(
                        x,
                        row_y,
                        radius=6,
                        fill="#ffffff",
                        stroke="#AAAAAA",
                        width=1.5,
                    )
                )
        elements.append(_text(left + 25, row_y + 37, "t0", size=13))
        elements.append(_text(left + 445, row_y + 37, "t6", size=13))
    for index in sorted(common):
        x = left + 25 + index * 70
        elements.append(_line(x, 220, x, 430, stroke=INK, width=1.5, dash="4,4"))
    elements.append(
        _text(
            left + 235,
            492,
            "black points = common clean anchors",
            size=16,
            fill=INK,
        )
    )
    elements.append(
        _text(
            left + 235,
            525,
            "|Sx ∩ Sz| ≥ M − 2k ≥ 1",
            size=19,
            weight="bold",
        )
    )

    left = panel_lefts[1]
    base_y = 430
    anchor_x1, anchor_x2 = left + 70, left + 385
    peak_x = (anchor_x1 + anchor_x2) / 2
    peak_y = 210
    elements.append(_line(left + 35, base_y, left + 445, base_y, stroke=INK))
    elements.append(_line(anchor_x1, 180, anchor_x1, base_y + 15, stroke=LIGHT_GRAY))
    elements.append(_line(anchor_x2, 180, anchor_x2, base_y + 15, stroke=LIGHT_GRAY))
    triangle = (
        f"{anchor_x1:.2f},{base_y:.2f} "
        f"{peak_x:.2f},{peak_y:.2f} "
        f"{anchor_x2:.2f},{base_y:.2f}"
    )
    elements.append(
        f'<polygon points="{triangle}" fill="{BLUE}" opacity="0.18" '
        f'stroke="{BLUE}" stroke-width="2.5"/>'
    )
    elements.append(_circle(anchor_x1, base_y, radius=7, fill=INK))
    elements.append(_circle(anchor_x2, base_y, radius=7, fill=INK))
    elements.append(
        _text(
            peak_x,
            peak_y - 18,
            "|x(t) − z(t)| ≤ 2L min(t−a, b−t)",
            size=16,
            fill=BLUE,
        )
    )
    elements.append(
        _text(
            peak_x,
            base_y + 43,
            "internal gap g = b − a",
            size=16,
        )
    )
    elements.append(
        _text(
            peak_x,
            505,
            "area contribution ≤ Lg²/2",
            size=19,
            weight="bold",
        )
    )
    elements.append(
        _text(
            peak_x,
            542,
            "boundary gaps are one-sided and cost up to La²",
            size=15,
            fill=GRAY,
        )
    )

    left = panel_lefts[2]
    plot_top, plot_height = 190, 290
    elements.extend(
        _axes(
            left + 25,
            plot_top,
            430,
            plot_height,
            x_label="Time",
            y_label="Path value",
        )
    )
    x0, x1 = left + 55, left + 425
    y_minus, y_plus, y_h = plot_top + 225, plot_top + 65, plot_top + 145
    elements.append(_line(x0, y_h, x1, y_h, stroke=ORANGE, width=2, dash="7,5"))
    elements.append(_line(x0, y_minus, x1, y_minus, stroke=BLUE, width=3))
    elements.append(_line(x0, y_plus, x1, y_plus, stroke=PINK, width=3))
    elements.append(_circle(x0 + 40, y_minus, radius=7, fill=BLUE))
    elements.append(_circle(x1 - 40, y_plus, radius=7, fill=PINK))
    elements.append(
        _text(x1 - 4, y_h - 9, "H = 0", size=14, anchor="end", fill="#A86200")
    )
    elements.append(
        _text(
            x1 - 6,
            y_minus - 10,
            "x(t) = −A; retain y0 = −A",
            size=15,
            anchor="end",
            fill=BLUE,
        )
    )
    elements.append(
        _text(
            x1 - 6,
            y_plus - 10,
            "z(t) = +A; retain y1 = +A",
            size=15,
            anchor="end",
            fill=PINK,
        )
    )
    elements.append(
        _text(
            left + panel_width / 2,
            548,
            "M = 2, k = 1, L = 0; deficit gap = AT",
            size=18,
            weight="bold",
        )
    )
    elements.append(
        _text(
            left + panel_width / 2,
            590,
            "Report-uniform failure only; x(t) ≥ B supplies a separate finite cap.",
            size=14,
            fill=GRAY,
        )
    )

    elements.append("</svg>")
    return _write_all(
        FIGURES / "Figure_3_common_anchor_geometry.svg",
        "\n".join(elements) + "\n",
    )


def write_caption_and_manifest(assets: list[Path]) -> tuple[Path, Path]:
    captions = OUTPUT / "CONCEPTUAL_FIGURE_CAPTIONS.md"
    captions.write_text(
        "\n".join(
            [
                "# Conceptual figure captions",
                "",
                "## Figure 1",
                "",
                "**From uncertain records to four sharp burden endpoints.** "
                "(A) An illustrative published DO record under the exact-inlier "
                "gross-replacement sensitivity model; the crossed report is a "
                "candidate corruption for this witness, not a verified error. "
                "(B) For one compatible retained set, every admissible "
                "Lipschitz path lies between the lower and upper envelopes. "
                "(C) Monotonicity makes the upper envelope attain both lower "
                "functional endpoints and the floor-adjusted lower envelope "
                "attain both upper endpoints. Four ordered dynamic programs "
                "optimize these costs over all admissible retained sets.",
                "",
                "## Figure 3",
                "",
                "**Common-anchor geometry behind the report-uniform deficit "
                "bound.** (A) When $M>2k$, any two candidate retained sets "
                "share at least one report. (B) Between adjacent common anchors, "
                "the integrated path difference is bounded by a triangular "
                "missing-run contribution; one-sided boundary runs are more "
                "costly and generate the regular-grid contamination coefficient. "
                "(C) Without an external state bound, $M\\leq2k$ permits "
                "candidate explanations with no common anchor, so no finite "
                "report-uniform bound independent of report magnitude exists. "
                "The two-point construction is a design counterexample, not a "
                "claim that every fixed record is unbounded; an explicit state "
                "floor supplies a separate finite cap.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    assets_with_caption = [*assets, captions]
    manifest = OUTPUT / "CONCEPTUAL_FIGURE_MANIFEST_20260727.json"
    payload = {
        "status": "conceptual_figures_built",
        "source_script": {
            str(HERE / "make_conceptual_paper_figures.py"): _sha256(
                HERE / "make_conceptual_paper_figures.py"
            )
        },
        "assets": {str(path): _sha256(path) for path in assets_with_caption},
        "scope": (
            "conceptual/theoretical figures only; no empirical values or "
            "parameter estimates are introduced"
        ),
    }
    manifest.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return captions, manifest


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure_1 = build_figure_1()
    figure_3 = build_figure_3()
    captions, manifest = write_caption_and_manifest([*figure_1, *figure_3])
    print(figure_1[0])
    print(figure_3[0])
    print(captions)
    print(manifest)


if __name__ == "__main__":
    main()
