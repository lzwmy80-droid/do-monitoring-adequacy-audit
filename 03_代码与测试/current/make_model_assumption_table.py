# -*- coding: utf-8 -*-
"""Build the manuscript table that separates model inputs from interpretations."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OUTPUT = PROJECT / "09_论文写作输出"
TABLES = OUTPUT / "tables"

ROWS = [
    {
        "item": r"Reports \((t_i,y_i)\)",
        "model_role": "Discrete published observations at strictly ordered times.",
        "do_stress_test": "Public MARACOOS/Maryland DNR DO records in mg/L; missing or nonnumeric values are not imputed.",
        "not_implied": "The published value is not assumed to be error-free unless its index is retained by a candidate explanation.",
    },
    {
        "item": r"Latent path \(x(t)\)",
        "model_role": "A continuous real-valued path over the observed horizon.",
        "do_stress_test": "Represents a candidate continuous DO history.",
        "not_implied": "It is not a fitted ecological process model or a probabilistic prediction.",
    },
    {
        "item": r"Lipschitz bound \(L\)",
        "model_role": r"Restricts \(\lvert x(t)-x(s)\rvert\le L\lvert t-s\rvert\).",
        "do_stress_test": "Frozen calibration-slope labels and a seven-value diagnostic grid.",
        "not_implied": "An empirical adjacent slope is not a validated physical bound on every continuous path.",
    },
    {
        "item": r"Replacement budget \(k\)",
        "model_role": r"Permits at most \(k\) arbitrarily corrupted report values.",
        "do_stress_test": r"\(k_{\min}(L)\) is used only for the 24-h feasibility diagnostic; \(k=0\) is the sampling baseline and fixed \(k=1\) is a paired sensitivity.",
        "not_implied": r"Neither \(k\) nor \(k_{\min}\) is a verified error count; fixed \(k=1\) is not a constant contamination rate.",
    },
    {
        "item": r"Retained set \(S\)",
        "model_role": r"A candidate-specific subset with \(\lvert S\rvert\ge M-k\) and \(x(t_i)=y_i\) for \(i\in S\).",
        "do_stress_test": "Optimized separately for each endpoint.",
        "not_implied": "Indices outside one endpoint witness are not thereby diagnosed as bad sensor records.",
    },
    {
        "item": r"Inlier tolerance \(\epsilon=0\)",
        "model_role": "Retained reports are exact path anchors.",
        "do_stress_test": "The current endpoint algorithm uses exact inliers.",
        "not_implied": "Ordinary measurement noise, rounding uncertainty, and a joint tolerance-plus-replacement model are not covered.",
    },
    {
        "item": r"Threshold \(H\)",
        "model_role": "Defines strict occupation and hinge deficit.",
        "do_stress_test": r"\(H=5\) mg/L is a policy-relevant demonstration threshold.",
        "not_implied": "It is not asserted to be the universal ecological or regulatory criterion for every station, season, or designated use.",
    },
    {
        "item": r"State floor \(B\)",
        "model_role": r"Optionally restricts \(x(t)\ge B\).",
        "do_stress_test": r"\(B=0\) mg/L enforces physical DO nonnegativity.",
        "not_implied": "It is a path-state constraint, not a deterministic bound on report error.",
    },
    {
        "item": r"\(O_H,D_H\)",
        "model_role": r"Strict-threshold occupation \(1\{x(t)<H\}\) and cumulative deficit \((H-x(t))_+\).",
        "do_stress_test": r"Widths are normalized by \(T\) and \(HT\).",
        "not_implied": r"The four extrema need not imply that every value between them is attainable; results are endpoint hulls unless connectedness is proved.",
    },
    {
        "item": "15-min operational reference",
        "model_role": "A retained high-frequency series used for deterministic checks.",
        "do_stress_test": r"Containment is checked only when the reference is compatible with the frozen \(L\).",
        "not_implied": "Containment is not statistical coverage, external validation, or proof of the true continuous burden.",
    },
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tex_escape_plain(text: str) -> str:
    return (
        text.replace("&", r"\&")
        .replace("%", r"\%")
        .replace("#", r"\#")
    )


def build() -> Path:
    TABLES.mkdir(parents=True, exist_ok=True)
    csv_path = TABLES / "Table_1_model_assumptions.csv"
    md_path = TABLES / "Table_1_model_assumptions.md"
    tex_path = TABLES / "Table_1_model_assumptions.tex"
    caption_path = OUTPUT / "TABLE_1_CAPTION.md"

    fields = ("item", "model_role", "do_stress_test", "not_implied")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(ROWS)

    md_lines = [
        "| Item | Model role | DO stress-test setting | Does not imply |",
        "|---|---|---|---|",
    ]
    for row in ROWS:
        md_lines.append(
            f"| {row['item']} | {row['model_role']} | "
            f"{row['do_stress_test']} | {row['not_implied']} |"
        )
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    tex_lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Model assumptions, stress-test settings, and interpretation limits.}",
        r"\label{tab:model_assumptions}",
        r"\begin{tabular}{p{0.13\textwidth}p{0.25\textwidth}p{0.27\textwidth}p{0.27\textwidth}}",
        r"\toprule",
        r"Item & Model role & DO stress-test setting & Does not imply \\",
        r"\midrule",
    ]
    for row in ROWS:
        tex_lines.append(
            " & ".join(_tex_escape_plain(row[field]) for field in fields)
            + r" \\"
        )
    tex_lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table*}"])
    tex_path.write_text("\n".join(tex_lines) + "\n", encoding="utf-8")

    caption_path.write_text(
        "\n".join(
            [
                "# Table 1 caption",
                "",
                "**Model assumptions, stress-test settings, and "
                "interpretation limits.** The general theory permits an "
                r"optional state floor; the DO application sets \(B=0\) "
                "mg/L. Replacement sets and their retained/deleted "
                "certificates are endpoint-specific. The exact-inlier model "
                "does not include ordinary measurement tolerance, and the "
                "15-min operational reference supplies deterministic "
                "consistency checks rather than statistical validation.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    assets = (csv_path, md_path, tex_path, caption_path)
    script = Path(__file__).resolve()
    manifest = OUTPUT / "TABLE_1_ASSET_MANIFEST_20260727.json"
    manifest.write_text(
        json.dumps(
            {
                "status": "model_assumption_table_frozen",
                "source_script": {str(script): _sha256(script)},
                "assets": {str(path): _sha256(path) for path in assets},
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


if __name__ == "__main__":
    print(build())
