#!/usr/bin/env python3
"""Run quick verification or the complete frozen MDDNR+USGS analysis."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CURRENT = ROOT / "03_代码与测试" / "current"
SNAPSHOT = ROOT / "03_代码与测试" / "原始快照"
LOG_DIR = ROOT / "verification_logs"

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def run(command: list[str], env: dict[str, str], log) -> None:
    rendered = " ".join(command)
    print(f"\n$ {rendered}", flush=True)
    log.write(f"\n$ {rendered}\n")
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="")
        log.write(line)
        log.flush()
    code = process.wait()
    if code:
        raise SystemExit(code)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full",
        action="store_true",
        help="recompute MDDNR and USGS analyses plus all frozen assets",
    )
    args = parser.parse_args()

    LOG_DIR.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(CURRENT), str(SNAPSHOT)))
    env["PYTHONIOENCODING"] = "utf-8"
    commands = [
        [
            sys.executable,
            "-m",
            "unittest",
            "discover",
            "-s",
            str(SNAPSHOT),
            "-p",
            "test_*.py",
            "-v",
        ],
        [
            sys.executable,
            "-m",
            "unittest",
            "discover",
            "-s",
            str(CURRENT),
            "-p",
            "test_*.py",
            "-v",
        ],
    ]
    if args.full:
        commands.extend(
            [
                [sys.executable, str(CURRENT / "run_maracoos_nonnegative_table2.py")],
                [sys.executable, str(CURRENT / "run_maracoos_nonnegative_trilemma.py")],
                [sys.executable, str(CURRENT / "run_maracoos_fixed_k1_sensitivity.py")],
                [sys.executable, str(CURRENT / "make_frozen_paper_assets.py")],
                [sys.executable, str(CURRENT / "make_conceptual_paper_figures.py")],
                [sys.executable, str(CURRENT / "make_model_assumption_table.py")],
                [sys.executable, str(CURRENT / "make_k1_supplement_table.py")],
                [sys.executable, str(CURRENT / "generate_table_s3_block_descriptives.py")],
                [sys.executable, str(CURRENT / "run_usgs_external_replication.py")],
                [sys.executable, str(CURRENT / "verify_usgs_external_replication.py")],
                [sys.executable, str(CURRENT / "make_revision_o_assets.py")],
            ]
        )
    else:
        commands.append(
            [sys.executable, str(CURRENT / "verify_usgs_external_replication.py")]
        )
    verify = [sys.executable, str(ROOT / "verify_reproduction.py")]
    if not args.full:
        verify.append("--package-only")
    commands.append(verify)

    log_path = LOG_DIR / ("full_reproduction.log" if args.full else "verification.log")
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        for command in commands:
            run(command, env, log)
    print(f"\nPASS; log: {log_path}")


if __name__ == "__main__":
    main()
