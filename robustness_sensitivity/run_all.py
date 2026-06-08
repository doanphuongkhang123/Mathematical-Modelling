#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def run(script: str, extra: list[str]) -> None:
    cmd = [sys.executable, os.path.join(THIS_DIR, script), *extra]
    print("\n" + "=" * 72)
    print("RUN:", " ".join(cmd))
    print("=" * 72)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(f"{script} failed with exit code {result.returncode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all robustness/sensitivity experiments.")
    parser.add_argument("--baseline", default="coexistence_P3",
                        help="Baseline scenario name (must exist in data/baseline_scenarios.csv).")
    parser.add_argument("--no-plots", action="store_true", help="Skip the plotting stage.")
    args = parser.parse_args()

    run("sensitivity.py", ["--baseline", args.baseline])
    run("robustness.py", ["--baseline", args.baseline])
    if not args.no_plots:
        run("plots.py", ["--baseline", args.baseline])

    print("\nAll experiments finished. Outputs are in robustness_sensitivity/outputs/.")


if __name__ == "__main__":
    main()
