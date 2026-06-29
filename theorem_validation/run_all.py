#!/usr/bin/env python3
"""Run the four theorem-validation experiments and then the report figures."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def run(script: str) -> None:
    cmd = [sys.executable, os.path.join(THIS_DIR, script)]
    print("\n" + "=" * 72)
    print("RUN:", " ".join(cmd))
    print("=" * 72)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(f"{script} failed with exit code {result.returncode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all theorem-validation experiments.")
    parser.add_argument("--no-plots", action="store_true", help="Skip the plotting stage.")
    args = parser.parse_args()

    run("validate_boundedness.py")        # Experiment 1 -- boundedness/dissipativity (Thm 2.1)
    run("validate_local_stability.py")    # Experiment 2 -- local stability (Thm 2.2/2.4/2.8)
    run("validate_bifurcation.py")        # Experiments 3-4 -- persistence (Thm 3.4) + Hopf (Thm 4.3)
    if not args.no_plots:
        run("plots.py")

    print("\nAll theorem-validation experiments finished. "
          "Outputs are in theorem_validation/outputs/.")


if __name__ == "__main__":
    main()
