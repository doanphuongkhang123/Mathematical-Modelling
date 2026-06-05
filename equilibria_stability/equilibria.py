#!/usr/bin/env python3
"""Command-line entry point for equilibrium and stability-threshold calculations."""

from __future__ import annotations

import argparse
from io_utils import print_summary, read_scenario_csv, write_equilibrium_csv


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute P0-P3 equilibria for Monod chemostat scenarios.")
    parser.add_argument("--input", default="equilibria_stability/data/scenario_parameters.csv", help="Input scenario CSV.")
    parser.add_argument("--output", default="equilibria_stability/outputs/equilibrium_examples.csv", help="Output equilibrium CSV.")
    args = parser.parse_args()

    examples = read_scenario_csv(args.input)
    print_summary(examples)
    write_equilibrium_csv(examples, args.output)
    print(f"\nWrote equilibrium table to {args.output}")


if __name__ == "__main__":
    main()
