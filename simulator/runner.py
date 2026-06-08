#!/usr/bin/env python3
"""CLI entry point for running chemostat simulations.

Usage:
    python3 -m simulator.runner --scenario coexistence_P3
    python3 -m simulator.runner --scenario prey_only_P1 --dt 0.005 --t-end 200
    python3 -m simulator.runner --scenario coexistence_P3 --euler
    python3 -m simulator.runner --scenario coexistence_P3 --init 0.5,0.2,0.1,0.05
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

# Make equilibria_stability importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "equilibria_stability"))

from model import PARAMETER_NAMES, params_from_dict
from analysis import compute_all_equilibria

from .rhs import rhs
from .solvers import euler, rk4
from .diagnostics import diagnose


def read_scenario(csv_path: str, scenario_name: str) -> dict[str, float] | None:
    """Read parameters for a named scenario from the CSV file."""
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["scenario"] == scenario_name:
                return {name: float(row[name]) for name in PARAMETER_NAMES}
    return None


def write_csv(path: str, times: list[float], states: list[list[float]]) -> None:
    """Write trajectory to CSV with format t,S,x,y,z."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["t", "S", "x", "y", "z"])
        for t, s in zip(times, states):
            writer.writerow([f"{t:.6f}"] + [f"{v:.6f}" for v in s])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run chemostat food-chain simulation."
    )
    parser.add_argument(
        "--scenario", required=True, help="Scenario name from scenario_parameters.csv"
    )
    parser.add_argument(
        "--input",
        default="equilibria_stability/data/scenario_parameters.csv",
        help="Path to scenario parameters CSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/main_scenarios",
        help="Output directory for CSV and diagnostics.",
    )
    parser.add_argument("--dt", type=float, default=0.01, help="Time step size.")
    parser.add_argument("--t-end", type=float, default=100.0, help="Final simulation time.")
    parser.add_argument(
        "--euler", action="store_true", help="Also run Euler solver (default: RK4 only)."
    )
    parser.add_argument(
        "--init",
        default=None,
        help="Initial conditions as comma-separated values: S0,x0,y0,z0",
    )
    args = parser.parse_args()

    # Resolve paths relative to simulator/ directory
    sim_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(sim_dir)  # Mathematical-Modelling/
    input_path = os.path.join(project_dir, args.input)
    output_dir = os.path.join(project_dir, args.output_dir)

    # Read scenario parameters
    raw_params = read_scenario(input_path, args.scenario)
    if raw_params is None:
        print(f"Error: scenario '{args.scenario}' not found in {input_path}")
        sys.exit(1)

    params = params_from_dict(raw_params)

    # Initial conditions
    if args.init:
        state0 = [float(v) for v in args.init.split(",")]
        if len(state0) != 4:
            print("Error: --init must have exactly 4 values: S0,x0,y0,z0")
            sys.exit(1)
    else:
        # Default initial conditions — Person 4 may override per scenario
        state0 = [0.5, 0.2, 0.1, 0.05]

    # Build RHS function
    def rhs_fn(state):
        return rhs(state, params)

    # Equilibrium info
    equilibria = compute_all_equilibria(params)
    print(f"Scenario: {args.scenario}")
    print(f"Parameters: {raw_params}")
    print(f"Initial state: {state0}")
    print(f"Known equilibria:")
    for name, eq in equilibria.items():
        if eq is not None:
            print(f"  {name} = ({', '.join(f'{v:.6f}' for v in eq)})")
        else:
            print(f"  {name} = does not exist")

    # Run RK4
    print(f"\nRunning RK4 (dt={args.dt}, t_end={args.t_end})...")
    t_rk4, s_rk4, info_rk4 = rk4(rhs_fn, state0, t_end=args.t_end, dt=args.dt)
    diag_rk4 = diagnose(t_rk4, s_rk4, info_rk4, params)

    rk4_csv = os.path.join(output_dir, f"{args.scenario}_rk4.csv")
    write_csv(rk4_csv, t_rk4, s_rk4)
    print(f"  Wrote {rk4_csv}")
    print(f"  Final state: {diag_rk4['final_state']}")
    print(f"  Converged: {diag_rk4['converged']}")
    if diag_rk4.get("equilibrium_proximity"):
        print(f"  Equilibrium distances:")
        for name, dist in diag_rk4["equilibrium_proximity"].items():
            if dist is not None:
                print(f"    {name}: {dist:.6e}")

    # Run Euler if requested
    if args.euler:
        print(f"\nRunning Euler (dt={args.dt}, t_end={args.t_end})...")
        t_euler, s_euler, info_euler = euler(rhs_fn, state0, t_end=args.t_end, dt=args.dt)
        diag_euler = diagnose(t_euler, s_euler, info_euler, params)

        euler_csv = os.path.join(output_dir, f"{args.scenario}_euler.csv")
        write_csv(euler_csv, t_euler, s_euler)
        print(f"  Wrote {euler_csv}")
        print(f"  Final state: {diag_euler['final_state']}")
        print(f"  Converged: {diag_euler['converged']}")

    # Write diagnostics JSON
    diag_path = os.path.join(output_dir, f"{args.scenario}_diagnostics.json")
    os.makedirs(os.path.dirname(diag_path) or ".", exist_ok=True)
    with open(diag_path, "w", encoding="utf-8") as f:
        json.dump({"rk4": diag_rk4}, f, indent=2)
    print(f"\nWrote diagnostics to {diag_path}")


if __name__ == "__main__":
    main()
