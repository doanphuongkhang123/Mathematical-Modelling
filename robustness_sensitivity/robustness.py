from __future__ import annotations
import argparse
import csv
import math
import os
from typing import Optional

import chemostat as cx
from chemostat import PARAMETER_NAMES, Params, State, params_from_dict
from sensitivity import DEFAULT_DATA, DEFAULT_STATE0, read_baselines


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTDIR = os.path.join(THIS_DIR, "outputs", "robustness")


def distance(a: State, b: State) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(4)))


def reference_equilibrium(params: Params) -> Optional[State]:
    """Theoretical target = richest existing equilibrium."""
    equilibria = cx.compute_all_equilibria(params)
    for key in ("P3", "P2", "P1", "P0"):
        if equilibria[key] is not None:
            return equilibria[key]
    return None



# Experiment 3: step-size / method robustness
def step_size_study(base: dict[str, float], state0: State, dt_values: list[float],
            t_acc: float, t_long: float,
            methods: tuple[str, ...] = ("rk4", "euler")) -> list[dict[str, object]]:
    """Integrate at each (method, dt) and measure self-convergence + stability."""

    params = params_from_dict(base)
    equilibrium = reference_equilibrium(params)
    ordered = sorted(dt_values, reverse=True)  # coarse -> fine

    rows: list[dict[str, object]] = []
    for method in methods:
        # Transient final states (at t_acc) for every dt, for the halving test.
        finals_acc: dict[float, Optional[State]] = {}
        for dt in ordered:
            sim_acc = cx.integrate(state0, params, dt=dt, t_end=t_acc, method=method)
            finals_acc[dt] = None if sim_acc.exploded else sim_acc.final_state

        previous_diff: Optional[float] = None
        for i, dt in enumerate(ordered):
            # Self-convergence: change in final state when the step is halved.
            diff_to_finer = float("nan")
            if i + 1 < len(ordered):
                here, finer = finals_acc[dt], finals_acc[ordered[i + 1]]
                if here is not None and finer is not None:
                    diff_to_finer = distance(here, finer)

            order = float("nan")
            if (previous_diff is not None and previous_diff > 0
                    and diff_to_finer == diff_to_finer and diff_to_finer > 0):
                order = math.log(previous_diff / diff_to_finer) / math.log(2.0)
            if diff_to_finer == diff_to_finer:
                previous_diff = diff_to_finer

            sim_long = cx.integrate(state0, params, dt=dt, t_end=t_long, method=method)
            fl = sim_long.final_state
            regime = "exploded" if sim_long.exploded else cx.classify_regime(fl)
            err_eq = (float("nan") if sim_long.exploded or equilibrium is None
                      else distance(fl, equilibrium))
            clamps = sim_long.clamps
            in_region = (not sim_long.exploded) and cx.in_invariant_region(fl, params)
            physical_ok = (not sim_long.exploded) and clamps == 0 and in_region

            rows.append(
                {
                    "method": method,
                    "dt": dt,
                    "diff_to_half_step": diff_to_finer,
                    "conv_order_est": order,
                    "regime_long": regime,
                    "final_S": fl[0], "final_x": fl[1], "final_y": fl[2], "final_z": fl[3],
                    "dist_to_equilibrium": err_eq,
                    "exploded": sim_long.exploded,
                    "clamps": clamps,
                    "inside_region": in_region,
                    "physical_ok": physical_ok,
                }
            )
    return rows


def write_step_size_csv(rows: list[dict[str, object]], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fields = [
        "method", "dt", "diff_to_half_step", "conv_order_est", "regime_long",
        "final_S", "final_x", "final_y", "final_z",
        "dist_to_equilibrium", "exploded", "clamps", "inside_region", "physical_ok",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_step_size(rows: list[dict[str, object]]) -> str:
    lines = ["  Step-size robustness (self-convergence under step halving):",
             "    method  dt        ||u(dt)-u(dt/2)||   order   regime  clamps  physical?"]
    for row in rows:
        order = row["conv_order_est"]
        order_str = f"{order:5.2f}" if isinstance(order, float) and order == order else "   - "
        diff = row["diff_to_half_step"]
        diff_str = f"{diff:.3e}" if isinstance(diff, float) and diff == diff else "    -    "
        if row["physical_ok"]:
            physical = "ok"
        elif row["exploded"]:
            physical = "FAIL(exploded)"
        elif row["clamps"]:
            physical = f"FAIL(clamped x{row['clamps']})"
        else:
            physical = "FAIL(outside Gamma)"
        lines.append(
            f"    {row['method']:<6} {row['dt']:<8.4g} {diff_str:>12}        "
            f"{order_str}    {str(row['regime_long']):<6}  {str(row['clamps']):>5}   {physical}"
        )
    return "\n".join(lines)

# Experiment 4: initial-condition robustness
DEFAULT_INITIAL_CONDITIONS: list[State] = [
    (1.0, 0.5, 0.3, 0.1),     # the standard start
    (0.2, 0.2, 0.2, 0.2),     # small, balanced
    (0.8, 1.0, 0.5, 0.3),     # large predators
    (0.5, 0.05, 0.05, 0.05),  # tiny organisms (near boundary)
    (1.0, 0.1, 0.6, 0.2),     # asymmetric
    (0.05, 0.8, 0.28, 0.08),  # already near the coexistence equilibrium
]


def initial_condition_study(base: dict[str, float], initial_conditions: list[State],
                    dt: float, t_end: float, method: str) -> list[dict[str, object]]:
    
    """Integrate from each initial condition and classify the outcome."""
    params = params_from_dict(base)
    equilibrium = reference_equilibrium(params)

    rows: list[dict[str, object]] = []
    for state0 in initial_conditions:
        sim = cx.integrate(state0, params, dt=dt, t_end=t_end, method=method)
        fs = sim.final_state
        regime = "exploded" if sim.exploded else cx.classify_regime(fs)
        err_eq = distance(fs, equilibrium) if equilibrium is not None else float("nan")
        rows.append(
            {
                "S0": state0[0], "x0": state0[1], "y0": state0[2], "z0": state0[3],
                "final_S": fs[0], "final_x": fs[1], "final_y": fs[2], "final_z": fs[3],
                "regime": regime,
                "dist_to_equilibrium": err_eq,
            }
        )
    return rows


def write_initial_condition_csv(rows: list[dict[str, object]], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fields = [
        "S0", "x0", "y0", "z0",
        "final_S", "final_x", "final_y", "final_z",
        "regime", "dist_to_equilibrium",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_initial_conditions(rows: list[dict[str, object]]) -> str:
    regimes = {row["regime"] for row in rows}
    max_err = max((row["dist_to_equilibrium"] for row in rows
                   if isinstance(row["dist_to_equilibrium"], float)), default=float("nan"))
    verdict = ("ROBUST: every initial condition reached the same regime"
               if len(regimes) == 1
               else f"NOT ROBUST: outcomes diverged into {sorted(regimes)}")
    lines = [f"  Initial-condition robustness ({len(rows)} starts):",
             f"    distinct outcomes: {sorted(regimes)}",
             f"    max distance to equilibrium across starts: {max_err:.3e}",
             f"    verdict: {verdict}"]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Robustness analysis for the chemostat food chain.")
    parser.add_argument("--baseline", default="coexistence_P3",
                        help="Baseline scenario name from the data CSV.")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Baseline scenarios CSV.")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, help="Output directory.")
    parser.add_argument("--dt0", type=float, default=0.4,
                        help="Coarsest step; the study uses dt0, dt0/2, dt0/4, ...")
    parser.add_argument("--levels", type=int, default=6,
                        help="Number of step-halving levels for Experiment C.")
    parser.add_argument("--t-acc", type=float, default=2.0,
                        help="Short transient horizon for the convergence-order estimate.")
    parser.add_argument("--t-end", type=float, default=400.0,
                        help="Long horizon for the long-run regime and IC study.")
    parser.add_argument("--method", default="rk4", choices=["rk4", "euler"],
                        help="Solver for the initial-condition study.")
    args = parser.parse_args()

    baselines = read_baselines(args.data)
    if args.baseline not in baselines:
        raise SystemExit(f"Unknown baseline {args.baseline!r}; have {sorted(baselines)}")
    base = baselines[args.baseline]

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Baseline scenario: {args.baseline}")
    print(f"  parameters: {base}")

    # ---- Experiment 3 --- #
    print("\n=== Experiment 3: step-size / method robustness ===")
    dt_values = [args.dt0 / (2 ** i) for i in range(args.levels)]
    rows = step_size_study(base, DEFAULT_STATE0, dt_values, t_acc=args.t_acc, t_long=args.t_end)
    path = os.path.join(args.outdir, "step_size.csv")
    write_step_size_csv(rows, path)
    print(summarize_step_size(rows))
    print(f"    -> wrote {path}")

    # ---- Experiment 4 --- #
    print("\n=== Experiment 4: initial-condition robustness ===")
    rows = initial_condition_study(base, DEFAULT_INITIAL_CONDITIONS,
                                   dt=min(dt_values), t_end=args.t_end, method=args.method)
    path = os.path.join(args.outdir, "initial_conditions.csv")
    write_initial_condition_csv(rows, path)
    print(summarize_initial_conditions(rows))
    print(f"    -> wrote {path}")


if __name__ == "__main__":
    main()
