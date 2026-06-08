
from __future__ import annotations
import argparse
import csv
import os
from typing import Optional
import chemostat as cx
from chemostat import PARAMETER_NAMES, Params, State, params_from_dict

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA = os.path.join(THIS_DIR, "data", "baseline_scenarios.csv")
DEFAULT_OUTDIR = os.path.join(THIS_DIR, "outputs", "sensitivity")
DEFAULT_STATE0: State = (1.0, 0.5, 0.3, 0.1)

def read_baselines(path: str = DEFAULT_DATA) -> dict[str, dict[str, float]]:
    baselines: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = row["scenario"]
            baselines[name] = {key: float(row[key]) for key in PARAMETER_NAMES}
    return baselines


def with_override(base: dict[str, float], name: str, value: float) -> dict[str, float]:
    updated = dict(base)
    updated[name] = value
    return updated


def linspace(low: float, high: float, num: int) -> list[float]:
    """Inclusive evenly spaced values (pure-Python, avoids a numpy dependency)."""
    if num < 2:
        return [low]
    step = (high - low) / (num - 1)
    return [low + step * i for i in range(num)]


def representative_state(params: Params) -> Optional[State]:
    """Return the richest existing equilibrium's coordinates (theory side)."""
    equilibria = cx.compute_all_equilibria(params)
    for key in ("P3", "P2", "P1", "P0"):
        if equilibria[key] is not None:
            return equilibria[key]
    return None

# Experiment 1: parameter sweeps
DEFAULT_SWEEPS: dict[str, tuple[float, float, int]] = {
    "D1": (0.10, 4.50, 45),
    "D2": (0.10, 4.20, 42),
    "D3": (0.05, 1.30, 51),
}


def sweep_parameter(base: dict[str, float], param: str, low: float, high: float,
    num: int, state0: State, dt: float, t_end: float, method: str) -> list[dict[str, object]]:

    rows: list[dict[str, object]] = []
    for value in linspace(low, high, num):
        params = params_from_dict(with_override(base, param, value))

        equilibria = cx.compute_all_equilibria(params)
        exists = {name: (equilibria[name] is not None) for name in ("P0", "P1", "P2", "P3")}
        predicted = cx.highest_existing_equilibrium(params)
        theory_state = representative_state(params)

        sim = cx.integrate(state0, params, dt=dt, t_end=t_end, method=method)
        sim_regime = "exploded" if sim.exploded else cx.classify_regime(sim.final_state)

        ts = theory_state if theory_state is not None else (float("nan"),) * 4
        fs = sim.final_state
        rows.append(
            {
                "param": param,
                "value": value,
                "P0_exists": exists["P0"],
                "P1_exists": exists["P1"],
                "P2_exists": exists["P2"],
                "P3_exists": exists["P3"],
                "predicted_regime": predicted,
                "theory_S": ts[0], "theory_x": ts[1], "theory_y": ts[2], "theory_z": ts[3],
                "sim_regime": sim_regime,
                "sim_S": fs[0], "sim_x": fs[1], "sim_y": fs[2], "sim_z": fs[3],
                "agree": predicted == sim_regime,
            }
        )
    return rows


def write_sweep_csv(rows: list[dict[str, object]], path: str):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fields = [
        "param", "value",
        "P0_exists", "P1_exists", "P2_exists", "P3_exists",
        "predicted_regime", "theory_S", "theory_x", "theory_y", "theory_z",
        "sim_regime", "sim_S", "sim_x", "sim_y", "sim_z", "agree",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_sweep(rows: list[dict[str, object]], param: str):
    lines = [f"  {param} sweep ({len(rows)} points):"]
    previous = None
    for row in rows:
        regime = row["predicted_regime"]
        if regime != previous:
            lines.append(f"    {param} >= {row['value']:.4g}: regime {previous or '-'} -> {regime}")
            previous = regime
    # Report the z-washout threshold specifically (top predator survival).
    survive = [r["value"] for r in rows if r["P3_exists"]]
    if survive:
        lines.append(
            f"    coexistence (P3) for {param} in "
            f"[{min(survive):.4g}, {max(survive):.4g}]"
        )
    else:
        lines.append(f"    coexistence (P3) never reached over the swept {param} range")
    return "\n".join(lines)

# Experiment 2: local sensitivity (elasticity) coefficients
OUTPUT_LABELS = ("S3", "x3", "y3", "z3")

def coexistence_outputs(params: Params) -> Optional[State]:
    """Coordinates of the interior coexistence equilibrium P3, or None."""
    return cx.compute_P3(params)


def elasticity_row(base: dict[str, float], param: str, rel_delta: float, central: bool) -> dict[str, object]:

    """
        Compute S(output, param) for all four coexistence outputs.
        Forward scheme:   S = (q(p+) - q0)/q0 / rel_delta
        Central scheme:         S = (q(p+) - q(p-)) / (2 q0 rel_delta)
    """

    base_params = params_from_dict(base)
    q0 = coexistence_outputs(base_params)
    result: dict[str, object] = {"param": param, "baseline": base[param]}

    if q0 is None:
        for label in OUTPUT_LABELS:
            result[label] = float("nan")
        result["note"] = "no P3 at baseline"
        return result

    p_value = base[param]
    plus = params_from_dict(with_override(base, param, p_value * (1.0 + rel_delta)))
    q_plus = coexistence_outputs(plus)

    q_minus = None
    if central:
        minus = params_from_dict(with_override(base, param, p_value * (1.0 - rel_delta)))
        q_minus = coexistence_outputs(minus)

    for i, label in enumerate(OUTPUT_LABELS):
        if q_plus is None or (central and q_minus is None):
            result[label] = float("nan")
            continue
        if central:
            dq_rel = (q_plus[i] - q_minus[i]) / (2.0 * q0[i])  # type: ignore[index]
        else:
            dq_rel = (q_plus[i] - q0[i]) / q0[i]
        result[label] = dq_rel / rel_delta
    result["note"] = "" if q_plus is not None else "P3 lost under perturbation"
    return result


def elasticity_table(base: dict[str, float], rel_delta: float = 0.01, central: bool = True) -> list[dict[str, object]]:
    return [elasticity_row(base, name, rel_delta, central) for name in PARAMETER_NAMES]

def write_elasticity_csv(rows: list[dict[str, object]], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fields = ["param", "baseline", *OUTPUT_LABELS, "note"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_elasticity(rows: list[dict[str, object]]):
    """Rank parameters by overall influence on the coexistence state."""
    def magnitude(row: dict[str, object]) -> float:
        vals = [abs(row[label]) for label in OUTPUT_LABELS
                if isinstance(row[label], float) and row[label] == row[label]]
        return max(vals) if vals else float("nan")

    ranked = sorted(rows, key=lambda r: (magnitude(r) != magnitude(r), -magnitude(r)))
    lines = ["  Elasticity S(output, param) at the coexistence equilibrium:",
             "    param   |  S(S3)   S(x3)   S(y3)   S(z3)"]
    for row in rows:
        cells = "  ".join(
            f"{row[label]:+6.3f}" if isinstance(row[label], float) and row[label] == row[label]
            else "   nan"
            for label in OUTPUT_LABELS
        )
        lines.append(f"    {row['param']:<7} | {cells}")
    lines.append("    most influential (by max |S|): " +
                 ", ".join(str(r["param"]) for r in ranked[:3]))
    return "\n".join(lines)



def main():
    parser = argparse.ArgumentParser(description="Sensitivity analysis for the chemostat food chain.")
    parser.add_argument("--baseline", default="coexistence_P3",
                        help="Baseline scenario name from the data CSV.")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Baseline scenarios CSV.")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, help="Output directory.")
    parser.add_argument("--params", default="D1,D2,D3",
                        help="Comma-separated parameters to sweep (Experiment A).")
    parser.add_argument("--rel-delta", type=float, default=0.01,
                        help="Relative perturbation for elasticity (default 1%%).")
    parser.add_argument("--forward", action="store_true",
                        help="Use forward differences for elasticity (course style) "
                             "instead of the default central differences.")
    parser.add_argument("--dt", type=float, default=0.02,
                        help="Integration step for the sweep cross-check. Kept small so the "
                             "non-negativity clamping of the solver does not spuriously wash out "
                             "stiff transients at extreme parameter values.")
    parser.add_argument("--t-end", type=float, default=400.0, help="Integration horizon for sweeps.")
    parser.add_argument("--method", default="rk4", choices=["rk4", "euler"])
    args = parser.parse_args()

    baselines = read_baselines(args.data)
    if args.baseline not in baselines:
        raise SystemExit(f"Unknown baseline {args.baseline!r}; have {sorted(baselines)}")
    base = baselines[args.baseline]

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Baseline scenario: {args.baseline}")
    print(f"  parameters: {base}")


    # ---- Experiment 1 --- #
    print("\n=== Experiment 1 ===")
    all_rows: list[dict[str, object]] = []
    for param in [p.strip() for p in args.params.split(",") if p.strip()]:
        low, high, num = DEFAULT_SWEEPS.get(param, (0.1, 2.0, 40))
        rows = sweep_parameter(base, param, low, high, num,
                               DEFAULT_STATE0, args.dt, args.t_end, args.method)
        path = os.path.join(args.outdir, f"sweep_{param}.csv")
        write_sweep_csv(rows, path)
        all_rows.extend(rows)
        print(summarize_sweep(rows, param))
        print(f"    -> wrote {path}")
    write_sweep_csv(all_rows, os.path.join(args.outdir, "sweep_all.csv"))

    # ---- Experiment 2 --- #
    print("\n=== Experiment 2 ===")
    central = not args.forward
    rows = elasticity_table(base, rel_delta=args.rel_delta, central=central)
    path = os.path.join(args.outdir, "elasticity_coefficients.csv")
    write_elasticity_csv(rows, path)
    print(f"  scheme: {'central' if central else 'forward'} differences, "
          f"rel_delta = {args.rel_delta:.3g}")
    print(summarize_elasticity(rows))
    print(f"    -> wrote {path}")


if __name__ == "__main__":
    main()
