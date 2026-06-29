
from __future__ import annotations

import argparse
import os

import _common as C


# A spread of starts: small, balanced, predator-heavy, and two extreme totals
# that begin well outside the absorbing band (one far above, one far below).
DEFAULT_STARTS: list[C.State] = [
    (1.0, 0.5, 0.3, 0.1),     # standard start
    (0.2, 0.2, 0.2, 0.2),     # small balanced
    (0.8, 1.0, 0.5, 0.3),     # predator-heavy
    (5.0, 5.0, 5.0, 5.0),     # very large total -> starts above the band
    (8.0, 0.0, 0.0, 0.0),     # large nutrient pulse, no organisms
    (0.02, 0.02, 0.02, 0.02), # tiny total -> starts below the band
]


def total_mass(state: C.State) -> float:
    return sum(state)


def boundedness_study(base: dict[str, float], starts: list[C.State],
                      dt: float, t_end: float, k: float) -> tuple[list[dict], dict]:
    params = C.params_from_dict(base)
    d_max = max(1.0, params.D1, params.D2, params.D3)
    d_min = min(1.0, params.D1, params.D2, params.D3)
    lower, upper = 1.0 / d_max, 1.0 / d_min          # absorbing interval (k -> 0)
    band_lo, band_hi = lower - k, upper + k          # Gamma_k of Theorem 2.1

    rows: list[dict] = []
    series: dict[str, list[tuple[float, float]]] = {}
    tol = 1e-9

    for start in starts:
        sim = C.integrate(start, params, dt=dt, t_end=t_end, method="rk4", sample_every=20)
        times, traj = sim.times, sim.trajectory
        masses = [total_mass(s) for s in traj]

        # Non-negativity over the whole sampled trajectory.  Theorem 2.1(i) is a
        # claim about the continuous flow, so the test is that no component goes
        # negative; the solver's non-negativity clamp count is reported
        # separately as a numerical-quality signal (a clamp of size ~1e-16 at the
        # boundary is consistent with a component staying >= 0).
        min_component = min(min(s) for s in traj)
        nonneg_ok = min_component >= -tol and not sim.exploded

        # First time T enters the absorbing interval and stays in Gamma_k afterwards.
        entry_time = None
        for t, m in zip(times, masses):
            if lower - tol <= m <= upper + tol:
                entry_time = t
                break
        stayed = all(band_lo - tol <= m <= band_hi + tol
                     for t, m in zip(times, masses)
                     if entry_time is not None and t >= entry_time)

        final_mass = masses[-1]
        rows.append({
            "S0": start[0], "x0": start[1], "y0": start[2], "z0": start[3],
            "T0": round(total_mass(start), 6),
            "min_component": f"{min_component:.3e}",
            "nonneg_ok": nonneg_ok,
            "solver_clamps": sim.clamps,
            "max_T": f"{max(masses):.6f}",
            "final_T": f"{final_mass:.6f}",
            "band_lower": f"{lower:.6f}", "band_upper": f"{upper:.6f}",
            "entered_band": entry_time is not None,
            "entry_time": "" if entry_time is None else f"{entry_time:.2f}",
            "stayed_inside": stayed,
            "final_in_band": bool(lower - tol <= final_mass <= upper + tol),
        })
        series[str(tuple(round(v, 2) for v in start))] = list(zip(times, masses))

    meta = {"lower": lower, "upper": upper, "band_lo": band_lo, "band_hi": band_hi,
            "d_max": d_max, "d_min": d_min, "series": series}
    return rows, meta


def summarize(rows: list[dict], meta: dict) -> str:
    n = len(rows)
    nonneg = sum(1 for r in rows if r["nonneg_ok"])
    entered = sum(1 for r in rows if r["entered_band"])
    stayed = sum(1 for r in rows if r["stayed_inside"])
    lines = [
        "  TV1 -- Boundedness & dissipativity (Theorem 2.1):",
        f"    absorbing interval for total mass T: "
        f"[1/D_max, 1/D_min] = [{meta['lower']:.4f}, {meta['upper']:.4f}]",
        f"    non-negativity preserved (no component < 0):       {nonneg}/{n} starts",
        f"    total mass entered the band:                       {entered}/{n} starts",
        f"    total mass stayed in Gamma_k afterwards:           {stayed}/{n} starts",
        "    verdict: " + ("CONFIRMED -- every start is attracted into one bounded set"
                           if nonneg == entered == stayed == n
                           else "CHECK -- at least one start did not behave as predicted"),
    ]
    return "\n".join(lines)


def write_series(meta: dict, path: str) -> None:
    """Write total-mass time series (one column per start) for plotting."""
    series = meta["series"]
    names = list(series)
    # All share the same sampled times.
    times = [t for t, _ in next(iter(series.values()))]
    fields = ["t"] + names
    rows = []
    for i, t in enumerate(times):
        row = {"t": f"{t:.3f}"}
        for name in names:
            row[name] = f"{series[name][i][1]:.6f}"
        rows.append(row)
    C.write_csv(rows, fields, path)


def main() -> None:
    parser = argparse.ArgumentParser(description="TV1: validate Theorem 2.1 (boundedness/dissipativity).")
    parser.add_argument("--scenario", default="coexistence_P3")
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--t-end", type=float, default=120.0)
    parser.add_argument("--k", type=float, default=0.5, help="Gamma_k slack constant.")
    args = parser.parse_args()

    C.ensure_dirs()
    scenarios = C.read_scenarios()
    if args.scenario not in scenarios:
        raise SystemExit(f"Unknown scenario {args.scenario!r}; have {sorted(scenarios)}")
    base = scenarios[args.scenario]

    print(f"=== TV1: Theorem 2.1 (boundedness / dissipativity) -- scenario {args.scenario} ===")
    rows, meta = boundedness_study(base, DEFAULT_STARTS, args.dt, args.t_end, args.k)

    fields = ["S0", "x0", "y0", "z0", "T0", "min_component", "nonneg_ok", "solver_clamps",
              "max_T", "final_T", "band_lower", "band_upper",
              "entered_band", "entry_time", "stayed_inside", "final_in_band"]
    path = os.path.join(C.DATA_DIR, "tv1_boundedness.csv")
    C.write_csv(rows, fields, path)
    write_series(meta, os.path.join(C.DATA_DIR, "tv1_total_mass_series.csv"))
    print(summarize(rows, meta))
    print(f"    -> wrote {path}")


if __name__ == "__main__":
    main()
