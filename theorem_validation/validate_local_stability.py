
from __future__ import annotations

import argparse
import os

import _common as C

EPS = 1e-7
EQ_LEVEL = {"P0": 1, "P1": 2, "P2": 3, "P3": 4}   # number of surviving coords (incl. S)


def block_and_invasions(eq: C.State, params: C.Params, k: int):
    """Return (internal kxk block eigenvalues, [invasion rates of missing species])."""
    J = C.jacobian(eq, params)
    block = [[J[i][j] for j in range(k)] for i in range(k)]
    block_eigs = C.poly_roots(C.charpoly(block)) if k > 0 else []
    invasions = [J[i][i] for i in range(k, 4)]
    return block_eigs, invasions


def analytic_verdict(eq: C.State, params: C.Params, name: str):
    """Paper-side prediction: internal block stable AND all invasion rates < 0."""
    k = EQ_LEVEL[name]
    block_eigs, invasions = block_and_invasions(eq, params, k)
    block_max_re = max((e.real for e in block_eigs), default=-1.0)
    block_stable = block_max_re < -EPS
    invasions_neg = all(r < -EPS for r in invasions)
    # The single threshold the paper highlights for this equilibrium (the rate
    # at which the next trophic level invades); None for the top equilibrium P3.
    decisive = invasions[0] if invasions else None
    stable = block_stable and invasions_neg
    if stable:
        label = "stable"
    elif block_stable and not invasions_neg:
        label = "saddle (invasion)"   # boundary point a higher species can invade
    else:
        label = "unstable (block)"    # the surviving block itself is unstable
    return {"label": label, "stable": stable, "block_max_re": block_max_re,
            "invasions": invasions, "decisive": decisive}


def numeric_verdict(eq: C.State, params: C.Params):
    eigs = C.eigenvalues(eq, params)
    kind, max_re = C.classify_eigs(eigs, EPS)
    return {"kind": kind, "max_re": max_re, "stable": max_re < -EPS}


def sim_verdict(eq: C.State, params: C.Params, dt: float, t_end: float, delta: float):
    """Seed every component with +delta and see whether the regime is recovered."""
    start = tuple(v + delta for v in eq)
    sim = C.integrate(start, params, dt=dt, t_end=t_end, method="rk4")
    final = sim.final_state
    eq_regime = C.classify_regime(eq)
    final_regime = "exploded" if sim.exploded else C.classify_regime(final)
    returns = (final_regime == eq_regime) and (C.distance(final, eq) < 1e-2)
    return {"returns": returns, "final_regime": final_regime,
            "eq_regime": eq_regime, "dist": C.distance(final, eq)}


def study(scenarios: dict[str, dict[str, float]], dt: float, t_end: float, delta: float):
    rows: list[dict] = []
    for sname, base in scenarios.items():
        params = C.params_from_dict(base)
        eqs = C.compute_all_equilibria(params)
        for name in ("P0", "P1", "P2", "P3"):
            eq = eqs[name]
            if eq is None:
                continue
            a = analytic_verdict(eq, params, name)
            n = numeric_verdict(eq, params)
            s = sim_verdict(eq, params, dt, t_end, delta)
            agree = (a["stable"] == n["stable"] == s["returns"])
            rows.append({
                "scenario": sname,
                "equilibrium": name,
                "coords": "(" + ", ".join(f"{v:.4f}" for v in eq) + ")",
                "decisive_threshold": ("" if a["decisive"] is None else f"{a['decisive']:+.4f}"),
                "analytic": a["label"],
                "numeric_kind": n["kind"],
                "max_Re_eig": f"{n['max_re']:+.4e}",
                "sim_outcome": ("returns" if s["returns"] else f"departs->{s['final_regime']}"),
                "all_agree": agree,
            })
    return rows


def summarize(rows: list[dict]) -> str:
    n = len(rows)
    agree = sum(1 for r in rows if r["all_agree"])
    stable = sum(1 for r in rows if r["analytic"] == "stable")
    saddle = sum(1 for r in rows if r["analytic"].startswith("saddle"))
    lines = [
        "  TV2 -- Local-stability theorems (2.2, 2.4/2.5, 2.8/2.9):",
        f"    equilibria examined: {n}   (stable: {stable}, saddle/unstable: {n - stable})",
        f"    analytic == numeric == simulation: {agree}/{n}",
        "    verdict: " + ("CONFIRMED -- all three methods agree at every equilibrium"
                           if agree == n else
                           f"CHECK -- {n - agree} equilibria disagree"),
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="TV2: validate local-stability theorems.")
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--t-end", type=float, default=300.0)
    parser.add_argument("--delta", type=float, default=1e-3,
                        help="Positive seeding perturbation for the simulation test.")
    args = parser.parse_args()

    C.ensure_dirs()
    scenarios = C.read_scenarios()
    print("=== TV2: local stability  (analytic vs eigenvalues vs simulation) ===")
    rows = study(scenarios, args.dt, args.t_end, args.delta)

    fields = ["scenario", "equilibrium", "coords", "decisive_threshold", "analytic",
              "numeric_kind", "max_Re_eig", "sim_outcome", "all_agree"]
    path = os.path.join(C.DATA_DIR, "tv2_local_stability.csv")
    C.write_csv(rows, fields, path)
    print(summarize(rows))
    print(f"    -> wrote {path}")


if __name__ == "__main__":
    main()
