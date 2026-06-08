"""Post-simulation diagnostics for chemostat trajectories."""

from __future__ import annotations

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "equilibria_stability"))
from model import Params
from analysis import compute_all_equilibria


def diagnose(
    times: list[float],
    states: list[list[float]],
    solver_info: dict | None = None,
    params: Params | None = None,
) -> dict:
    """Produce a diagnostics dictionary from a completed simulation.

    Parameters
    ----------
    times : list[float]
        Time stamps from the solver.
    states : list[list[float]]
        State vectors from the solver.
    solver_info : dict, optional
        Info dict returned by the solver (clamps, exploded, steps).
    params : Params, optional
        If provided, compute distance to each known equilibrium.

    Returns
    -------
    dict
        Diagnostic summary.
    """
    if not states:
        return {"error": "empty trajectory"}

    final = states[-1]

    # Per-component min/max across trajectory
    n_comp = len(final)
    mins = [min(s[i] for s in states) for i in range(n_comp)]
    maxs = [max(s[i] for s in states) for i in range(n_comp)]

    # Convergence: std of last 10% of trajectory < 1e-6 for all components
    tail_start = max(1, len(states) - len(states) // 10)
    tail = states[tail_start:]
    converged = True
    if len(tail) < 2:
        converged = len(states) <= 2  # trivially converged if very short
    else:
        for i in range(n_comp):
            vals = [s[i] for s in tail]
            mean = sum(vals) / len(vals)
            variance = sum((v - mean) ** 2 for v in vals) / len(vals)
            if math.sqrt(variance) > 1e-6:
                converged = False
                break

    # Equilibrium proximity
    eq_distances = {}
    if params is not None:
        equilibria = compute_all_equilibria(params)
        for name, eq_state in equilibria.items():
            if eq_state is None:
                eq_distances[name] = None
            else:
                dist = math.sqrt(sum((final[i] - eq_state[i]) ** 2 for i in range(4)))
                eq_distances[name] = dist

    result = {
        "final_state": final,
        "min_values": mins,
        "max_values": maxs,
        "negative_clamps": solver_info.get("total_clamps", 0) if solver_info else 0,
        "exploded": solver_info.get("exploded", False) if solver_info else False,
        "converged": converged,
        "total_time": times[-1] if times else 0.0,
        "total_steps": len(states) - 1,
    }
    if eq_distances:
        result["equilibrium_proximity"] = eq_distances

    return result
