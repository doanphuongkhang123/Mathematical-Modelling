from __future__ import annotations
import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_THIS_DIR)
EQUILIBRIA_STABILITY_DIR = os.path.join(_PROJECT_ROOT, "equilibria_stability")
for path in (_PROJECT_ROOT, EQUILIBRIA_STABILITY_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)

from model import PARAMETER_NAMES, Monod, Params, State, params_from_dict
from analysis import compute_all_equilibria, compute_P0, compute_P1, compute_P2, compute_P3, p0_stability

from simulator.rhs import rhs
from simulator.solvers import euler, rk4

__all__ = [
    # re-exported model / equilibrium API
    "PARAMETER_NAMES", "Monod", "Params", "State", "params_from_dict",
    "compute_all_equilibria", "compute_P0", "compute_P1", "compute_P2",
    "compute_P3", "p0_stability",
    # re-exported simulator API
    "rhs", "euler", "rk4",
    # helpers added here
    "integrate", "SimulationResult", "classify_regime",
    "highest_existing_equilibrium", "in_invariant_region", "eps_nonneg",
]
eps_nonneg = 1e-6
_SOLVERS = {"euler": euler, "rk4": rk4}

class SimulationResult:
    def __init__(self, method: str, dt: float, t_end: float, state0: State, final_state: State,
        times: list[float], trajectory: list[State], exploded: bool, clamps: int, 
        min_state: State, max_state: State):

        self.method = method
        self.dt = dt
        self.t_end = t_end
        self.state0 = state0
        self.final_state = final_state
        self.times = times           
        self.trajectory = trajectory  
        self.exploded = exploded      
        self.clamps = clamps         
        self.min_state = min_state    
        self.max_state = max_state    


def integrate(state0: State, params: Params, dt: float = 0.05, t_end: float = 600.0,
        method: str = "rk4",sample_every: int = 0) -> SimulationResult:

    if method not in _SOLVERS:
        raise ValueError(f"Unknown method {method!r}; choose from {sorted(_SOLVERS)}")

    rhs_fn = lambda s: rhs(s, params)  # noqa: E731 - simple adapter closure
    times, states, info = _SOLVERS[method](rhs_fn, list(state0), t_end=t_end, dt=dt)

    final_state: State = tuple(states[-1])  # type: ignore[assignment]
    min_state: State = tuple(min(s[i] for s in states) for i in range(4))  # type: ignore[assignment]
    max_state: State = tuple(max(s[i] for s in states) for i in range(4))  # type: ignore[assignment]

    traj_times: list[float] = []
    trajectory: list[State] = []
    if sample_every > 0:
        indices = list(range(0, len(states), sample_every))
        if indices[-1] != len(states) - 1:
            indices.append(len(states) - 1)
        traj_times = [times[i] for i in indices]
        trajectory = [tuple(states[i]) for i in indices]  # type: ignore[misc]

    return SimulationResult(
        method=method,
        dt=dt,
        t_end=t_end,
        state0=tuple(state0),  # type: ignore[arg-type]
        final_state=final_state,
        times=traj_times,
        trajectory=trajectory,
        exploded=bool(info.get("exploded", False)),
        clamps=int(info.get("total_clamps", 0)),
        min_state=min_state,
        max_state=max_state,
    )

def classify_regime(state: State, eps: float = 1e-3) -> str:
    _S, x, y, z = state
    if x <= eps:
        return "P0"
    if y <= eps:
        return "P1"
    if z <= eps:
        return "P2"
    return "P3"


def highest_existing_equilibrium(params: Params) -> str:
    """Theory-side prediction: the richest equilibrium that exists for these params."""
    equilibria = compute_all_equilibria(params)
    for name in ("P3", "P2", "P1", "P0"):
        if equilibria[name] is not None:
            return name
    return "P0"


def in_invariant_region(state: State, params: Params, k: float = 0.5) -> bool:
    """
        Check in bounded region Gamma of Theorem 2.1.
        Gamma = { (S,x,y,z) >= 0 : 1/Dmax - k <= S+x+y+z <= 1/Dmin + k },
    """
    if any(v < -eps_nonneg for v in state):
        return False
    total = sum(state)
    d_max = max(1.0, params.D1, params.D2, params.D3)
    d_min = min(1.0, params.D1, params.D2, params.D3)
    lower = 1.0 / d_max - k
    upper = 1.0 / d_min + k
    return lower <= total <= upper
