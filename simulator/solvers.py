"""Numerical ODE solvers: Forward Euler and classic RK4.

Both solvers clamp state variables to >= 0 after each step (enforcing
biological non-negativity as guaranteed by Theorem 2.1) and stop early
if any component explodes above 1e10.
"""

from __future__ import annotations

import math
from typing import Callable


def _validate_inputs(state0: list[float], t_end: float, dt: float) -> None:
    """Validate common solver inputs before integration starts."""
    if not state0:
        raise ValueError("state0 must contain at least one component")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError(f"dt must be finite and positive, got {dt}")
    if not math.isfinite(t_end) or t_end < 0.0:
        raise ValueError(f"t_end must be finite and non-negative, got {t_end}")
    if any(not math.isfinite(value) for value in state0):
        raise ValueError("state0 values must be finite")


def _clamp(state: list[float]) -> tuple[list[float], int]:
    """Clamp negative components to 0. Returns (clamped_state, clamp_count)."""
    clamped = 0
    out = []
    for v in state:
        if v < 0.0:
            out.append(0.0)
            clamped += 1
        else:
            out.append(v)
    return out, clamped


def _check_explode(state: list[float], limit: float = 1e10) -> bool:
    """Return True if any component exceeds the explosion limit."""
    return any(abs(v) > limit for v in state)


def euler(
    rhs_fn: Callable[[list[float]], list[float]],
    state0: list[float],
    t_end: float = 100.0,
    dt: float = 0.01,
) -> tuple[list[float], list[list[float]], dict]:
    """Forward Euler solver.

    Parameters
    ----------
    rhs_fn : callable
        Function taking a state vector and returning derivatives.
        Typically ``lambda s: rhs(s, params)``.
    state0 : list[float]
        Initial state [S, x, y, z].
    t_end : float
        Final simulation time.
    dt : float
        Time step size.

    Returns
    -------
    times : list[float]
        Time stamps.
    states : list[list[float]]
        State vectors at each time stamp.
    info : dict
        Diagnostics: total_clamps, exploded, steps.
    """
    _validate_inputs(state0, t_end, dt)
    n = len(state0)
    n_steps = int(math.ceil(t_end / dt))
    times = [0.0]
    states = [list(state0)]
    total_clamps = 0
    exploded = False

    for i in range(n_steps):
        t = times[-1]
        y = states[-1]
        step = min(dt, t_end - t)
        dydt = rhs_fn(y)

        new_y = [y[j] + step * dydt[j] for j in range(n)]

        new_y, clamps = _clamp(new_y)
        total_clamps += clamps

        if _check_explode(new_y):
            exploded = True
            times.append(t + step)
            states.append(new_y)
            break

        times.append(t + step)
        states.append(new_y)

    return times, states, {
        "total_clamps": total_clamps,
        "exploded": exploded,
        "steps": len(times) - 1,
    }


def rk4(
    rhs_fn: Callable[[list[float]], list[float]],
    state0: list[float],
    t_end: float = 100.0,
    dt: float = 0.01,
) -> tuple[list[float], list[list[float]], dict]:
    """Classic 4th-order Runge-Kutta solver.

    Parameters
    ----------
    rhs_fn : callable
        Function taking a state vector and returning derivatives.
    state0 : list[float]
        Initial state [S, x, y, z].
    t_end : float
        Final simulation time.
    dt : float
        Time step size.

    Returns
    -------
    times : list[float]
        Time stamps.
    states : list[list[float]]
        State vectors at each time stamp.
    info : dict
        Diagnostics: total_clamps, exploded, steps.
    """
    _validate_inputs(state0, t_end, dt)
    n = len(state0)
    n_steps = int(math.ceil(t_end / dt))
    times = [0.0]
    states = [list(state0)]
    total_clamps = 0
    exploded = False

    for i in range(n_steps):
        t = times[-1]
        y = states[-1]
        step = min(dt, t_end - t)

        k1 = rhs_fn(y)
        y_k2 = [y[j] + 0.5 * step * k1[j] for j in range(n)]
        k2 = rhs_fn(y_k2)
        y_k3 = [y[j] + 0.5 * step * k2[j] for j in range(n)]
        k3 = rhs_fn(y_k3)
        y_k4 = [y[j] + step * k3[j] for j in range(n)]
        k4 = rhs_fn(y_k4)

        new_y = [
            y[j]
            + (step / 6.0)
            * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
            for j in range(n)
        ]

        new_y, clamps = _clamp(new_y)
        total_clamps += clamps

        if _check_explode(new_y):
            exploded = True
            times.append(t + step)
            states.append(new_y)
            break

        times.append(t + step)
        states.append(new_y)

    return times, states, {
        "total_clamps": total_clamps,
        "exploded": exploded,
        "steps": len(times) - 1,
    }
