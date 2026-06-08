"""Right-hand side of the chemostat food-chain ODE system.

Model:
    dS/dt = 1 - S - f1(S)*x
    dx/dt = x*(f1(S) - D1) - f2(x)*y
    dy/dt = y*(f2(x) - D2) - f3(y)*z
    dz/dt = z*(f3(y) - D3)

State vector: [S, x, y, z]
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "equilibria_stability"))
from model import Params


def rhs(state: list[float], params: Params) -> list[float]:
    """Compute derivatives [dS/dt, dx/dt, dy/dt, dz/dt].

    Parameters
    ----------
    state : list[float]
        Current state [S, x, y, z].  Values may be negative (clamping
        is handled by the solver, not here).
    params : Params
        Model parameters containing Monod functions f1, f2, f3 and
        removal rates D1, D2, D3.

    Returns
    -------
    list[float]
        Derivatives [dS/dt, dx/dt, dy/dt, dz/dt].
    """
    S, x, y, z = state

    f1_S = params.f1.value(S)
    f2_x = params.f2.value(x)
    f3_y = params.f3.value(y)

    dS = 1.0 - S - f1_S * x
    dx = x * (f1_S - params.D1) - f2_x * y
    dy = y * (f2_x - params.D2) - f3_y * z
    dz = z * (f3_y - params.D3)

    return [dS, dx, dy, dz]
