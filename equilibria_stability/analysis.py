"""Equilibrium formulas and stability-threshold checks."""

from __future__ import annotations

from model import Params, State
from numerics import bisection


def solve_S_for_x(params: Params, x_value: float) -> float | None:
    """Solve S + f1(S)*x = 1 for S in [0,1]."""
    return bisection(lambda S: S + params.f1.value(S) * x_value - 1.0, 0.0, 1.0)


def compute_P0(params: Params) -> State:
    return (1.0, 0.0, 0.0, 0.0)


def compute_P1(params: Params) -> State | None:
    S1 = params.f1.inverse(params.D1)
    if S1 is None or not (0.0 < S1 < 1.0):
        return None
    x1 = (1.0 - S1) / params.D1
    if x1 <= 0.0:
        return None
    return (S1, x1, 0.0, 0.0)


def compute_P2(params: Params) -> State | None:
    x2 = params.f2.inverse(params.D2)
    if x2 is None or x2 <= 0.0:
        return None
    S2 = solve_S_for_x(params, x2)
    if S2 is None:
        return None
    y2 = x2 * (params.f1.value(S2) - params.D1) / params.D2
    if y2 <= 0.0:
        return None
    return (S2, x2, y2, 0.0)


def find_P3_x(params: Params, y3: float) -> float | None:
    """Find x3 from x(f1(S(x))-D1)-f2(x)y3=0."""

    def residual(x_value: float) -> float | None:
        S_value = solve_S_for_x(params, x_value)
        if S_value is None:
            return None
        return x_value * (params.f1.value(S_value) - params.D1) - params.f2.value(x_value) * y3

    previous_x = 1e-8
    previous_value = residual(previous_x)
    for i in range(1, 900):
        x_value = 10.0 ** (-8.0 + i * 12.0 / 899.0)
        value = residual(x_value)
        if value is None or previous_value is None:
            previous_x = x_value
            previous_value = value
            continue
        if previous_value * value < 0.0:
            return bisection(lambda xx: residual(xx) or 0.0, previous_x, x_value)
        previous_x = x_value
        previous_value = value
    return None


def compute_P3(params: Params) -> State | None:
    y3 = params.f3.inverse(params.D3)
    if y3 is None or y3 <= 0.0:
        return None
    x3 = find_P3_x(params, y3)
    if x3 is None or x3 <= 0.0:
        return None
    S3 = solve_S_for_x(params, x3)
    if S3 is None:
        return None
    z3 = y3 * (params.f2.value(x3) - params.D2) / params.D3
    if z3 <= 0.0:
        return None
    return (S3, x3, y3, z3)


def compute_all_equilibria(params: Params) -> dict[str, State | None]:
    return {
        "P0": compute_P0(params),
        "P1": compute_P1(params),
        "P2": compute_P2(params),
        "P3": compute_P3(params),
    }


def p0_stability(params: Params) -> tuple[bool, float]:
    threshold = params.f1.value(1.0) - params.D1
    return threshold < 0.0, threshold
