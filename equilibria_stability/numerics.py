"""Small numerical helpers for solving equilibrium equations."""

from __future__ import annotations

from typing import Callable


def bisection(
    func: Callable[[float], float],
    low: float,
    high: float,
    tol: float = 1e-11,
    max_iter: int = 120,
) -> float | None:
    """Solve func(x)=0 on [low, high] when the endpoints bracket a root."""
    f_low = func(low)
    f_high = func(high)
    if abs(f_low) < tol:
        return low
    if abs(f_high) < tol:
        return high
    if f_low * f_high > 0.0:
        return None
    for _ in range(max_iter):
        mid = 0.5 * (low + high)
        f_mid = func(mid)
        if abs(f_mid) < tol or high - low < tol:
            return mid
        if f_low * f_mid <= 0.0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid
    return 0.5 * (low + high)
