
from __future__ import annotations

from model import Params, State


def jacobian(state: State, params: Params) -> list[list[float]]:
    """Analytic Jacobian of the model at ``state`` (matches simulator/rhs.py)."""
    S, x, y, z = state
    f1, f1p = params.f1.value(S), params.f1.derivative(S)
    f2, f2p = params.f2.value(x), params.f2.derivative(x)
    f3, f3p = params.f3.value(y), params.f3.derivative(y)
    return [
        [-1.0 - f1p * x, -f1,                       0.0,                      0.0],
        [ f1p * x,        f1 - params.D1 - f2p * y, -f2,                      0.0],
        [ 0.0,            f2p * y,                   f2 - params.D2 - f3p * z, -f3],
        [ 0.0,            0.0,                       f3p * z,                  f3 - params.D3],
    ]


def charpoly(A: list[list[float]]) -> list[float]:
    """Characteristic-polynomial coefficients [1, c1, ..., cn] (Faddeev--LeVerrier)."""
    n = len(A)
    I = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    def matmul(P, Q):
        return [[sum(P[i][k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    def add_diag(P, s):
        return [[P[i][j] + (s if i == j else 0.0) for j in range(n)] for i in range(n)]

    def trace(P):
        return sum(P[i][i] for i in range(n))

    coeffs = [1.0]
    M = I
    for k in range(1, n + 1):
        AM = matmul(A, M)
        ck = -trace(AM) / k
        coeffs.append(ck)
        M = add_diag(AM, ck)
    return coeffs


def poly_roots(coeffs: list[float]) -> list[complex]:
    """All (complex) roots of a monic real polynomial via Durand--Kerner."""
    n = len(coeffs) - 1
    if n <= 0:
        return []

    def evalp(z):
        r = 0j
        for c in coeffs:
            r = r * z + c
        return r

    roots = [(0.4 + 0.9j) ** k for k in range(n)]
    for _ in range(500):
        moved = 0.0
        for i in range(n):
            num = evalp(roots[i])
            den = 1 + 0j
            for j in range(n):
                if j != i:
                    den *= (roots[i] - roots[j])
            if den == 0:
                den = 1e-12
            delta = num / den
            roots[i] -= delta
            moved = max(moved, abs(delta))
        if moved < 1e-12:
            break
    return roots


def eigenvalues(state: State, params: Params) -> list[complex]:
    """Eigenvalues of the Jacobian at ``state``."""
    return poly_roots(charpoly(jacobian(state, params)))


def classify_eigs(eigs: list[complex], eps: float = 1e-7) -> tuple[str, float]:
    """Plain-language local-stability label + max real part."""
    max_re = max(e.real for e in eigs)
    has_pos = any(e.real > eps for e in eigs)
    has_neg = any(e.real < -eps for e in eigs)
    spiral = any(abs(e.imag) > 1e-6 for e in eigs)
    if max_re < -eps:
        kind = "stable focus" if spiral else "stable node"
    elif max_re > eps:
        if has_pos and has_neg:
            kind = "saddle"
        else:
            kind = "unstable focus" if spiral else "unstable node"
    else:
        kind = "marginal"
    return kind, max_re
