"""Model parameter objects and Monod response functions."""

from __future__ import annotations

from dataclasses import dataclass


PARAMETER_NAMES = ("a1", "k1", "a2", "k2", "a3", "k3", "D1", "D2", "D3")
State = tuple[float, float, float, float]


@dataclass(frozen=True)
class Monod:
    """Monod response f(u)=a*u/(k+u)."""

    a: float
    k: float

    def value(self, u: float) -> float:
        u = max(0.0, u)
        if u == 0.0:
            return 0.0
        return self.a * u / (self.k + u)

    def derivative(self, u: float) -> float:
        u = max(0.0, u)
        return self.a * self.k / ((self.k + u) ** 2)

    def inverse(self, target: float) -> float | None:
        """Return u satisfying f(u)=target, or None if no finite positive root exists."""
        if target < 0.0 or target >= self.a:
            return None
        if target == 0.0:
            return 0.0
        return target * self.k / (self.a - target)


@dataclass(frozen=True)
class Params:
    f1: Monod
    f2: Monod
    f3: Monod
    D1: float
    D2: float
    D3: float


def params_from_dict(data: dict[str, float]) -> Params:
    missing = [name for name in PARAMETER_NAMES if name not in data]
    if missing:
        raise ValueError(f"Missing parameters: {', '.join(missing)}")
    for name in PARAMETER_NAMES:
        if data[name] <= 0.0:
            raise ValueError(f"{name} must be positive, got {data[name]}")
    return Params(
        f1=Monod(data["a1"], data["k1"]),
        f2=Monod(data["a2"], data["k2"]),
        f3=Monod(data["a3"], data["k3"]),
        D1=data["D1"],
        D2=data["D2"],
        D3=data["D3"],
    )
