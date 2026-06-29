"""Shared helpers for the theorem-validation experiments.

This module deliberately reuses the *existing* project code so that the
validation experiments describe exactly the same system as the rest of the
repository:

* the model, equilibrium solver and integrator come from the
  ``robustness_sensitivity/chemostat.py`` re-export shim
  (which itself re-exports ``equilibria_stability`` and ``simulator``);
* the analytic Jacobian and its eigenvalues are computed here with the same
  pure-stdlib routines used by ``demo/app.py`` (Faddeev--LeVerrier for the
  characteristic polynomial, Durand--Kerner for its roots).  These realise the
  Routh--Hurwitz / Hopf eigenvalue analysis of Sections 2 and 4 of the paper.

Nothing in the existing project is modified; this file only *imports* it.
"""

from __future__ import annotations

import csv
import math
import os
import sys

# --- Make the project's own modules importable ---------------------------------
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(THIS_DIR)
ROBUSTNESS_DIR = os.path.join(PROJECT_ROOT, "robustness_sensitivity")
EQUILIBRIA_DIR = os.path.join(PROJECT_ROOT, "equilibria_stability")
for _p in (PROJECT_ROOT, ROBUSTNESS_DIR, EQUILIBRIA_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import chemostat as cx                                     # noqa: E402
from chemostat import (                                    # noqa: E402,F401
    PARAMETER_NAMES, Params, State, params_from_dict,
    compute_all_equilibria, compute_P0, compute_P1, compute_P2, compute_P3,
    p0_stability, integrate, classify_regime, in_invariant_region,
)
# Local-stability primitives are shared with demo/app.py (single source of truth).
from stability import (                                    # noqa: E402,F401
    jacobian, charpoly, poly_roots, eigenvalues, classify_eigs,
)

OUT_DIR = os.path.join(THIS_DIR, "outputs")
DATA_DIR = os.path.join(OUT_DIR, "data")
FIG_DIR = os.path.join(OUT_DIR, "figures")

# Full catalogue of fixed parameter sets (no random data) shipped with the
# project.  The richer 13-scenario file is the most convenient source.
SCENARIO_CSV = os.path.join(EQUILIBRIA_DIR, "data", "scenario_parameters.csv")


# --- Scenario input ------------------------------------------------------------

def read_scenarios(path: str = SCENARIO_CSV) -> dict[str, dict[str, float]]:
    """Return {scenario_name: {param: value}} from a parameter CSV."""
    scenarios: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = row["scenario"]
            scenarios[name] = {key: float(row[key]) for key in PARAMETER_NAMES}
    return scenarios


def with_override(base: dict[str, float], name: str, value: float) -> dict[str, float]:
    updated = dict(base)
    updated[name] = value
    return updated


def linspace(low: float, high: float, num: int) -> list[float]:
    """Inclusive evenly spaced values (pure-Python, no numpy dependency)."""
    if num < 2:
        return [low]
    step = (high - low) / (num - 1)
    return [low + step * i for i in range(num)]


# Jacobian / charpoly / poly_roots / eigenvalues / classify_eigs are imported
# from equilibria_stability/stability.py above (single shared source of truth).


# --- small numeric utilities ---------------------------------------------------

def distance(a: State, b: State) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(4)))


def write_csv(rows: list[dict[str, object]], fields: list[str], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
