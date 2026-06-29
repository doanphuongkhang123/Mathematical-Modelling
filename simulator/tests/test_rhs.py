from __future__ import annotations

import csv
import unittest
from pathlib import Path

from simulator.rhs import rhs

from analysis import compute_all_equilibria
from model import PARAMETER_NAMES, Monod, Params, params_from_dict


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCENARIO_CSV = PROJECT_ROOT / "equilibria_stability" / "data" / "scenario_parameters.csv"


def load_scenario(name: str) -> Params:
    with SCENARIO_CSV.open(newline="", encoding="utf-8") as handle:
        row = next(row for row in csv.DictReader(handle) if row["scenario"] == name)
    return params_from_dict({key: float(row[key]) for key in PARAMETER_NAMES})


class RhsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.params = Params(
            f1=Monod(2.0, 1.0),
            f2=Monod(3.0, 2.0),
            f3=Monod(4.0, 3.0),
            D1=0.5,
            D2=0.6,
            D3=0.7,
        )

    def test_known_derivative(self) -> None:
        derivative = rhs([1.0, 0.5, 0.25, 0.125], self.params)
        expected = [-0.5, 0.1, -0.038461538461538464, -0.04903846153846154]
        for actual, target in zip(derivative, expected):
            self.assertAlmostEqual(actual, target, places=12)

    def test_rhs_does_not_mutate_state(self) -> None:
        state = [0.7, 0.4, 0.2, 0.1]
        original = list(state)
        rhs(state, self.params)
        self.assertEqual(state, original)

    def test_vector_field_points_into_nonnegative_region(self) -> None:
        self.assertGreater(rhs([0.0, 0.4, 0.2, 0.1], self.params)[0], 0.0)
        self.assertEqual(rhs([0.5, 0.0, 0.2, 0.1], self.params)[1], 0.0)
        self.assertEqual(rhs([0.5, 0.2, 0.0, 0.1], self.params)[2], 0.0)
        self.assertEqual(rhs([0.5, 0.2, 0.1, 0.0], self.params)[3], 0.0)

    def test_computed_equilibria_are_stationary(self) -> None:
        cases = {
            "washout_P0": "P0",
            "prey_only_P1": "P1",
            "prey_first_predator_P2": "P2",
            "coexistence_P3": "P3",
        }
        for scenario, equilibrium_name in cases.items():
            with self.subTest(scenario=scenario):
                params = load_scenario(scenario)
                equilibrium = compute_all_equilibria(params)[equilibrium_name]
                self.assertIsNotNone(equilibrium)
                derivative = rhs(list(equilibrium), params)  # type: ignore[arg-type]
                self.assertLess(max(abs(value) for value in derivative), 1e-8)


if __name__ == "__main__":
    unittest.main()
