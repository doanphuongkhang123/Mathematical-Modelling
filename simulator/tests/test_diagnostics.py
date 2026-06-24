from __future__ import annotations

import unittest

from simulator.diagnostics import diagnose

from model import Monod, Params


class DiagnosticsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.washout_params = Params(
            f1=Monod(0.6, 0.5),
            f2=Monod(1.5, 0.5),
            f3=Monod(1.2, 0.4),
            D1=0.8,
            D2=0.7,
            D3=0.7,
        )

    def test_empty_trajectory_returns_error(self) -> None:
        self.assertEqual(diagnose([], []), {"error": "empty trajectory"})

    def test_summary_contains_final_min_max_and_solver_info(self) -> None:
        times = [0.0, 0.5, 1.0]
        states = [
            [0.5, 0.2, 0.1, 0.05],
            [0.7, 0.1, 0.08, 0.04],
            [0.9, 0.05, 0.06, 0.03],
        ]
        info = {"total_clamps": 2, "exploded": False, "steps": 2}
        result = diagnose(times, states, info)
        self.assertEqual(result["final_state"], states[-1])
        self.assertEqual(result["min_values"], [0.5, 0.05, 0.06, 0.03])
        self.assertEqual(result["max_values"], [0.9, 0.2, 0.1, 0.05])
        self.assertEqual(result["negative_clamps"], 2)
        self.assertFalse(result["exploded"])
        self.assertEqual(result["total_time"], 1.0)
        self.assertEqual(result["total_steps"], 2)

    def test_constant_tail_is_classified_as_converged(self) -> None:
        times = [float(i) for i in range(20)]
        states = [[1.0, 0.0, 0.0, 0.0] for _ in times]
        self.assertTrue(diagnose(times, states)["converged"])

    def test_varying_tail_is_not_classified_as_converged(self) -> None:
        times = [float(i) for i in range(20)]
        states = [[1.0, 0.0, 0.0, 0.0] for _ in times]
        states[-2] = [0.5, 0.0, 0.0, 0.0]
        self.assertFalse(diagnose(times, states)["converged"])

    def test_equilibrium_distance_is_zero_at_washout(self) -> None:
        result = diagnose(
            [0.0, 1.0],
            [[0.8, 0.1, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]],
            params=self.washout_params,
        )
        distances = result["equilibrium_proximity"]
        self.assertEqual(distances["P0"], 0.0)
        self.assertIsNone(distances["P1"])
        self.assertIsNone(distances["P2"])
        self.assertIsNone(distances["P3"])


if __name__ == "__main__":
    unittest.main()
