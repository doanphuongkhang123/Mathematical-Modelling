from __future__ import annotations

import csv
import math
import tempfile
import unittest
from pathlib import Path

from simulator.runner import write_csv
from simulator.solvers import euler, rk4


class SolverTests(unittest.TestCase):
    @staticmethod
    def exponential_rhs(state: list[float]) -> list[float]:
        return [-state[0]]

    def test_euler_matches_its_recurrence(self) -> None:
        times, states, info = euler(self.exponential_rhs, [1.0], t_end=1.0, dt=0.1)
        self.assertEqual(info["steps"], 10)
        self.assertFalse(info["exploded"])
        self.assertEqual(info["total_clamps"], 0)
        self.assertAlmostEqual(times[-1], 1.0)
        self.assertAlmostEqual(states[-1][0], 0.9**10, places=14)

    def test_rk4_is_more_accurate_than_euler(self) -> None:
        _, euler_states, _ = euler(self.exponential_rhs, [1.0], t_end=1.0, dt=0.1)
        _, rk4_states, _ = rk4(self.exponential_rhs, [1.0], t_end=1.0, dt=0.1)
        exact = math.exp(-1.0)
        euler_error = abs(euler_states[-1][0] - exact)
        rk4_error = abs(rk4_states[-1][0] - exact)
        self.assertLess(rk4_error, euler_error)
        self.assertLess(rk4_error, 1e-6)

    def test_rk4_exhibits_fourth_order_convergence(self) -> None:
        errors = []
        for dt in (0.2, 0.1):
            _, states, _ = rk4(self.exponential_rhs, [1.0], t_end=1.0, dt=dt)
            errors.append(abs(states[-1][0] - math.exp(-1.0)))
        self.assertGreater(errors[0] / errors[1], 12.0)

    def test_final_partial_step_stops_exactly_at_t_end(self) -> None:
        for solver in (euler, rk4):
            with self.subTest(solver=solver.__name__):
                times, _, info = solver(self.exponential_rhs, [1.0], t_end=1.0, dt=0.3)
                self.assertEqual(info["steps"], 4)
                self.assertAlmostEqual(times[-1], 1.0, places=14)
                self.assertTrue(all(left < right for left, right in zip(times, times[1:])))

    def test_negative_values_are_clamped(self) -> None:
        times, states, info = euler(lambda state: [-10.0], [1.0], t_end=0.2, dt=0.2)
        self.assertEqual(times, [0.0, 0.2])
        self.assertEqual(states[-1], [0.0])
        self.assertEqual(info["total_clamps"], 1)

    def test_explosion_stops_integration(self) -> None:
        times, states, info = euler(lambda state: [2e10], [0.0], t_end=10.0, dt=1.0)
        self.assertTrue(info["exploded"])
        self.assertEqual(info["steps"], 1)
        self.assertEqual(len(times), 2)
        self.assertGreater(states[-1][0], 1e10)

    def test_initial_state_is_not_mutated(self) -> None:
        for solver in (euler, rk4):
            with self.subTest(solver=solver.__name__):
                state0 = [1.0]
                solver(self.exponential_rhs, state0, t_end=0.2, dt=0.1)
                self.assertEqual(state0, [1.0])

    def test_invalid_inputs_raise_value_error(self) -> None:
        for solver in (euler, rk4):
            with self.subTest(solver=solver.__name__):
                with self.assertRaises(ValueError):
                    solver(self.exponential_rhs, [1.0], t_end=1.0, dt=0.0)
                with self.assertRaises(ValueError):
                    solver(self.exponential_rhs, [1.0], t_end=-1.0, dt=0.1)
                with self.assertRaises(ValueError):
                    solver(self.exponential_rhs, [], t_end=1.0, dt=0.1)

    def test_csv_output_has_required_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.csv"
            times = [0.0, 0.1]
            states = [[1.0, 0.2, 0.1, 0.05], [0.9, 0.3, 0.08, 0.04]]
            write_csv(str(path), times, states)
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
        self.assertEqual(rows[0], ["t", "S", "x", "y", "z"])
        self.assertEqual(len(rows), 3)
        self.assertEqual(
            rows[-1],
            ["0.100000", "0.900000", "0.300000", "0.080000", "0.040000"],
        )


if __name__ == "__main__":
    unittest.main()
