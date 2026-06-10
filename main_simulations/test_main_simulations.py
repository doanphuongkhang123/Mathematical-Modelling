"""Focused tests for the main simulation workflow."""

from __future__ import annotations

import unittest

from main_simulations.run_all import (
    count_peaks,
    euclidean,
    oscillation_criteria,
    sample_indices,
    window_stats,
)


class SimulationHelperTests(unittest.TestCase):
    def test_sample_indices_keep_first_and_last(self):
        self.assertEqual(sample_indices(12, 5), [0, 5, 10, 11])
        self.assertEqual(sample_indices(11, 5), [0, 5, 10])

    def test_euclidean_distance(self):
        self.assertAlmostEqual(euclidean((0, 0, 0, 0), (1, 2, 2, 0)), 3.0)

    def test_window_stats(self):
        means, amplitudes = window_stats(
            [0.0, 1.0, 2.0],
            [[0.0, 1.0, 2.0, 3.0], [2.0, 3.0, 4.0, 5.0], [4.0, 5.0, 6.0, 7.0]],
            1.0,
            2.0,
        )
        self.assertEqual(means, [3.0, 4.0, 5.0, 6.0])
        self.assertEqual(amplitudes, [2.0, 2.0, 2.0, 2.0])

    def test_peak_count(self):
        values = [0.0, 1.0, 0.0, 2.0, 0.0, 1.5, 0.0]
        self.assertEqual(count_peaks(values), 3)

    def test_oscillation_criteria(self):
        passing = oscillation_criteria(4e-4, 8, 7, 0.5)
        self.assertTrue(all(passing.values()))
        failing = oscillation_criteria(2e-3, 3, 2, 0.9)
        self.assertFalse(any(failing.values()))


if __name__ == "__main__":
    unittest.main()
