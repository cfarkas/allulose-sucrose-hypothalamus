#!/usr/bin/env python3
"""Numerical checks for the scale-resolved curves and the global envelope test."""
from __future__ import annotations

import sys
import unittest
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import f_oneway

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig3_scale_resolved as scale


def brute_force_moments(xy, n_positive, radii):
    """Mean and variance of positive-pair counts over every labelling."""
    n = len(xy)
    distance = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=-1)
    means, variances = [], []
    for radius in radii:
        values = []
        for chosen in combinations(range(n), n_positive):
            mark = np.zeros(n, dtype=bool)
            mark[list(chosen)] = True
            pairs = sum(1 for i, j in combinations(range(n), 2)
                        if mark[i] and mark[j] and distance[i, j] <= radius)
            values.append(pairs)
        means.append(np.mean(values))
        variances.append(np.var(values))
    return np.array(means), np.array(variances)


class ScaleResolvedTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20260924)
        self.xy = rng.uniform(0, 120, size=(9, 2))
        self.radii = np.array([30.0, 60.0, 90.0, 150.0])

    def test_exact_moments_match_enumeration(self):
        for n_positive in (2, 3, 4, 5):
            labels = np.zeros(len(self.xy), dtype=bool)
            labels[:n_positive] = True
            _, expected, variance, _ = scale.pair_statistics(self.xy, labels, self.radii)
            mean_ref, variance_ref = brute_force_moments(self.xy, n_positive, self.radii)
            np.testing.assert_allclose(expected, mean_ref, rtol=0, atol=1e-9)
            np.testing.assert_allclose(variance, variance_ref, rtol=0, atol=1e-9)

    def test_observed_count_is_the_actual_pair_count(self):
        labels = np.zeros(len(self.xy), dtype=bool)
        labels[[0, 1, 2, 3]] = True
        observed, _, _, pairs = scale.pair_statistics(self.xy, labels, self.radii)
        distance = np.linalg.norm(self.xy[:, None, :] - self.xy[None, :, :], axis=-1)
        for index, radius in enumerate(self.radii):
            direct = sum(1 for i, j in combinations(range(len(self.xy)), 2)
                         if labels[i] and labels[j] and distance[i, j] <= radius)
            self.assertEqual(observed[index], direct)
            total = sum(1 for i, j in combinations(range(len(self.xy)), 2) if distance[i, j] <= radius)
            self.assertEqual(pairs[index], total)

    def test_fewer_than_two_positives_gives_no_variance(self):
        labels = np.zeros(len(self.xy), dtype=bool)
        labels[0] = True
        observed, expected, variance, _ = scale.pair_statistics(self.xy, labels, self.radii)
        self.assertTrue(np.all(observed == 0) and np.all(expected == 0) and np.all(variance == 0))

    def test_one_way_f_matches_scipy(self):
        rng = np.random.default_rng(7)
        values = rng.normal(size=(9, 4))
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3)
        mine = scale.one_way_f(values, labels)
        for index in range(values.shape[1]):
            reference = f_oneway(values[labels == "Water", index], values[labels == "Sucrose", index],
                                 values[labels == "Allulose", index])
            self.assertAlmostEqual(mine[index], float(reference.statistic), places=10)

    def test_envelope_enumerates_every_allocation(self):
        rng = np.random.default_rng(3)
        values = rng.normal(size=(9, 4))
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3)
        outcome = scale.envelope_test(values, labels, np.arange(4))
        self.assertEqual(outcome["allocations"], comb(9, 3) * comb(6, 3))
        self.assertGreaterEqual(outcome["p_value"], outcome["minimum_attainable_p"])
        self.assertLessEqual(outcome["p_value"], 1.0)

    def test_envelope_is_calibrated_under_exchangeability(self):
        """Exchangeable curves must not be flagged; the p is uniform by construction."""
        rng = np.random.default_rng(11)
        values = rng.normal(size=(9, 5))
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3)
        self.assertGreater(scale.envelope_test(values, labels, np.arange(5))["p_value"], 0.05)

    def test_envelope_detects_a_separated_group(self):
        """Three mutually separated groups make the observed allocation uniquely extreme."""
        values = np.vstack([np.full((3, 5), 12.0), np.full((3, 5), 6.0), np.zeros((3, 5))])
        values += np.random.default_rng(5).normal(scale=.01, size=values.shape)
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3)
        outcome = scale.envelope_test(values, labels, np.arange(5))
        self.assertAlmostEqual(outcome["p_value"], outcome["minimum_attainable_p"], places=12)

    def test_tied_allocations_are_not_counted_as_evidence(self):
        """Two identical groups leave many allocations as extreme as the observed one."""
        values = np.vstack([np.full((3, 5), 6.0), np.zeros((3, 5)), np.zeros((3, 5))])
        values += np.random.default_rng(5).normal(scale=.01, size=values.shape)
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3)
        outcome = scale.envelope_test(values, labels, np.arange(5))
        self.assertGreater(outcome["p_value"], outcome["minimum_attainable_p"])
        self.assertLess(outcome["p_value"], .05)

    def test_critical_curve_separates_rejection(self):
        """The observed curve exits its critical curve exactly when the test rejects."""
        values = np.vstack([np.full((3, 5), 6.0), np.zeros((3, 5)), np.zeros((3, 5))])
        values += np.random.default_rng(9).normal(scale=.01, size=values.shape)
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3)
        outcome = scale.envelope_test(values, labels, np.arange(5))
        exits = np.any(outcome["observed_f"] > outcome["critical_curve"])
        self.assertEqual(bool(exits), outcome["p_value"] <= outcome["alpha"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
