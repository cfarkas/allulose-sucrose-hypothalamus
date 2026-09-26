#!/usr/bin/env python3
"""Numerical checks for the unequal-variance Figure 3 abundance tests."""
from __future__ import annotations

import sys
import unittest
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import levene, ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig3_robust_abundance as robust

WATER = np.array([8.58, 6.66, 1.64, 0.34, 0.87])
SUCROSE = np.array([1.40, 2.50, 1.42])
ALLULOSE = np.array([5.67, 5.54, 9.05])


class RobustAbundanceTests(unittest.TestCase):
    def test_mann_whitney_null_matches_enumeration(self):
        for left, right in ((2, 2), (3, 3), (4, 3), (5, 3), (3, 5), (2, 5), (6, 4)):
            counts = robust.mann_whitney_null_counts(left, right)
            brute = np.zeros(left * right + 1)
            pooled = list(range(left + right))
            for chosen in combinations(pooled, left):
                selected = set(chosen)
                remainder = [value for value in pooled if value not in selected]
                brute[sum(1 for a in selected for b in remainder if a > b)] += 1
            np.testing.assert_array_equal(counts, brute)
            self.assertEqual(counts.sum(), comb(left + right, left))

    def test_welch_anova_reduces_to_welch_t_squared(self):
        result = robust.welch_anova([WATER, SUCROSE])
        reference = ttest_ind(WATER, SUCROSE, equal_var=False)
        self.assertAlmostEqual(result["statistic"], reference.statistic ** 2, places=10)
        self.assertAlmostEqual(result["p_value"], float(reference.pvalue), places=12)

    def test_brown_forsythe_matches_levene_on_medians(self):
        result = robust.brown_forsythe([WATER, SUCROSE, ALLULOSE])
        reference = levene(WATER, SUCROSE, ALLULOSE, center="median")
        self.assertAlmostEqual(result["statistic"], float(reference.statistic), places=10)
        self.assertAlmostEqual(result["p_value"], float(reference.pvalue), places=12)

    def test_welch_anova_is_insensitive_to_a_variance_only_change(self):
        """Equal means with very different spreads must not produce a small p."""
        left = np.array([1.0, 9.0, 1.0, 9.0])
        right = np.array([4.9, 5.0, 5.1])
        self.assertGreater(robust.welch_anova([left, right])["p_value"], 0.5)

    def test_welch_anova_reports_a_degenerate_group_instead_of_a_number(self):
        """A group with no within-group variance gives an infinite Welch weight."""
        result = robust.welch_anova([np.array([1.0, 2.0, 3.0]), np.array([5.0, 5.0, 5.0])])
        self.assertTrue(np.isnan(result["statistic"]))
        self.assertEqual(result["status"], "not_estimable_zero_variance_or_singleton_group")

    def test_studentized_permutation_enumerates_every_allocation(self):
        result = robust.exact_studentized_permutation(WATER, SUCROSE)
        self.assertEqual(result["enumerated_labelings"], comb(len(WATER) + len(SUCROSE), len(WATER)))
        self.assertGreaterEqual(result["extreme_labelings"], 1)
        self.assertLessEqual(result["p_value"], 1.0)

    def test_studentized_permutation_is_exchangeable_under_relabelling(self):
        forward = robust.exact_studentized_permutation(WATER, SUCROSE)
        backward = robust.exact_studentized_permutation(SUCROSE, WATER)
        self.assertAlmostEqual(forward["p_value"], backward["p_value"], places=12)

    def test_studentized_permutation_detects_a_clear_separation(self):
        left = np.array([10.0, 10.1, 9.9, 10.2])
        right = np.array([1.0, 1.1, 0.9])
        self.assertLessEqual(robust.exact_studentized_permutation(left, right)["p_value"], 1 / comb(7, 4) + 1e-12)

    def test_hodges_lehmann_interval_brackets_the_shift(self):
        result = robust.hodges_lehmann(WATER, SUCROSE)
        self.assertAlmostEqual(result["estimate"], float(np.median(
            (WATER[:, None] - SUCROSE[None, :]).ravel())), places=12)
        self.assertEqual(result["ci_status"], "estimable")
        self.assertLessEqual(result["ci_low"], result["estimate"])
        self.assertGreaterEqual(result["ci_high"], result["estimate"])
        self.assertGreater(result["ci_level"], 0.9)

    def test_three_against_three_admits_no_exact_interval(self):
        """With 20 allocations the smallest attainable two-sided tail is 0.10."""
        result = robust.hodges_lehmann(ALLULOSE, SUCROSE)
        self.assertEqual(result["ci_status"], "not_resolvable_at_this_group_size")
        self.assertTrue(np.isnan(result["ci_low"]) and np.isnan(result["ci_high"]))
        self.assertFalse(np.isnan(result["estimate"]))

    def test_hodges_lehmann_shifts_with_the_data(self):
        base = robust.hodges_lehmann(ALLULOSE, SUCROSE)["estimate"]
        shifted = robust.hodges_lehmann(ALLULOSE + 3.0, SUCROSE)["estimate"]
        self.assertAlmostEqual(shifted - base, 3.0, places=12)

    def test_cliffs_delta_bounds_and_sign(self):
        self.assertAlmostEqual(robust.cliffs_delta(np.array([2.0, 3.0]), np.array([0.0, 1.0])), 1.0)
        self.assertAlmostEqual(robust.cliffs_delta(np.array([0.0, 1.0]), np.array([2.0, 3.0])), -1.0)
        self.assertAlmostEqual(robust.cliffs_delta(ALLULOSE, ALLULOSE), 0.0)

    def test_permutation_omnibus_enumerates_every_allocation(self):
        result = robust.exact_permutation_welch_anova([WATER, SUCROSE, ALLULOSE])
        expected = comb(11, 5) * comb(6, 3)
        self.assertEqual(result["enumerated_labelings"], expected)
        self.assertGreaterEqual(result["p_value"], result["minimum_attainable_p"])
        self.assertLessEqual(result["p_value"], 1.0)

    def test_permutation_omnibus_is_not_smaller_than_the_resolution_limit(self):
        """Three animals per group cannot resolve a p below the enumeration floor."""
        tight_a = np.array([1.0, 1.01, 0.99])
        tight_b = np.array([50.0, 50.1, 49.9])
        tight_c = np.array([100.0, 100.1, 99.9])
        result = robust.exact_permutation_welch_anova([tight_a, tight_b, tight_c])
        self.assertAlmostEqual(result["p_value"], result["minimum_attainable_p"], places=12)
        self.assertGreater(result["p_value"], 0.0)

    def test_permutation_omnibus_is_calibrated_on_exchangeable_data(self):
        """Under full exchangeability the enumerated p is uniform, so it is rarely small."""
        values = np.array([3.0, 1.0, 4.0, 1.5, 5.0, 9.0, 2.0, 6.0, 5.5])
        result = robust.exact_permutation_welch_anova([values[:3], values[3:6], values[6:]])
        self.assertGreater(result["p_value"], 0.05)

    def test_studentized_permutation_reports_its_resolution_limit(self):
        three = robust.exact_studentized_permutation(SUCROSE, ALLULOSE)
        self.assertAlmostEqual(three["minimum_attainable_p"], 2 / comb(6, 3), places=12)
        five = robust.exact_studentized_permutation(WATER, SUCROSE)
        self.assertAlmostEqual(five["minimum_attainable_p"], 1 / comb(8, 5), places=12)

    def test_small_groups_are_reported_as_not_estimable(self):
        result = robust.exact_studentized_permutation(np.array([1.0]), np.array([2.0, 3.0]))
        self.assertEqual(result["status"], "fewer_than_two_animals_in_a_group")


if __name__ == "__main__":
    unittest.main(verbosity=2)
