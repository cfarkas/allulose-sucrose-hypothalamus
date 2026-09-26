#!/usr/bin/env python3
"""Numerical and replication safeguards for the Figure 3 statistics module."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

_spec = spec_from_file_location("fig3_stats", Path(__file__).with_name("fig3_keyence_permanova.py"))
stats = module_from_spec(_spec)
_spec.loader.exec_module(stats)


class SpatialStatisticsTests(unittest.TestCase):
    def test_exact_allocation_counts_and_block_preservation(self):
        labels = np.array(["Water"] * 3 + ["Sucrose"] * 4 + ["Allulose"] * 5)
        self.assertEqual(stats.allocation_count(labels), 27720)
        self.assertEqual(sum(1 for _ in stats.allocations(labels)), 27720)
        sugar_labels = np.array(["Allulose"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 2 + ["Sucrose"])
        blocks = np.array(["Zeiss_single"] * 6 + ["Keyence_two"] * 3)
        assigned = list(stats.allocations(sugar_labels, blocks))
        self.assertEqual(len(assigned), 60)
        self.assertEqual(sum(np.array_equal(row, sugar_labels) for row in assigned), 1)
        self.assertEqual(len({tuple(row) for row in assigned}), 60)
        for row in assigned:
            for block in np.unique(blocks):
                np.testing.assert_array_equal(np.sort(row[blocks == block]), np.sort(sugar_labels[blocks == block]))

    def test_pseudo_f_matches_distance_matrix_definition(self):
        matrix = np.array([[0, 1], [1, 2], [3, 1], [4, 2], [2, 5], [3, 6]], dtype=float)
        labels = np.array(["W", "W", "S", "S", "A", "A"])
        squared_distances = np.square(matrix[:, None] - matrix[None, :]).sum(axis=2)
        total = squared_distances.sum() / (2 * len(matrix))
        within = sum(squared_distances[np.ix_(labels == group, labels == group)].sum() / (2 * (labels == group).sum()) for group in np.unique(labels))
        expected_f = ((total - within) / 2) / (within / 3)
        statistic, rsquared = stats.pseudo_f(matrix, labels)
        self.assertAlmostEqual(statistic, expected_f, places=12)
        self.assertAlmostEqual(rsquared, (total - within) / total, places=12)

    def test_partial_statistic_matches_independent_residualized_regression(self):
        frame = pd.DataFrame({"platform": ["Zeiss"] * 6 + ["Keyence"] * 3,
                              "biological_cohort": ["single"] * 6 + ["two"] * 3})
        labels = np.array(["Allulose"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 2 + ["Sucrose"])
        treatment = (labels == "Allulose").astype(float)
        noise = np.array([[.1, .2], [-.4, .5], [.7, -.6], [.8, .1], [-.3, .9], [.2, -.7], [.5, -.3], [-.6, .4], [.3, .8]])
        matrix = noise + treatment[:, None] * np.array([.7, -.2]) + (frame.platform == "Keyence").to_numpy()[:, None] * np.array([100., -30.])
        reduced = stats.nuisance_design(frame)
        result = stats.partial_f(matrix, reduced, labels)
        # Independent Frisch-Waugh residualization by the known two blocks.
        residual_matrix = matrix.copy()
        residual_treatment = treatment.copy()
        for platform in frame.platform.unique():
            selected = frame.platform.eq(platform).to_numpy()
            residual_matrix[selected] -= matrix[selected].mean(axis=0)
            residual_treatment[selected] -= treatment[selected].mean()
        beta = residual_treatment @ residual_matrix / (residual_treatment @ residual_treatment)
        fitted = residual_treatment[:, None] * beta
        effect_ss = np.square(fitted).sum()
        full_ss = np.square(residual_matrix - fitted).sum()
        expected_f = effect_ss / (full_ss / (len(frame) - 3))
        self.assertAlmostEqual(result["partial_pseudo_f"], expected_f, places=10)
        self.assertAlmostEqual(result["partial_r_squared"], effect_ss / np.square(residual_matrix).sum(), places=10)
        self.assertEqual(result["df_condition"], 1)
        self.assertEqual(result["reduced_rank"], 2)  # platform/cohort columns alias
        self.assertEqual(result["full_rank"], 3)
        self.assertNotAlmostEqual(result["partial_pseudo_f"], stats.pseudo_f(matrix, labels)[0], places=3)

    def test_scaling_is_invariant_to_feature_units_and_origin(self):
        matrix = np.array([[.7, -.2], [.8, -.1], [.9, .4], [1.1, .5]])
        scaled, _, _, _ = stats.standardize(matrix)
        rescaled, _, _, _ = stats.standardize(matrix * [1000, 100] + [52, -4])
        np.testing.assert_allclose(scaled, rescaled, atol=1e-12)
        np.testing.assert_allclose(scaled.std(axis=0, ddof=1), [1, 1])

    def test_geometric_median_handles_vertex_and_symmetric_cases(self):
        points = np.array([[0., 0.], [1., 0.], [0., 1.]])
        median = stats.geometric_median(points)
        expected = (3 - np.sqrt(3)) / 6
        np.testing.assert_allclose(median, [expected, expected], atol=1e-9)
        np.testing.assert_allclose(stats.geometric_median(np.array([[0., 0.], [0., 0.], [4., 5.]])), [0, 0], atol=1e-9)
        np.testing.assert_allclose(stats.geometric_median(np.array([[0., 0.], [2., 4.]])), [1, 2])

    def test_platform_alias_is_not_presented_as_adjusted_water_effect(self):
        frame = pd.DataFrame({"condition": ["Water"] * 3 + ["Sucrose"] * 3 + ["Allulose"] * 3,
                              "platform": ["Leica"] * 3 + ["Zeiss"] * 6,
                              "biological_cohort": ["water"] * 3 + ["single"] * 6})
        diagnostic = stats.design_diagnostic(frame)
        self.assertFalse(diagnostic["all_three_condition_contrasts_identifiable"])
        self.assertEqual(diagnostic["estimable_condition_degrees_of_freedom"], 1)

    def test_holm_retains_unavailable_planned_endpoint(self):
        np.testing.assert_allclose(stats.holm([.01, .03]), [.02, .03])
        result = stats.holm([.02, np.nan])
        self.assertAlmostEqual(result[0], .04)
        self.assertTrue(np.isnan(result[1]))

    def test_end_to_end_missing_positive_animal_and_roster_consistency(self):
        animals = ["W1", "W2", "W3", "S1", "S2", "A1", "A2"]
        conditions = ["Water"] * 3 + ["Sucrose"] * 2 + ["Allulose"] * 2
        roster = pd.DataFrame({"animal_id": animals, "condition": conditions,
                              "platform": ["Leica_SP8"] * 3 + ["Zeiss_legacy"] * 4,
                              "biological_cohort": ["water"] * 3 + ["single"] * 4,
                              "cage": ["unknown"] * 7})
        rows = []
        values = [(1.0, -.1), (1.3, -.2), (1.1, .05), (.8, .3), (1.0, .15), (1.2, .1), (1.1, .4)]
        for endpoint in stats.ENDPOINTS:
            for animal, condition, (core, excess) in zip(animals, conditions, values):
                missing = animal == "W3" and endpoint == "marker_cfos"
                rows.append(dict(cohort="NPY", region="FIELD", animal_id=animal, condition=condition,
                                 endpoint=endpoint, core_enrichment=np.nan if missing else core,
                                 cluster_excess=np.nan if missing else excess, m=0 if missing else 5))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            roster.to_csv(path / "animal_roster.csv", index=False)
            pd.DataFrame(rows).to_csv(path / "animal_metrics.csv", index=False)
            receipt = stats.run(path)
            self.assertEqual(receipt["status"], "PASS")
            result = pd.read_csv(path / "spatial_permanova.csv").set_index(["endpoint", "feature"])
            # One test per measure, so each endpoint contributes two rows.
            self.assertEqual(sorted(result.index.get_level_values("feature").unique()),
                             sorted(stats.FEATURES))
            self.assertEqual(len(result), len(stats.ENDPOINTS) * len(stats.FEATURES))
            for feature in stats.FEATURES:
                self.assertEqual(result.loc[("cfos", feature), "n_water"], 3)
                self.assertEqual(result.loc[("marker_cfos", feature), "n_water"], 2)
                self.assertEqual(result.loc[("marker_cfos", feature), "n_recorded_water"], 3)
                self.assertEqual(result.loc[("cfos", feature), "enumerated_labelings"], 210)
                self.assertEqual(result.loc[("marker_cfos", feature), "enumerated_labelings"], 90)
            self.assertEqual(int(result.family_size.iloc[0]), len(stats.ENDPOINTS) * len(stats.FEATURES))
            self.assertNotIn("p_holm_two_endpoint_family", result.columns)
            self.assertFalse((path / "spatial_dispersion.csv").exists())
            eligibility = pd.read_csv(path / "spatial_feature_eligibility.csv")
            row = eligibility.loc[eligibility.animal_id.eq("W3") & eligibility.endpoint.eq("marker_cfos")].iloc[0]
            self.assertFalse(row.joint_eligible)
            self.assertTrue(np.isnan(row.core_enrichment))
            sensitivity = pd.read_csv(path / "spatial_sugar_platform_sensitivity.csv")
            self.assertTrue(sensitivity.enumerated_labelings.eq(6).all())
            self.assertTrue(sensitivity.n_water.eq(0).all())
            # Conflicting treatment assignments fail instead of silently merging.
            roster.loc[roster.animal_id.eq("A1"), "condition"] = "Sucrose"
            roster.to_csv(path / "animal_roster.csv", index=False)
            with self.assertRaisesRegex(ValueError, "conflicts"):
                stats.read_inputs(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
