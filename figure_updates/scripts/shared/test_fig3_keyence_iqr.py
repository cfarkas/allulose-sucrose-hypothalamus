#!/usr/bin/env python3
"""Fence, missingness, and endpoint-specific eligibility safeguards."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import fig3_keyence_iqr as iqr
import fig3_keyence_permanova as spatial


class IQRTests(unittest.TestCase):
    def test_strict_fence_keeps_exact_boundary(self):
        flags, summary = iqr.fence_flags([0, 0, 1, 2, 5])
        self.assertEqual(summary["upper_fence"], 5)
        self.assertFalse(flags.any())
        flags, _ = iqr.fence_flags([0, 0, 1, 2, 5.01])
        np.testing.assert_array_equal(flags, [False, False, False, False, True])

    def test_three_finite_animals_and_missingness(self):
        flags, summary = iqr.fence_flags([0, 1, 1000000, np.nan])
        self.assertEqual(summary["n_finite"], 3)
        self.assertFalse(flags.any())
        flags, summary = iqr.fence_flags([np.nan, np.inf])
        self.assertFalse(flags.any())
        self.assertEqual(summary["n_finite"], 0)
        self.assertTrue(np.isnan(summary["q1"]))

    def test_endpoint_specific_union_does_not_exclude_other_endpoint(self):
        frame = pd.DataFrame({"animal_id": ["a", "b", "c", "d", "e"], "condition": ["Allulose"] * 5,
                              "core_enrichment": [1, 1, 1, 1, 5], "cluster_excess": [0, 0, 0, 1, 0]})
        ledger = iqr.feature_ledger(frame, "spatial", "cfos", iqr.SPATIAL_FEATURES)
        other = frame.assign(core_enrichment=1, cluster_excess=0)
        ledger = pd.concat([ledger, iqr.feature_ledger(other, "spatial", "marker_cfos", iqr.SPATIAL_FEATURES)])
        eligibility = iqr.endpoint_eligibility(ledger)
        rejected = eligibility.loc[eligibility.iqr_excluded]
        self.assertEqual(set(rejected.animal_id), {"d", "e"})
        self.assertEqual(set(rejected.endpoint), {"cfos"})
        self.assertTrue(eligibility.loc[eligibility.endpoint.eq("marker_cfos"), "include_in_main"].all())
        self.assertEqual(rejected.set_index("animal_id").loc["d", "flagged_features"], "cluster_excess")
        self.assertEqual(rejected.set_index("animal_id").loc["e", "flagged_features"], "core_enrichment")

    def test_structural_na_is_nonestimable_not_outlier(self):
        frame = pd.DataFrame({"animal_id": ["W1", "W2", "W3"], "condition": ["Water"] * 3,
                              "core_enrichment": [1., 1.1, np.nan], "cluster_excess": [.1, .2, np.nan]})
        eligibility = iqr.endpoint_eligibility(iqr.feature_ledger(frame, "spatial", "marker_cfos", iqr.SPATIAL_FEATURES))
        water3 = eligibility.set_index("animal_id").loc["W3"]
        self.assertFalse(water3.iqr_excluded)
        self.assertTrue(water3.nonestimable)
        self.assertFalse(water3.include_in_main)

    def test_explicit_stats_eligibility_keeps_observed_values(self):
        roster = pd.DataFrame({"animal_id": ["W1", "W2", "S1", "S2", "A1", "A2"],
                               "condition": ["Water"] * 2 + ["Sucrose"] * 2 + ["Allulose"] * 2,
                               "platform": ["Leica"] * 2 + ["Zeiss"] * 4,
                               "biological_cohort": ["water"] * 2 + ["legacy"] * 4, "cage": "unknown"})
        values, flags = [], []
        for endpoint in spatial.ENDPOINTS:
            for i, row in roster.iterrows():
                excluded = endpoint == "cfos" and row.animal_id == "A2"
                values.append(dict(animal_id=row.animal_id, condition=row.condition, cohort="NPY", region="FIELD",
                                   endpoint=endpoint, core_enrichment=1 + i / 10, cluster_excess=i / 20))
                flags.append(dict(animal_id=row.animal_id, endpoint=endpoint, include_in_main=not excluded,
                                  exclusion_reason="IQR exclusion:core_enrichment" if excluded else "included"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            roster.to_csv(path / "animal_roster.csv", index=False)
            pd.DataFrame(values).to_csv(path / "animal_metrics.csv", index=False)
            pd.DataFrame(flags).to_csv(path / "spatial_eligibility.csv", index=False)
            _, result = spatial.read_inputs(path, path / "spatial_eligibility.csv")
            excluded = result.loc[result.animal_id.eq("A2") & result.endpoint.eq("cfos")].iloc[0]
            self.assertFalse(excluded.joint_eligible)
            self.assertEqual(excluded.core_enrichment, 1.5)
            self.assertTrue(result.loc[result.animal_id.eq("A2") & result.endpoint.eq("marker_cfos"), "joint_eligible"].iloc[0])
            self.assertIn("IQR exclusion", excluded.joint_exclusion_reason)

    def test_abundance_exact_mwu_preserves_tie_enumeration(self):
        result = iqr.exact_mwu(np.array([0., 0.]), np.array([1., 1.]))
        self.assertEqual(result["enumerated_labelings"], 6)
        self.assertEqual(result["extreme_labelings"], 2)
        self.assertAlmostEqual(result["p_value"], 1 / 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
