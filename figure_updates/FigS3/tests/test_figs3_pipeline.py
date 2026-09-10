#!/usr/bin/env python3
"""Fast structural tests for Figure S3; the full data gate is 07_validate_figure_s3.py."""

from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


FIGS3 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIGS3))


class Figure6StructureTests(unittest.TestCase):
    def test_metadata_and_water_control(self):
        slides = pd.read_csv(FIGS3 / "source_data" / "Figure_S3_slide_manifest.csv")
        organs = pd.read_csv(
            FIGS3 / "source_data" / "Figure_S3_sample_organ_manifest.csv"
        )
        self.assertEqual(len(slides), 24)
        self.assertEqual(len(organs), 72)
        self.assertEqual(
            slides.treatment.value_counts().to_dict(),
            {"Sucrose": 9, "Allulose": 8, "Water": 7},
        )
        control = slides.loc[slides.sample_id == "m26-016"].iloc[0]
        self.assertEqual(control.treatment, "Water")
        self.assertEqual(control.sample_note, "brain_not_collected")

    def test_export_receipts(self):
        summary = json.loads(
            (FIGS3 / "exports" / "Figure_S3_L0_L2_export_summary.json").read_text()
        )
        manifest = pd.read_csv(
            FIGS3 / "exports" / "Figure_S3_L0_L2_export_manifest.csv"
        )
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(summary["complete_slides"], 24)
        self.assertEqual(len(manifest), 48)
        self.assertEqual(set(manifest.level), {0, 2})
        self.assertTrue((manifest.output_sha256.str.len() == 64).all())

    def test_segmentation_only_gate(self):
        summary = json.loads(
            (FIGS3 / "organ_segmentation" / "segmentation_summary.json").read_text()
        )
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(summary["slides"], 24)
        self.assertIs(summary["image_registration_performed"], False)
        self.assertFalse((FIGS3 / "registration").exists())
        self.assertTrue((FIGS3 / "03_segment_organs.py").is_file())
        self.assertFalse((FIGS3 / "03_segment_register_organs.py").exists())

    def test_representative_is_median_tissue_water_animal(self):
        slides = pd.read_csv(FIGS3 / "source_data" / "Figure_S3_slide_manifest.csv")
        water_ids = sorted(slides.loc[slides.treatment == "Water", "sample_id"])
        fractions = {
            sample_id: json.loads(
                (
                    FIGS3
                    / "organ_segmentation"
                    / sample_id
                    / "segmentation_receipt.json"
                ).read_text()
            )["tissue_fraction"]
            for sample_id in water_ids
        }
        median = float(np.median(list(fractions.values())))
        selected = min(
            water_ids,
            key=lambda sample_id: (abs(fractions[sample_id] - median), sample_id),
        )

        self.assertEqual(selected, "m26-015")
        analysis = json.loads(
            (FIGS3 / "provenance" / "Figure_S3_analysis_receipt.json").read_text()
        )
        selection = analysis["representative_sample_selection"]
        self.assertEqual(analysis["analysis_version"], "figs3-animal-composition-1.8")
        self.assertIs(analysis["image_registration_performed"], False)
        self.assertEqual(analysis["representative_sample_id"], selected)
        self.assertEqual(selection["selected_sample_id"], selected)
        self.assertEqual(selection["eligible_treatment"], "Water")
        self.assertEqual(selection["eligible_sample_count"], 7)
        self.assertAlmostEqual(selection["median_tissue_fraction"], median)

    def test_crop_vectorization_matches_lazyslide(self):
        try:
            from lazyslide.cv import InstanceMap
        except ImportError:
            self.skipTest("Run this equivalence test in the lazyslide311 environment")
        spec = importlib.util.spec_from_file_location(
            "figs3_gpu", FIGS3 / "04_run_histoplus_gpu.py"
        )
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        labels = np.zeros((840, 840), dtype=np.int32)
        labels[20:29, 30:43] = 1
        labels[200:211, 300:315] = 2
        labels[500:520, 505:525] = 3
        labels[520:525, 515:518] = 3
        generator = np.random.default_rng(3)
        probabilities = generator.random((15, 840, 840), dtype=np.float32)
        probabilities /= probabilities.sum(axis=0, keepdims=True)
        old = InstanceMap(
            labels,
            prob_map=probabilities,
            class_names={0: "Background", **module.CLASS_NAMES},
        ).to_polygons(detect_holes=False)
        new = module.vectorize_instance_tile(labels, probabilities, 0, 0, 1)
        self.assertEqual(len(old), len(new))
        for index in range(len(old)):
            self.assertTrue(old.iloc[index].geometry.equals(new.iloc[index].geometry))
            self.assertEqual(old.iloc[index]["class"], new.iloc[index]["class"])
            self.assertAlmostEqual(
                old.iloc[index]["prob"], new.iloc[index]["prob"], places=7
            )


if __name__ == "__main__":
    unittest.main()
