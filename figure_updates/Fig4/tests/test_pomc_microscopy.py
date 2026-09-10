"""Regression checks that ROI views display measured microscopy, not label paint."""
import importlib.util
import os
from pathlib import Path
import sys
import unittest

import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(os.environ.get('POMC_TEST_PAPER_ROOT', Path(__file__).resolve().parents[2]))
old_root = os.environ.get('FIG4_PAPER_ROOT')
os.environ['FIG4_PAPER_ROOT'] = str(PAPER)
spec = importlib.util.spec_from_file_location('fig4_microscopy', PAPER/'Fig4/03_make_figure_4_pomc_cfos.py')
renderer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = renderer
spec.loader.exec_module(renderer)
if old_root is None:
    os.environ.pop('FIG4_PAPER_ROOT', None)
else:
    os.environ['FIG4_PAPER_ROOT'] = old_root


class MicroscopyROIView(unittest.TestCase):
    def setUp(self):
        self.channel = np.tile(np.arange(256, dtype=np.uint8), (8, 1))
        self.labels = np.zeros_like(self.channel, dtype=np.int32)
        self.labels[2:6, 0:220] = 17

    def test_dark_roi_stays_dark_and_native_texture_is_retained(self):
        original = self.channel.copy()
        labels = self.labels.copy()
        rgb = renderer.roi_microscopy_view(self.channel, self.labels, 'POMC')
        self.assertTrue(np.all(rgb[2:6, 0] == 0))
        self.assertGreater(len(np.unique(rgb[3, :, 0])), 200)
        self.assertGreater(rgb[3, 190, 0], rgb[3, 90, 0])
        self.assertTrue(np.all(rgb[self.labels == 0] == 0))
        np.testing.assert_array_equal(rgb[self.labels > 0], renderer.single_channel_view(self.channel, 'POMC')[self.labels > 0])
        np.testing.assert_array_equal(self.channel, original)
        np.testing.assert_array_equal(self.labels, labels)

    def test_another_roi_cannot_change_the_display_of_existing_pixels(self):
        first = renderer.roi_microscopy_view(self.channel, self.labels, 'POMC')
        expanded = self.labels.copy()
        expanded[0:2, 200:] = 62
        second = renderer.roi_microscopy_view(self.channel, expanded, 'POMC')
        np.testing.assert_array_equal(first[self.labels > 0], second[self.labels > 0])
        relabeled = np.where(self.labels > 0, 999, 0)
        np.testing.assert_array_equal(first, renderer.roi_microscopy_view(self.channel, relabeled, 'POMC'))

    def test_inset_uses_source_intensity_instead_of_filling_a_dark_roi(self):
        # All-zero channels must produce a black image even with nonempty ROIs.
        zero = np.zeros((32, 256), dtype=np.uint8)
        labels = np.zeros_like(zero, dtype=np.int32)
        labels[10:20, 20:180] = 7
        tissue = np.ones_like(zero, dtype=bool)
        tissue[:3] = False
        lumen = np.zeros_like(tissue); lumen[:20, 120:130] = True
        fig, ax = plt.subplots()
        try:
            receipt = renderer.add_marker_distribution_inset(ax, zero, zero, labels, labels, tissue, lumen, 1.0, 20)
            image = np.asarray(fig.axes[-1].images[0].get_array())
            self.assertTrue(np.all(image == 0))
            self.assertFalse(receipt['uniform_roi_fill'])
        finally:
            plt.close(fig)

    def test_main_merge_retains_pomc_variation_over_a_bright_base(self):
        # A bright base formerly clipped moderately bright POMC pixels to 1.
        dapi = self.channel.copy()
        cfos = self.channel.copy()
        npy = np.zeros_like(self.channel)
        got = renderer.composite_with_roi_masked_pomc_signal_and_npy_gfp(dapi, cfos, self.channel, npy, self.labels)
        red = got[3, 140:210, 0]
        self.assertTrue(np.all(red < 1))
        self.assertTrue(np.all(np.diff(red) > 0))
        base = np.clip(renderer.single_channel_view(dapi, 'DAPI') + renderer.single_channel_view(cfos, 'c-FOS'), 0, 1)
        np.testing.assert_allclose(got[self.labels == 0], base[self.labels == 0], atol=1e-15)

    def test_mismatched_microscopy_and_masks_fail(self):
        with self.assertRaises(ValueError):
            renderer.roi_microscopy_view(self.channel, self.labels[:, :-1], 'POMC')


if __name__ == '__main__':
    unittest.main()
