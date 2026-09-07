"""Regression checks for training-only normalization on small human review sets."""
import importlib.util
import logging
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from PIL import Image
import torch
from sklearn.model_selection import GroupShuffleSplit

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/01_analyze_gfap_iba1_microglia.py').is_file())
spec = importlib.util.spec_from_file_location('microglia_calibration_test_analyzer', PAPER / 'Fig5/01_analyze_gfap_iba1_microglia.py')
analyzer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = analyzer
spec.loader.exec_module(analyzer)
analyzer.import_science_stack()


class CalibrationTests(unittest.TestCase):
    def test_validation_cells_do_not_enter_normalization_and_optimization_is_unchanged(self):
        torch.set_num_threads(2)
        with tempfile.TemporaryDirectory(prefix='fig5-bn-test-') as tmp:
            root = Path(tmp)
            args = analyzer.build_argparser_gfap_iba1().parse_args([])
            args.microglia_classifier_mode = 'supervised'
            args.cnn_min_training_cells = 12
            args.cnn_min_class_cells = 2
            args.cnn_device = 'cpu'
            args.cnn_epochs = 2
            args.cnn_batch_size = 8
            args.cnn_patch_size = 32
            states = list(analyzer.MICROGLIA_STATES)
            groups = np.repeat(['fixture_a', 'fixture_b', 'fixture_c'], 8)
            targets = np.tile(np.repeat(states, 2), 3)
            tr, va = next(GroupShuffleSplit(n_splits=1, test_size=.2, random_state=args.cnn_seed).split(groups, targets, groups))
            rows, images = [], []
            for i, (group, state) in enumerate(zip(groups, targets)):
                uid = f'FIXTURE_ONLY_{i:03d}'
                a = np.zeros((32, 32, 3), np.uint8)
                # A large validation-only brightness shift makes leakage visible
                # directly in the saved first-layer normalization statistics.
                a[..., 0] = 240 if i in va else 30 + i % 4
                a[10:22, 12:20, 1] = 255
                a[14:18, 12:20, 2] = 150
                path = root / f'{uid}.png'
                Image.fromarray(a).save(path)
                rows.append(dict(cell_uid=uid, animal_key=group, patch_path=str(path), cell_qc_pass=True))
                images.append(a)
            cells = pd.DataFrame(rows)
            manual = dict(zip(cells.cell_uid, targets))
            for calibrated in [False, True]:
                args.cnn_calibrate_batchnorm = calibrated
                output = root / ('calibrated' if calibrated else 'legacy')
                output.mkdir()
                _, meta = analyzer._gi_train_and_predict_cnn(cells, targets, np.ones(len(cells)), manual, args, output, logging.getLogger('fixture'))
            checkpoint = torch.load(root / 'calibrated/models/microglia_morphology_cnn.pt', map_location='cpu', weights_only=False)
            state = checkpoint['model_state']
            x = torch.from_numpy(np.stack(images).transpose(0, 3, 1, 2).astype(np.float32) / 255)
            activations = torch.nn.functional.conv2d(x, state['features.0.weight'], state['features.0.bias'], padding=1)
            expected = activations[tr].mean(dim=(0, 2, 3))
            leaked = activations.mean(dim=(0, 2, 3))
            torch.testing.assert_close(state['features.1.running_mean'], expected, atol=1e-6, rtol=1e-5)
            self.assertGreater(float(torch.linalg.vector_norm(expected - leaked)), .01)
            self.assertEqual(int(state['features.1.num_batches_tracked']), math.ceil(len(tr) / args.cnn_batch_size))
            self.assertEqual(meta['batchnorm_calibration'], 'training_cells_before_validation')
            self.assertIn(meta['best_epoch'], [1, 2])
            legacy = pd.read_csv(root / 'legacy/microglia_cnn_training_history.csv')
            calibrated = pd.read_csv(root / 'calibrated/microglia_cnn_training_history.csv')
            np.testing.assert_allclose(legacy.train_loss, calibrated.train_loss, atol=1e-8, rtol=0)


if __name__ == '__main__':
    unittest.main()
