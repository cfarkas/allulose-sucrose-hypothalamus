"""Regressions for target preservation, image identity and animal separation."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from PIL import Image

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/09_improve_microglia_classifier.py').is_file())
spec = importlib.util.spec_from_file_location('fig5_transfer_tests', PAPER / 'Fig5/09_improve_microglia_classifier.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class TransferTests(unittest.TestCase):
    def test_crops_preserve_target_branches_and_do_not_add_annotation_colors(self):
        signal = np.tile(np.arange(96, dtype=np.uint8), (96, 1))
        cell = np.zeros((96, 96), bool)
        cell[1:9, 2:11] = True
        cell[4:6, 8:79] = True
        soma = np.zeros_like(cell)
        soma[2:7, 3:9] = True
        before = signal.copy()
        native = m.native_views(signal, cell, soma)
        self.assertEqual([int(mask.sum()) for _, mask in native], [int(cell.sum())] * 2)
        for rgb, weights in m.resized_views(signal, cell, soma):
            self.assertEqual(rgb.shape, (224, 224, 3))
            np.testing.assert_array_equal(rgb[..., 0], rgb[..., 1])
            np.testing.assert_array_equal(rgb[..., 1], rgb[..., 2])
            self.assertGreater(float(weights.sum()), 0)
        np.testing.assert_array_equal(signal, before)

    def test_mask_tampering_and_inconsistent_masks_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix='fig5-transfer-test-') as tmp:
            root = Path(tmp)
            cell = np.zeros((96, 96), bool)
            cell[42:54, 44:52] = True
            soma = np.zeros_like(cell)
            soma[45:50, 46:50] = True
            patch = np.zeros((96, 96, 3), np.uint8)
            patch[..., 0] = 80
            patch[..., 1] = cell * 255
            for name, array in [('patch', patch), ('cell', cell.astype(np.uint8)*255), ('soma', soma.astype(np.uint8)*255)]:
                Image.fromarray(array).save(root / f'{name}.png')
            row = {'cell_uid': 'DISPOSABLE_TEST', 'patch_path': str(root/'patch.png'), 'cell_mask_path': str(root/'cell.png'), 'soma_mask_path': str(root/'soma.png')}
            for name in ['patch', 'cell_mask', 'soma_mask']:
                row[name+'_sha256'] = m.sha256(row[name+'_path'])
            m.read_cell(row)
            changed = cell.copy()
            changed[1, 1] = True
            Image.fromarray(changed.astype(np.uint8)*255).save(root/'cell.png')
            with self.assertRaisesRegex(ValueError, 'mask changed'):
                m.read_cell(row)
            row['cell_mask_sha256'] = m.sha256(root/'cell.png')
            with self.assertRaisesRegex(ValueError, 'channel differ'):
                m.read_cell(row)

    def test_inner_validation_keeps_whole_animals_separate(self):
        groups = np.repeat([f'fixture_{i}' for i in range(15)], 4)
        y = np.tile(np.arange(4), 15)
        tested = []
        for train, validation in m.grouped_inner_splits(y, groups):
            self.assertFalse(set(groups[train]) & set(groups[validation]))
            self.assertEqual(set(y[train]), {0, 1, 2, 3})
            tested.extend(validation)
        self.assertEqual(sorted(tested), list(range(len(y))))

    def test_scaling_and_pca_fit_only_training_inputs(self):
        from threadpoolctl import threadpool_limits
        rng = np.random.default_rng(101)
        x = rng.normal(size=(60, 3 + m.EMBEDDING_WIDTH))
        x[48:] += 1000
        y = np.tile(np.arange(4), 12)
        with threadpool_limits(limits=2):
            model = m.make_model({'family': 'hybrid_logistic', 'C': 1.}, 3, 48).fit(x[:48], y)
        transforms = model.named_steps['features'].named_transformers_
        np.testing.assert_allclose(transforms['morphology'].mean_, x[:48, :3].mean(axis=0))
        np.testing.assert_allclose(transforms['dino'].named_steps['scale'].mean_, x[:48, 3:].mean(axis=0))
        self.assertLess(abs(transforms['dino'].named_steps['scale'].mean_).max(), 1.)
        self.assertEqual(transforms['dino'].named_steps['pca'].n_samples_, 48)
        self.assertEqual(m.probabilities(model, x[48:]).shape, (12, 4))


class TransferAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('fig5_transfer_audit_tests', PAPER / 'Fig5/10_audit_microglia_transfer.py')
        cls.audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.audit)

    def test_next_review_excludes_existing_labels_and_never_invents_human_labels(self):
        import pandas as pd
        candidates = pd.DataFrame([
            {'cell_uid': f'{animal}_{state}_{i}', 'animal_key': animal,
             'proposed_state': state, 'review_priority': float(i)}
            for animal in ['animal_a', 'animal_b'] for state in m.STATES for i in range(6)
        ])
        # These high-priority cells already have human labels and must stay out.
        existing = {'animal_a_Amoeboid_5', 'animal_b_Ramified_5'}
        original = candidates.copy(deep=True)
        queue = self.audit.select_review_queue(candidates, existing)
        self.assertEqual(len(queue), 40)
        self.assertEqual(set(queue.groupby('animal_key').size()), {20})
        self.assertFalse(set(queue.cell_uid) & existing)
        self.assertTrue(queue.human_state.eq('').all())
        self.assertTrue(queue.mask_review_status.eq('pending').all())
        pd.testing.assert_frame_equal(candidates, original)

    def test_validation_audit_rejects_same_animal_in_fit_and_test(self):
        import pandas as pd
        manifest = pd.DataFrame({'cell_uid': ['a1', 'a2', 'b1'], 'animal_key': ['a', 'a', 'b']})
        bad = [{'held_out_animal': 'a', 'training_uids': ['a2', 'b1'], 'test_uids': ['a1'], 'inner_folds': []}]
        with self.assertRaisesRegex(ValueError, 'Animal leakage'):
            self.audit.validate_splits(bad, manifest)


if __name__ == '__main__':
    unittest.main()
