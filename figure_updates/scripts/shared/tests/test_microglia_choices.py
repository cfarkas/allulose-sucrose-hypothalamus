"""Defaults preserve genuine saved choices; independent reviews start blank."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

PAPER = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('microglial_choices', PAPER / 'Fig5/12_microglial_choices.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class ChoiceInitializationTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.frozen = self.root / 'frozen'
        self.review = self.frozen / 'review'
        for directory in ['review/patches', 'review/masks', 'source']:
            (self.frozen / directory).mkdir(parents=True)
        (self.review / 'patches/testcell.png').write_bytes(b'synthetic test image fixture')
        (self.review / 'masks/testcell_soma.png').write_bytes(b'synthetic test mask fixture')
        (self.frozen / 'source/per_microglia_cell_measurements.csv').write_text('cell_uid\ntestcell\n')
        (self.frozen / 'source/microglia_annotation_template.csv').write_text('cell_uid\ntestcell\n')
        m.write_rows(self.review / 'microglia_annotation_template.csv', [{'cell_uid':'testcell','sample':'testsample','section':'S01','patch_path':'/old/location.png','state':'Ramified'}])
        (self.review / 'review_batch.json').write_text(json.dumps({'selected_cells':[{'cell_uid':'testcell','previous_state':'Ramified','previous_cell_uid':'oldcell'}], 'previous_labels_preserved':1}))
        self.saved = b'cell_uid,state,reviewer\ntestcell,Ramified,SYNTHETIC TEST FIXTURE\n'
        (self.review / 'microglia_reviewed_labels.csv').write_bytes(self.saved)

    def test_default_preserves_choice_bytes_and_rebases_working_paths(self):
        target = self.root / 'default'
        m.initialize_review(self.frozen, target)
        self.assertEqual((target / 'microglia_reviewed_labels.csv').read_bytes(), self.saved)
        self.assertEqual((target / 'labels_used_for_training.csv').read_bytes(), self.saved)
        rows = m.read_rows(target / 'microglia_annotation_template.csv')
        self.assertEqual(rows[0]['patch_path'], str(target / 'patches/testcell.png'))
        self.assertEqual(json.loads((target / 'review_batch.json').read_text())['source_analysis'], str(self.frozen / 'source'))
        self.assertEqual((self.review / 'microglia_reviewed_labels.csv').read_bytes(), self.saved)

    def test_optional_independent_review_has_zero_prefilled_choices(self):
        target = self.root / 'new'
        m.initialize_review(self.frozen, target, fresh=True)
        self.assertEqual(m.read_rows(target / 'microglia_reviewed_labels.csv'), [])
        self.assertFalse((target / 'labels_used_for_training.csv').exists())
        meta = json.loads((target / 'review_batch.json').read_text())
        self.assertEqual(meta['previous_labels_preserved'], 0)
        self.assertEqual(meta['selected_cells'][0]['previous_state'], '')
        self.assertEqual(meta['selected_cells'][0]['previous_cell_uid'], '')
        self.assertEqual((self.review / 'microglia_reviewed_labels.csv').read_bytes(), self.saved)

    def test_existing_review_is_not_erased_by_initialization(self):
        target = self.root / 'in_progress'; target.mkdir()
        (target / 'choices.csv').write_text('preserve this work')
        with self.assertRaisesRegex(ValueError, 'absent'):
            m.initialize_review(self.frozen, target, fresh=True)
        self.assertEqual((target / 'choices.csv').read_text(), 'preserve this work')

    def test_cache_paths_reject_traversal_and_links(self):
        for name in ['../outside', '/absolute', 'a/../../outside']:
            with self.assertRaises(ValueError): m.safe(self.root, name)
        (self.root / 'link').symlink_to(self.review, target_is_directory=True)
        with self.assertRaises(ValueError): m.safe(self.root, 'link/patches/testcell.png')


if __name__ == '__main__': unittest.main()
