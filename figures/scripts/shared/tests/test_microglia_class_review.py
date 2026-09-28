"""Cell-review persistence and human-only training handoff, in disposable fixtures."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch, Mock
from urllib.request import Request, urlopen
from urllib.error import HTTPError

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/05_review_microglia_classes.py').is_file())
spec = importlib.util.spec_from_file_location('fig5_review_test_module', PAPER / 'Fig5/05_review_microglia_classes.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='fig5-review-test-')
        self.root = Path(self.tmp.name)
        self.source = self.root / 'source'
        self.source.mkdir()
        rows = []
        for animal in range(3):
            for state in m.base.STATES:
                for n in range(4):
                    uid = f'test_{animal}_{state}_{n}'
                    row = {'cell_uid': uid, 'sample': str(animal), 'animal_id': str(animal),
                           'section': 'S01', 'region': 'ARC', 'patch_path': '/deleted/old/patch.png',
                           'state': '', 'pseudo_state': state, 'pseudo_confidence': '.8'}
                    p = self.source / 'microglia_patches' / str(animal) / 'S01' / (uid + '.png')
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_bytes(b'isolated-test-patch-' + uid.encode())
                    rows.append(row)
        with (self.source / m.base.TEMPLATE_NAME).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (self.source / 'per_microglia_cell_measurements.csv').write_text('fixture_only\n')
        self.batch = self.root / 'batch'
        self.meta = m.prepare(self.source, self.batch, 48, 42)
        self.review = m.base.Review(self.batch, self.batch / m.base.DEFAULT_OUTPUT_NAME, 42)
        self.training = m.Training(self.review, '/test/python', 48)

    def tearDown(self):
        self.tmp.cleanup()

    def commit(self, uid):
        return self.review.save({'cell_uid': uid, 'state': self.review.by_uid[uid]['pseudo_state'],
                                 'reviewer': 'UNIT_TEST_ONLY', 'session': 'DISPOSABLE_TEST_FIXTURE'})

    def test_balanced_batch_relocates_paths_without_prefilling_labels(self):
        self.assertEqual(self.meta['proposed_class_counts'], dict.fromkeys(m.base.STATES, 12))
        self.assertEqual(set(self.meta['animal_counts'].values()), {16})
        self.assertEqual(self.review.counts()['labelled'], 0)
        self.assertFalse(self.training.info()['ready'])
        self.assertFalse(self.review.output_csv.exists())
        for row in self.review.rows:
            self.assertTrue(Path(row['patch_path']).is_file())
            self.assertEqual(row['state'], '')

    def test_human_choice_persists_and_resumes(self):
        uid = self.review.order[0]
        self.commit(uid)
        reloaded = m.base.Review(self.batch, self.review.output_csv, 42)
        self.assertEqual(reloaded.counts()['labelled'], 1)
        self.assertEqual(reloaded.labels[uid]['reviewer'], 'UNIT_TEST_ONLY')
        self.assertEqual(reloaded.labels[uid]['patch_sha256'], m.base.sha256_file(Path(reloaded.by_uid[uid]['patch_path'])))
        with self.assertRaises(ValueError):
            self.review.save({'cell_uid': uid, 'state': 'unknown', 'reviewer': 'test', 'session': 'test'})
        with self.assertRaises(ValueError):
            self.review.save({'cell_uid': uid, 'state': 'Ramified'})

    def test_reopen_detects_changed_patch(self):
        row = self.meta['selected_cells'][0]
        (self.batch / 'patches' / (row['cell_uid'] + '.png')).write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'patch changed'):
            m.prepare(self.source, self.batch, 48, 42)

    def test_monitor_reads_the_versioned_run_report(self):
        original = self.batch / 'retrained'
        current = self.batch / 'retrained_second_run'
        original.mkdir()
        current.mkdir()
        (original / 'retraining_report.json').write_text(json.dumps({'run': 'old'}))
        (current / 'retraining_report.json').write_text(json.dumps({'run': 'current'}))
        self.training.state = {'status': 'running', 'output_dir': str(current)}
        self.training._resume_monitor()
        self.assertEqual(self.training.state['status'], 'complete')
        self.assertEqual(self.training.state['report']['run'], 'current')

    def test_training_cannot_start_from_proposals(self):
        with patch.object(m.subprocess, 'Popen') as launch:
            with self.assertRaisesRegex(ValueError, 'Still needed'):
                self.training.start()
            launch.assert_not_called()
        self.assertFalse((self.batch / 'labels_used_for_training.csv').exists())


    def test_training_split_keeps_animals_separate_and_all_classes_present(self):
        import pandas as pd
        spec = importlib.util.spec_from_file_location('fig5_training_test', PAPER / 'Fig5/06_retrain_microglia_classes.py')
        trainer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(trainer)
        for uid in self.review.order:
            self.commit(uid)
        cells = pd.DataFrame([{'cell_uid':r['cell_uid'], 'animal_key':r['animal_id']} for r in self.review.rows])
        labels, lab, groups, tr, va, seed = trainer.validate_labels(self.review.output_csv, self.meta, cells, 48)
        self.assertFalse(set(groups[tr]) & set(groups[va]))
        self.assertEqual(set(lab.iloc[tr].state), set(m.base.STATES))
        self.assertEqual(set(lab.iloc[va].state), set(m.base.STATES))
        labels.loc[0, 'patch_sha256'] = 'wrong'
        tampered = self.batch / 'tampered_test_labels.csv'
        labels.to_csv(tampered, index=False)
        with self.assertRaisesRegex(ValueError, 'different patch'):
            trainer.validate_labels(tampered, self.meta, cells, 48)

    def test_last_human_label_triggers_exactly_one_frozen_handoff(self):
        for uid in self.review.order[:-1]:
            self.commit(uid)
        server = m.base.ThreadingHTTPServer(('127.0.0.1', 0), m.Handler)
        server.review, server.training = self.review, self.training
        server.mutation_lock = threading.Lock()
        serving = threading.Thread(target=server.serve_forever, daemon=True)
        serving.start()
        url = f'http://127.0.0.1:{server.server_port}'
        stopped = threading.Event()
        process = Mock(pid=123456)
        process.wait.side_effect = lambda: (stopped.wait(3) or True) and 1
        try:
            with patch.object(m.subprocess, 'Popen', return_value=process) as launch:
                uid = self.review.order[-1]
                payload = {'cell_uid': uid, 'state': self.review.by_uid[uid]['pseudo_state'],
                           'reviewer': 'UNIT_TEST_ONLY', 'session': 'DISPOSABLE_TEST_FIXTURE'}
                req = Request(url + '/api/label', data=json.dumps(payload).encode(), headers={'Content-Type':'application/json'})
                with urlopen(req) as response:
                    self.assertEqual(response.status, 200)
                    json.load(response)
                # Serializing the next mutation ensures the automatic launch finished.
                with self.assertRaises(HTTPError) as context:
                    urlopen(Request(url + '/api/retrain', data=b'{}', headers={'Content-Type':'application/json'}))
                self.assertEqual(context.exception.code, 422)
                launch.assert_called_once()
                snapshot = self.batch / 'labels_used_for_training.csv'
                self.assertEqual(snapshot.read_bytes(), self.review.output_csv.read_bytes())
                with snapshot.open() as handle:
                    self.assertEqual(len(list(csv.DictReader(handle))), 48)
                argv = launch.call_args.args[0]
                self.assertIn('06_retrain_microglia_classes.py', argv[1])
                self.assertEqual(argv[argv.index('--labels-csv')+1], str(snapshot))
                with self.assertRaises(HTTPError) as context:
                    urlopen(req)
                self.assertEqual(context.exception.code, 409)
                stopped.set()
        finally:
            stopped.set()
            server.shutdown()
            server.server_close()
            serving.join()


if __name__ == '__main__':
    unittest.main()
