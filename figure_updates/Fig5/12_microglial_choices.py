#!/usr/bin/env python3
"""Reuse the 203 saved choices, or explicitly run a new 300-cell review/training session."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile
import zipfile

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/06_retrain_microglia_classes.py').is_file())
DEFAULT_BUNDLE = PAPER / 'Fig5/microglial_review_data'
STATES = ('Ramified', 'Rod-like', 'Activated', 'Amoeboid')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe(root, name):
    rel = PurePosixPath(name)
    if not name or rel.is_absolute() or rel.as_posix() != name or '..' in rel.parts or '\\' in name:
        raise ValueError(f'Unsafe bundle path: {name!r}')
    path = root
    for part in rel.parts:
        path = path / part
        if path.is_symlink():
            raise ValueError(f'Linked bundle path: {path}')
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Bundle path escapes its root')
    return path


def read_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows, fields=None):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def verify_bundle(bundle):
    manifest_path = safe(bundle, 'manifest.json')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema') != 'apotome-microglial-choices-v1' or manifest.get('human_choices') != 203 or manifest.get('review_targets') != 300:
        raise ValueError('Expected the frozen 203-choice, 300-target review bundle')
    if sha(safe(bundle, manifest['human_choices_csv'])) != manifest['human_choices_sha256']:
        raise ValueError('Saved human choices CSV differs from the frozen reference')
    seen = set()
    for part in manifest['archives']:
        archive = safe(bundle, part['name'])
        if archive.stat().st_size != part['size'] or sha(archive) != part['sha256']:
            raise ValueError(f'Frozen review archive changed: {archive.name}')
        for row in part['members']:
            safe(bundle, row['path'])
            if row['path'] in seen or PurePosixPath(row['path']).parts[0] not in {'source', 'review', 'candidate'}:
                raise ValueError('Duplicate or unexpected review bundle member')
            seen.add(row['path'])
    return manifest


def unpack(bundle, cache, manifest):
    if cache.is_symlink():
        raise ValueError('Cache cannot be a symlink')
    cache.mkdir(parents=True, exist_ok=True)
    frozen = cache / 'frozen'
    signature = sha(bundle / 'manifest.json')
    if frozen.exists():
        if frozen.is_symlink() or json.loads((frozen / 'IMPORT_RECEIPT.json').read_text())['manifest_sha256'] != signature:
            raise ValueError('Existing review cache belongs to a different bundle')
        for part in manifest['archives']:
            for row in part['members']:
                path = safe(frozen, row['path'])
                if not path.is_file() or path.stat().st_size != row['size'] or sha(path) != row['sha256']:
                    raise ValueError(f'Frozen cached review input changed: {row["path"]}')
        return frozen
    with tempfile.TemporaryDirectory(prefix='.import-', dir=cache) as temporary:
        stage = Path(temporary)
        for part in manifest['archives']:
            with zipfile.ZipFile(bundle / part['name']) as archive:
                expected = {row['path']: row for row in part['members']}
                if len(archive.infolist()) != len(expected) or set(archive.namelist()) != set(expected):
                    raise ValueError('ZIP members differ from the reviewed manifest')
                for member in archive.infolist():
                    row = expected[member.filename]
                    if (member.external_attr >> 16) & 0o170000 == 0o120000 or member.file_size != row['size']:
                        raise ValueError('Linked or wrong-sized ZIP member')
                    data = archive.read(member)
                    if hashlib.sha256(data).hexdigest() != row['sha256']:
                        raise ValueError(f'ZIP member hash differs: {member.filename}')
                    target = safe(stage, member.filename)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(data)
        (stage / 'IMPORT_RECEIPT.json').write_text(json.dumps({'manifest_sha256': signature, 'status': 'verified'}, indent=2) + '\n')
        # An existing destination is never replaced, including concurrent imports.
        if frozen.exists():
            raise ValueError('Another process initialized this cache; retry')
        stage.rename(frozen)
    return frozen


def initialize_review(frozen, root, fresh=False):
    if root.exists() or root.is_symlink():
        raise ValueError('Choose an absent working review directory')
    root.mkdir(parents=True)
    for name in ['patches', 'masks']:
        shutil.copytree(frozen / 'review' / name, root / name)
    for name in ['mask_history', 'mask_edits']:
        (root / name).mkdir()
    rows = read_rows(frozen / 'review/microglia_annotation_template.csv')
    for row in rows:
        row['patch_path'] = str(root / 'patches' / (row['cell_uid'] + '.png'))
        row['state'] = ''
    write_rows(root / 'microglia_annotation_template.csv', rows)
    meta = json.loads((frozen / 'review/review_batch.json').read_text())
    meta['source_analysis'] = str(frozen / 'source')
    meta['portable_origin_sha256'] = sha(frozen / 'review/review_batch.json')
    meta['source_measurements_sha256'] = sha(frozen / 'source/per_microglia_cell_measurements.csv')
    meta['source_template_sha256'] = sha(frozen / 'source/microglia_annotation_template.csv')
    if fresh:
        meta['previous_labels_preserved'] = 0
        meta['previous_targets_in_batch'] = 0
        meta['sampling'] = 'The same frozen 300 targets, independently reviewed with no prefilled human choices.'
    for row in meta['selected_cells']:
        template = next(r for r in rows if r['cell_uid'] == row['cell_uid'])
        row['source_patch'] = str(frozen / 'source/microglia_patches' / template['sample'] / template['section'] / (row['cell_uid'] + '.png'))
        if fresh:
            row['previous_state'] = ''
            row['previous_cell_uid'] = ''
    (root / 'review_batch.json').write_text(json.dumps(meta, indent=2) + '\n')
    labels = frozen / 'review/microglia_reviewed_labels.csv'
    if fresh:
        with labels.open(newline='') as stream:
            fields = next(csv.reader(stream))
        write_rows(root / 'microglia_reviewed_labels.csv', [], fields)
    else:
        shutil.copyfile(labels, root / labels.name)
        shutil.copyfile(labels, root / 'labels_used_for_training.csv')
    return root


def check_saved_choices(frozen, manifest):
    labels = frozen / 'review/microglia_reviewed_labels.csv'
    if sha(labels) != manifest['human_choices_sha256']:
        raise ValueError('The current choices differ from the frozen human labels')
    rows = read_rows(labels)
    batch = json.loads((frozen / 'review/review_batch.json').read_text())
    audit = {r['cell_uid']: r for r in batch['selected_cells']}
    if len(rows) != 203 or len({r['cell_uid'] for r in rows}) != 203 or len(audit) != 300:
        raise ValueError('The saved review count or identities differ')
    for row in rows:
        if row['state'] not in STATES or row['cell_uid'] not in audit or not all(row.get(k, '').strip() for k in ['reviewer', 'reviewed_at', 'session']):
            raise ValueError('Invalid or uncommitted human choice')
        if row['patch_sha256'] != audit[row['cell_uid']]['patch_sha256'] or sha(frozen / 'review/patches' / (row['cell_uid'] + '.png')) != row['patch_sha256']:
            raise ValueError('Human choice refers to a changed cell crop')
        for kind in ['cell', 'soma']:
            key = kind + '_mask_sha256'
            if row[key] != audit[row['cell_uid']][key] or sha(frozen / 'review/masks' / (row['cell_uid'] + '_' + kind + '.png')) != row[key]:
                raise ValueError('Human choice refers to a changed mask')
    counts = {state: sum(r['state'] == state for r in rows) for state in STATES}
    if counts != manifest['class_counts'] or sha(frozen / 'candidate/human_labels_used.csv') != sha(labels):
        raise ValueError('Saved candidate was not trained with these choices')
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle-dir', type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument('--cache-dir', type=Path, default=PAPER / 'analyses/Fig5/microglial_choices_cache_20260907')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--do_microglial_choices', '--do-microglial-choices', action='store_true', help='Optional: start a blank 300-cell human review and retrain; default reuses the current 203 choices and saved candidate.')
    parser.add_argument('--review-dir', type=Path, help='Fresh or resumable independent review directory (only with --do_microglial_choices).')
    parser.add_argument('--port', type=int, default=0)
    parser.add_argument('--prepare-only', action='store_true', help='Prepare the selected mode without opening a server or training.')
    parser.add_argument('--validate-only', action='store_true', help='Verify the frozen review bundle without writing or opening a server.')
    args = parser.parse_args()
    if args.review_dir and not args.do_microglial_choices:
        parser.error('--review-dir is only used with --do_microglial_choices')
    bundle = args.bundle_dir.resolve()
    manifest = verify_bundle(bundle)
    if args.validate_only:
        print('PASS: frozen 203 human choices, 300 review targets and 14,415 classifier inputs verified.', flush=True)
        return
    if args.output_dir is None:
        parser.error('--output-dir is required except with --validate-only')
    output = args.output_dir.resolve()
    if output.exists():
        parser.error('Use a fresh output directory; existing choices/results are preserved')
    cache = args.cache_dir.resolve()
    frozen = unpack(bundle, cache, manifest)
    counts = check_saved_choices(frozen, manifest)
    if args.do_microglial_choices:
        if args.review_dir:
            review = args.review_dir.resolve()
        else:
            sessions = cache / 'sessions'; sessions.mkdir(exist_ok=True)
            review = Path(tempfile.mkdtemp(prefix='review_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_', dir=sessions)) / 'choices'
        if not review.exists():
            initialize_review(frozen, review, fresh=True)
        elif review.is_symlink() or json.loads((review / 'review_batch.json').read_text()).get('source_analysis') != str(frozen / 'source'):
            raise ValueError('Existing review does not belong to this source bundle')
        candidate = review / 'retrained'
        mode = 'new_human_choices'
        print(f'MICROGLIAL_REVIEW_DIRECTORY: {review}', flush=True)
        if not args.prepare_only:
            subprocess.run([sys.executable, str(PAPER / 'Fig5/08_review_microglia_masks.py'),
                            '--analysis-dir', str(frozen / 'source'), '--review-dir', str(review),
                            '--port', str(args.port), '--training-python', sys.executable, '--exit-after-training'], check=True)
            state = json.loads((review / 'training_status.json').read_text())
            if state['status'] != 'complete' or state.get('returncode') != 0:
                raise ValueError('The new review has not completed training')
    else:
        review = cache / 'current_203_choices'
        if not review.exists():
            initialize_review(frozen, review)
        if sha(review / 'microglia_reviewed_labels.csv') != manifest['human_choices_sha256']:
            raise ValueError('Working default choices changed; use a separate independent review')
        candidate = frozen / 'candidate'
        mode = 'reuse_current_203_choices'
        print('MICROGLIAL_CHOICES: Reusing all 203 saved human choices and their trained candidate; no review prompt.', flush=True)
    output.mkdir(parents=True)
    labels = review / 'microglia_reviewed_labels.csv'
    shutil.copyfile(labels, output / 'human_microglial_choices.csv')
    if not args.prepare_only:
        shutil.copytree(candidate, output / 'candidate')
    used = read_rows(labels)
    receipt = {'mode': mode, 'labels': len(used), 'human_choices_sha256': sha(labels),
               'class_counts': {state: sum(r['state'] == state for r in used) for state in STATES},
               'review_targets': 300, 'source_cells': 14415, 'review_dir': str(review),
               'bundle_manifest_sha256': sha(bundle / 'manifest.json'), 'prepared_only': args.prepare_only,
               'candidate_saved': not args.prepare_only,
               'publication_status': 'Candidate/reference branch. Published Figure 5 and the frozen S6 comparison retain their established inputs.'}
    (output / 'MICROGLIAL_CHOICES_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    main()
