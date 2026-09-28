#!/usr/bin/env python3
"""Freeze the genuine 203 choices, 300 review targets and classifier inputs for GitHub."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile

PAPER = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def build(output):
    if output.exists():
        raise ValueError('Choose an absent output directory; an existing freeze is never replaced')
    source = PAPER / 'analyses/Fig5/results/soma_refined_20260907_v1'
    review = PAPER / 'Fig5/hil_review/microglia_masks_20260907_v2'
    candidate = review / 'retrained_bncal_20260907_v3'
    with (source / 'microglia_annotation_template.csv').open() as stream:
        cells = list(csv.DictReader(stream))
    with (review / 'microglia_reviewed_labels.csv').open() as stream:
        labels = list(csv.DictReader(stream))
    if len(cells) != 14415 or len(labels) != 203:
        raise ValueError('Expected the frozen 14,415 source cells and 203 genuine human choices')
    if (review / 'microglia_reviewed_labels.csv').read_bytes() != (review / 'labels_used_for_training.csv').read_bytes():
        raise ValueError('Current choices differ from the accepted training snapshot')
    inputs = []
    for name in ['per_microglia_cell_measurements.csv', 'microglia_annotation_template.csv', 'mask_refinement_report.json']:
        inputs.append((source / name, 'source/' + name))
    for row in cells:
        for kind in ['microglia_patches', 'soma_patches']:
            rel = f"{kind}/{row['sample']}/{row['section']}/{row['cell_uid']}.png"
            inputs.append((source / rel, 'source/' + rel))
    for name in ['review_batch.json', 'microglia_annotation_template.csv', 'microglia_reviewed_labels.csv', 'labels_used_for_training.csv']:
        inputs.append((review / name, 'review/' + name))
    for name in ['patches', 'masks']:
        inputs.extend((p, 'review/' + p.relative_to(review).as_posix()) for p in sorted((review / name).glob('*.png')))
    for p in sorted(candidate.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            inputs.append((p, 'candidate/' + p.relative_to(candidate).as_posix()))
    if any(not p.is_file() or p.is_symlink() for p, _ in inputs):
        raise ValueError('Missing or linked frozen input')
    groups, current, size = [], [], 0
    for item in inputs:
        n = item[0].stat().st_size
        if current and size + n > 8 * 1024 * 1024:
            groups.append(current); current = []; size = 0
        current.append(item); size += n
    if current:
        groups.append(current)
    output.mkdir(parents=True)
    archives = []
    for index, group in enumerate(groups, 1):
        path = output / f'microglial-review-{index:03d}.zip'
        members = []
        with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for src, name in group:
                data = src.read_bytes()
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, data)
                members.append({'path': name, 'size': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
        if path.stat().st_size > 10 * 1024 * 1024:
            raise ValueError(f'GitHub shard unexpectedly large: {path}')
        archives.append({'name': path.name, 'size': path.stat().st_size, 'sha256': sha(path), 'members': members})
        print(path.name, path.stat().st_size, 'bytes', flush=True)
    result = {'schema': 'apotome-microglial-choices-v1', 'date': '2026-09-07',
              'human_choices': 203, 'review_targets': 300, 'source_cells': 14415, 'animals': 15,
              'human_choices_sha256': sha(review / 'microglia_reviewed_labels.csv'),
              'original_review_batch_sha256': sha(review / 'review_batch.json'),
              'class_counts': {state: sum(r['state'] == state for r in labels) for state in ['Ramified','Rod-like','Activated','Amoeboid']},
              'status': 'Genuine saved human choices and their matching image/mask bytes; existing candidate, not a replacement of published Figure 5.',
              'archives': archives}
    result['human_choices_csv'] = 'human_choices_203.csv'
    (output / result['human_choices_csv']).write_bytes((review / 'microglia_reviewed_labels.csv').read_bytes())
    (output / 'manifest.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    (output / 'README.txt').write_text('''FROZEN MICROGLIAL CHOICES AND CLASSIFIER INPUTS\n\n203 genuine saved human classifications, 300 review crops with cell/soma/DAPI\nmasks, 14,415 original CNN patches and soma masks, source measurement/proposal\ntables, and the saved human-trained CNN candidate. No new microscopy was acquired.\nAll files are original bytes; the manifest binds each ZIP and member by SHA-256.\n\nUse Fig5/12_microglial_choices.py to verify/reuse these choices or initialize a\nseparate blank review with --do_microglial_choices. Portable initialization\nrewrites only working-copy paths, with a receipt, and preserves the frozen labels.\nThe archived candidate has 59.5% balanced accuracy on its 40-cell development\nvalidation split. Figure S6 separately reports the 203-cell morphology/DINOv2\ncomparison; these validation designs must not be conflated.\n\nCode: MIT. Original data, microscopy crops, masks and figures: CC BY 4.0.\n''')
    print(f'PASS: {len(inputs)} frozen files; {sum(a["size"] for a in archives)/2**20:.1f} MiB archives', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=PAPER / 'Fig5/microglial_review_data')
    build(parser.parse_args().output_dir)
