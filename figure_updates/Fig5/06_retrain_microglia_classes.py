#!/usr/bin/env python3
"""Retrain the existing Figure 5 CNN using only committed, hash-verified human labels.

Reuses the measured cells and their exact image patches. No segmentation or
anatomy changes are made. Outputs are versioned beside the review batch.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import shutil
import sys

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/01_analyze_gfap_iba1_microglia.py').is_file())
STATES = ('Ramified', 'Rod-like', 'Activated', 'Amoeboid')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def load_source(review_dir):
    import pandas as pd
    meta = json.loads((review_dir / 'review_batch.json').read_text())
    analysis = Path(meta['source_analysis'])
    measurement = analysis / 'per_microglia_cell_measurements.csv'
    template = analysis / 'microglia_annotation_template.csv'
    if sha256(measurement) != meta['source_measurements_sha256']:
        raise ValueError('Source measurements changed after this review batch was prepared.')
    if sha256(template) != meta['source_template_sha256']:
        raise ValueError('Source annotation template changed after batch preparation.')
    cells = pd.read_csv(measurement)
    proposals = pd.read_csv(template).set_index('cell_uid')
    if cells['cell_uid'].duplicated().any() or proposals.index.duplicated().any():
        raise ValueError('Cell identifiers must be unique.')
    if set(cells.cell_uid) != set(proposals.index):
        raise ValueError('Measurement and annotation cell identifiers differ.')
    cells['patch_path'] = [str(analysis / 'microglia_patches' / str(r.sample) / str(r.section) / (r.cell_uid + '.png')) for r in cells.itertuples()]
    # Relocated review bundles retain original CSV bytes for provenance. Resolve
    # reusable mask paths by stable cell identity, just as for CNN patches.
    for column, directory in [('soma_patch_path', 'soma_patches'), ('nucleus_patch_path', 'nucleus_patches')]:
        if column in cells:
            cells[column] = [str(analysis / directory / str(r.sample) / str(r.section) / (r.cell_uid + '.png')) for r in cells.itertuples()]
    missing = [p for p in cells.patch_path if not Path(p).is_file()]
    if missing:
        raise ValueError(f'{len(missing)} source patches are missing; first: {missing[0]}')
    if not cells['cell_qc_pass'].fillna(False).eq(True).all():
        raise ValueError('Source contains objects that failed biological-cell QC.')
    cell_lookup = cells.set_index('cell_uid')
    for row in meta['selected_cells']:
        uid = row['cell_uid']
        if uid not in cell_lookup.index:
            raise ValueError(f'Reviewed cell is absent from source: {uid}')
        source_patch = Path(cell_lookup.loc[uid, 'patch_path'])
        reviewed_patch = review_dir / 'patches' / (uid + '.png')
        if sha256(source_patch) != row.get('original_patch_sha256', row['patch_sha256']):
            raise ValueError(f'Source patch changed after review preparation: {source_patch}')
        if sha256(reviewed_patch) != row['patch_sha256']:
            raise ValueError(f'Reviewed patch changed: {reviewed_patch}')
        if meta.get('mask_review'):
            for kind in ['cell', 'soma']:
                mask_path = review_dir / 'masks' / f'{uid}_{kind}.png'
                if sha256(mask_path) != row[f'{kind}_mask_sha256']:
                    raise ValueError(f'Reviewed {kind} mask changed: {uid}')
            # The CNN must use the exact corrected mask channels the human saw.
            cells.loc[cells.cell_uid == uid, 'patch_path'] = str(reviewed_patch)
    return meta, cells, proposals.loc[cells.cell_uid]


def validate_labels(path, meta, cells, minimum):
    import pandas as pd
    from sklearn.model_selection import GroupShuffleSplit
    labels = pd.read_csv(path, keep_default_na=False)
    required = {'cell_uid', 'state', 'reviewer', 'reviewed_at', 'session', 'patch_sha256'}
    if not required.issubset(labels.columns):
        raise ValueError(f'Human review fields are missing: {required - set(labels.columns)}')
    if labels.cell_uid.duplicated().any():
        raise ValueError('Duplicate human label IDs.')
    audit = {r['cell_uid']:r for r in meta['selected_cells']}
    for row in labels.itertuples():
        if row.cell_uid not in audit or row.state not in STATES:
            raise ValueError(f'Invalid reviewed cell/state: {row.cell_uid}, {row.state}')
        if not all(str(getattr(row, k)).strip() for k in ('reviewer', 'reviewed_at', 'session')):
            raise ValueError(f'Uncommitted label: {row.cell_uid}')
        if row.patch_sha256 != audit[row.cell_uid]['patch_sha256']:
            raise ValueError(f'Label was committed for a different patch: {row.cell_uid}')
        if meta.get('mask_review'):
            if str(getattr(row, 'mask_reviewed', '')).lower() not in {'true', '1'}:
                raise ValueError(f'Cell mask has not been accepted with this classification: {row.cell_uid}')
            for kind in ['cell', 'soma']:
                if getattr(row, f'{kind}_mask_sha256', '') != audit[row.cell_uid][f'{kind}_mask_sha256']:
                    raise ValueError(f'Classification refers to an earlier {kind} mask: {row.cell_uid}')
    counts = Counter(labels.state)
    if len(labels) < minimum or any(counts[s] < 12 for s in STATES):
        raise ValueError(f'Need {minimum} committed labels and >=12 in each class; received {len(labels)} / {dict(counts)}')
    lookup = cells.set_index('cell_uid')
    lab = labels.set_index('cell_uid').loc[cells.loc[cells.cell_uid.isin(labels.cell_uid), 'cell_uid']]
    groups = lookup.loc[lab.index, 'animal_key'].astype(str).to_numpy()
    if len(set(groups)) < 3:
        raise ValueError('Need reviewed cells from at least three animals.')
    # Select solely by class coverage, before training or evaluating a model.
    # No accuracy-based split selection and no random cell-level split.
    for seed in range(20260620, 20260720):
        splitter = GroupShuffleSplit(n_splits=1, test_size=.2, random_state=seed)
        tr, va = next(splitter.split(lab, lab.state, groups))
        if set(lab.iloc[tr].state) == set(STATES) and set(lab.iloc[va].state) == set(STATES):
            return labels, lab, groups, tr, va, seed
    raise ValueError('Could not hold out whole animals with all four classes represented in both sets. Add class labels across more animals.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--review-dir', type=Path, required=True)
    ap.add_argument('--labels-csv', type=Path)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--min-labels', type=int, default=200)
    ap.add_argument('--epochs', type=int, default=100)
    ap.add_argument('--check-inputs', action='store_true')
    ap.add_argument('--validate-labels-only', action='store_true')
    ap.add_argument('--calibrate-batchnorm', action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args()
    if args.epochs < 1:
        ap.error('--epochs must be positive.')
    review = args.review_dir.resolve()
    meta, cells, proposals = load_source(review)
    print(f"Verified {len(cells)} QC-accepted source cells and {meta['cells']} frozen review patches.", flush=True)
    if args.check_inputs:
        return
    if args.labels_csv is None:
        ap.error('--labels-csv is required for retraining.')
    labels, lab, groups, tr, va, seed = validate_labels(args.labels_csv, meta, cells, args.min_labels)
    if args.validate_labels_only:
        print(f'Labels valid: {len(labels)}; training animals={sorted(set(groups[tr]))}; validation animals={sorted(set(groups[va]))}')
        return
    if args.output_dir is None:
        ap.error('--output-dir is required for retraining.')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    logger = logging.getLogger('fig5_supervised_retraining')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    import torch
    torch.set_num_threads(4)
    import pandas as pd
    import numpy as np
    from sklearn.metrics import balanced_accuracy_score, confusion_matrix, classification_report
    source = PAPER / 'Fig5/01_analyze_gfap_iba1_microglia.py'
    spec = importlib.util.spec_from_file_location('fig5_cnn_analyzer', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.import_science_stack()
    config = module.build_argparser_gfap_iba1().parse_args([])
    config.microglia_classifier_mode = 'supervised'
    config.microglia_labels_csv = str(args.labels_csv.resolve())
    config.cnn_min_training_cells = args.min_labels
    config.cnn_seed = seed
    config.cnn_device = 'cpu'
    config.cnn_model_in = None
    config.cnn_model_out = None
    config.cnn_calibrate_batchnorm = args.calibrate_batchnorm
    config.cnn_epochs = args.epochs
    (out / 'source_code').mkdir()
    for script in (source, Path(__file__).resolve()):
        shutil.copy2(script, out / 'source_code' / script.name)
    manual = dict(zip(labels.cell_uid, labels.state))
    prediction, training_meta = module._gi_train_and_predict_cnn(
        cells, proposals.pseudo_state.to_numpy(), proposals.pseudo_confidence.to_numpy(),
        manual, config, out, logger)
    if not training_meta['cnn_used'] or training_meta['training_candidates_n'] != len(labels):
        raise RuntimeError('Supervised CNN did not train using exactly the reviewed labels.')
    prediction.to_csv(out / 'per_microglia_classifier_predictions.csv', index=False)
    # The analyzer restores the best checkpoint before inference. Evaluate that
    # restored model, rather than retaining confusion counts from its last epoch.
    validation = lab.iloc[va]
    predicted = prediction.set_index('cell_uid').loc[validation.index, 'microglia_state']
    cm = confusion_matrix(validation.state, predicted, labels=list(STATES))
    pd.DataFrame(cm, index=STATES, columns=STATES).to_csv(out / 'microglia_cnn_validation_confusion_matrix.csv')
    split = lab[['state']].copy()
    split['animal_key'] = groups
    split['split'] = 'training'
    split.loc[validation.index, 'split'] = 'validation'
    split.to_csv(out / 'reviewed_label_split.csv')
    report = {**training_meta, 'source_analysis': meta['source_analysis'],
              'source_measurements_sha256': meta['source_measurements_sha256'],
              'analyzer_sha256': sha256(source), 'review_batch_sha256': sha256(review / 'review_batch.json'),
              'reviewed_labels_sha256': sha256(args.labels_csv),
              'training_script_sha256': sha256(Path(__file__).resolve()),
              'reviewed_masks_n': len(labels) if meta.get('mask_review') else 0,
              'training_configuration': {key: getattr(config, key) for key in ('cnn_epochs', 'cnn_patience', 'cnn_batch_size', 'cnn_learning_rate', 'cnn_weight_decay', 'cnn_calibrate_batchnorm')},
              'software': {'python': sys.version, 'torch': torch.__version__, 'numpy': np.__version__},
              'reviewed_mask_edits_n': sum(int(r.get('mask_revision', 0)) > 0 for r in meta['selected_cells']),
              'training_mode': 'supervised; only human labels; existing CNN architecture trained from scratch',
              'seed': seed, 'training_cells': len(tr), 'validation_cells': len(va),
              'fit_animals': sorted(set(groups[tr])), 'validation_animals': sorted(set(groups[va])),
              'training_class_counts': dict(Counter(lab.iloc[tr].state)),
              'validation_class_counts': dict(Counter(validation.state)),
              'validation_balanced_accuracy': float(balanced_accuracy_score(validation.state, predicted)),
              'validation_classification_report': classification_report(validation.state, predicted, labels=list(STATES), output_dict=True, zero_division=0),
              'interpretation': 'Animal-held-out validation on a proposed-class-balanced review sample; not an unbiased population prevalence estimate. Validation guides early stopping, so this is not a separate final test set.',
              'publication_status': 'Retrained candidate; performance and figure outputs require scientific assessment.'}
    (out / 'microglia_classifier_metadata.json').write_text(json.dumps(report, indent=2) + '\n')
    shutil.copy2(args.labels_csv, out / 'human_labels_used.csv')
    completion = out / 'retraining_report.json.tmp'
    completion.write_text(json.dumps(report, indent=2) + '\n')
    os.replace(completion, out / 'retraining_report.json')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
