#!/usr/bin/env python3
"""Audit Figure 5 classifier candidates and prepare a separate, unlabeled review queue."""
from __future__ import annotations
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from threadpoolctl import threadpool_limits

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/09_improve_microglia_classifier.py').is_file())
spec = importlib.util.spec_from_file_location('fig5_transfer', PAPER / 'Fig5/09_improve_microglia_classifier.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)


def validate_splits(design, manifest):
    uid_to_animal = manifest.set_index('cell_uid').animal_key.to_dict()
    universe = set(uid_to_animal)
    seen = []
    for fold in design:
        fit, test = set(fold['training_uids']), set(fold['test_uids'])
        if fit & test or fit | test != universe:
            raise ValueError('Outer fold does not partition all reviewed cells.')
        fit_animals = {uid_to_animal[u] for u in fit}
        test_animals = {uid_to_animal[u] for u in test}
        if test_animals != {fold['held_out_animal']} or fit_animals & test_animals:
            raise ValueError('Animal leakage in outer validation.')
        inner_seen = []
        for inner in fold['inner_folds']:
            a, b = set(inner['fit_uids']), set(inner['selection_uids'])
            if a & b or a | b != fit:
                raise ValueError('Inner fold does not partition the outer training cells.')
            if {uid_to_animal[u] for u in a} & {uid_to_animal[u] for u in b}:
                raise ValueError('Animal leakage in inner selection.')
            inner_seen.extend(b)
        if Counter(inner_seen) != Counter({u: 1 for u in fit}):
            raise ValueError('Each training cell must be selected out once in inner CV.')
        seen.extend(test)
    if Counter(seen) != Counter({u: 1 for u in universe}):
        raise ValueError('Each reviewed cell must be held out exactly once.')


def select_review_queue(candidates, previously_reviewed, per_animal=20):
    """Balance animals and proposed classes, then rank disagreement/uncertainty.

    This produces acquisition targets only. It never supplies a human class.
    """
    pool = candidates.loc[~candidates.cell_uid.isin(previously_reviewed)].copy()
    if pool.cell_uid.duplicated().any():
        raise ValueError('Duplicate cell identity in review candidates.')
    if not set(pool.proposed_state).issubset(m.STATES):
        raise ValueError('Unknown proposed class.')
    selected = []
    for animal, group in pool.groupby('animal_key', sort=True):
        group = group.sort_values(['review_priority', 'cell_uid'], ascending=[False, True])
        chosen = group.groupby('proposed_state', sort=False).head(per_animal // 4)
        if len(chosen) < per_animal:
            chosen = pd.concat([chosen, group.loc[~group.cell_uid.isin(chosen.cell_uid)].head(per_animal-len(chosen))])
        if len(chosen) != per_animal:
            raise ValueError(f'Insufficient new review targets for {animal}.')
        selected.append(chosen)
    queue = pd.concat(selected, ignore_index=True)
    queue = queue.sort_values(['review_priority', 'animal_key', 'cell_uid'], ascending=[False, True, True]).reset_index(drop=True)
    queue.insert(0, 'review_order', np.arange(1, len(queue)+1))
    queue['human_state'] = ''
    queue['reviewer'] = ''
    queue['mask_review_status'] = 'pending'
    if set(queue.cell_uid) & set(previously_reviewed):
        raise ValueError('An existing human label was selected for replacement.')
    return queue


def render_error_atlas(run, errors, manifest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    lookup = manifest.set_index('cell_uid')
    # At most three per true class keeps this diagnostic atlas readable.
    chosen = errors.groupby('human_state', sort=False).head(3).reset_index(drop=True)
    path = run / 'held_out_error_examples.pdf'
    with PdfPages(path) as pdf:
        for start in range(0, len(chosen), 6):
            sub = chosen.iloc[start:start+6]
            fig, axes = plt.subplots(len(sub), 2, figsize=(10, 2.25*len(sub)), squeeze=False)
            for i, row in enumerate(sub.itertuples(index=False)):
                record = lookup.loc[row.cell_uid].to_dict() | {'cell_uid': row.cell_uid}
                signal, cell, soma = m.read_cell(record)
                axes[i, 0].imshow(signal, cmap='gray', vmin=0, vmax=255)
                axes[i, 0].contour(cell, levels=[.5], colors=['#00dbe0'], linewidths=.6)
                axes[i, 0].contour(soma, levels=[.5], colors=['#ffcf40'], linewidths=.6)
                axes[i, 0].set_title(f'{row.cell_uid}\nHuman: {row.human_state}', fontsize=9)
                close, _ = m.native_views(signal, cell, soma)[0]
                axes[i, 1].imshow(close, cmap='gray', vmin=0, vmax=255)
                axes[i, 1].set_title(f'Held-out model: {row.predicted_state} ({row.model_score:.2f})\nExisting proposal: {row.proposed_state}', fontsize=9)
                for ax in axes[i]:
                    ax.axis('off')
            fig.suptitle('Figure 5 held-out errors: actual signal and cell/soma contours', fontsize=12)
            fig.text(.5, .018, 'Cyan: cell mask; yellow: soma. Human labels are retained. Model scores are uncalibrated.\nExamples chosen by descending model score within each human class; the whole animal was excluded from training.', ha='center', fontsize=8)
            fig.tight_layout(rect=[0, .05, 1, .97])
            pdf.savefig(fig)
            if start == 0:
                fig.savefig(run / 'held_out_error_examples.png', dpi=200)
            plt.close(fig)
    chosen.to_csv(run / 'held_out_error_examples.csv', index=False)


def audit(args):
    run = args.run_dir.resolve()
    protocol, reviewed, xr = m.load_matrix(run, reviewed_only=True)
    _, manifest, x = m.load_matrix(run, reviewed_only=False)
    report = json.loads((run / 'transfer_training_report.json').read_text())
    review_root = Path(protocol['review_dir'])
    if m.sha256(review_root / 'microglia_reviewed_labels.csv') != protocol['human_labels_sha256']:
        raise ValueError('Current human labels changed since this run was frozen.')
    source = Path(protocol['source_analysis'])
    if m.sha256(source / 'per_microglia_cell_measurements.csv') != protocol['source_measurements_sha256']:
        raise ValueError('Source measurements changed.')
    design_path = run / 'validation/validation_design.json'
    if m.sha256(design_path) != report['validation_design_sha256']:
        raise ValueError('Validation design changed after training.')
    design = json.loads(design_path.read_text())
    validate_splits(design, reviewed)
    labels = pd.read_csv(run / 'human_labels_used.csv').set_index('cell_uid').state
    y = np.asarray([m.STATES.index(s) for s in labels.loc[reviewed.cell_uid]], dtype=int)
    pcols = ['prob_' + s.lower().replace('-', '') for s in m.STATES]
    audited = {}
    for method, expected in report['nested_validation'].items():
        frame = pd.read_csv(run / f'validation/{method}_out_of_fold.csv')
        if not frame.cell_uid.equals(reviewed.cell_uid) or not frame.human_state.equals(reviewed.cell_uid.map(labels)):
            raise ValueError(f'Misaligned held-out predictions for {method}.')
        p = frame[pcols].to_numpy()
        np.testing.assert_allclose(p.sum(axis=1), 1., atol=1e-7)
        actual = m.scores(y, p)
        for metric in ('balanced_accuracy', 'accuracy'):
            np.testing.assert_allclose(actual[metric], expected[metric], atol=1e-12)
        stored_cm = pd.read_csv(run / f'validation/{method}_confusion.csv', index_col=0)
        np.testing.assert_array_equal(stored_cm.loc[list(m.STATES), list(m.STATES)], confusion_matrix(y, p.argmax(axis=1), labels=list(range(4))))
        audited[method] = actual['balanced_accuracy']
    # Independently refit every family for one predetermined outer animal.
    first = design[0]
    choices = json.loads((run / 'validation/model_selection_audit.json').read_text())
    ranked = sorted(choices[0]['inner_results'], key=lambda r: (-r['balanced_accuracy'], r['log_loss'], r['ordinal']))
    fit = np.flatnonzero(reviewed.cell_uid.isin(first['training_uids']))
    test = np.flatnonzero(reviewed.cell_uid.isin(first['test_uids']))
    with threadpool_limits(limits=4):
        for family in m.FAMILIES:
            config = next(r['configuration'] for r in ranked if r['configuration']['family'] == family)
            fresh = m.make_model(config, len(protocol['shape_features']), len(fit)).fit(xr[fit], y[fit])
            p = m.probabilities(fresh, xr[test])
            stored = pd.read_csv(run / f'validation/{family}_out_of_fold.csv').iloc[test][pcols].to_numpy()
            np.testing.assert_allclose(p, stored, atol=1e-10)
    predicted_dir = run / 'model_predictions'
    predicted_dir.mkdir(exist_ok=False)
    receipts, family_predictions = {}, {}
    for family in m.FAMILIES:
        model_path = run / f'models/{family}.joblib'
        bundle = joblib.load(model_path)
        if bundle['recipe_sha256'] != protocol['recipe_sha256'] or bundle['human_labels_sha256'] != protocol['human_labels_sha256']:
            raise ValueError('Saved model provenance mismatch.')
        if bundle['training_uids'] != reviewed.cell_uid.tolist():
            raise ValueError('Saved model training identities mismatch.')
        # Confirm the full-data scalers include exactly the reviewed inputs.
        transformers = bundle['model'].named_steps['features'].named_transformers_
        nshape = len(protocol['shape_features'])
        if 'morphology' in transformers:
            np.testing.assert_allclose(transformers['morphology'].mean_, xr[:, :nshape].mean(axis=0), atol=1e-12)
        if 'dino' in transformers:
            np.testing.assert_allclose(transformers['dino'].named_steps['scale'].mean_, xr[:, nshape:].mean(axis=0), atol=1e-12)
        with threadpool_limits(limits=4):
            p = m.probabilities(bundle['model'], x)
        if not np.isfinite(p).all():
            raise ValueError('Nonfinite candidate prediction.')
        np.testing.assert_allclose(p.sum(axis=1), 1., atol=1e-7)
        frame = manifest[['cell_uid']].copy()
        frame['predicted_state'] = [m.STATES[i] for i in p.argmax(axis=1)]
        frame['model_score'] = p.max(axis=1)
        for i, col in enumerate(pcols):
            frame[col] = p[:, i]
        target = predicted_dir / f'{family}.csv'
        frame.to_csv(target, index=False)
        family_predictions[family] = frame
        receipts[family] = {'model_sha256': m.sha256(model_path), 'predictions_sha256': m.sha256(target),
                            'cells_n': len(frame), 'class_counts': frame.predicted_state.value_counts().to_dict()}
    selected_family = report['selected_final_configuration']['family']
    original = pd.read_csv(run / 'per_microglia_classifier_predictions.csv')
    selected = family_predictions[selected_family]
    if not original.cell_uid.equals(selected.cell_uid):
        raise ValueError('Final prediction row order changed.')
    np.testing.assert_allclose(original[pcols], selected[pcols], atol=1e-12)
    # New labels remain blank; candidates are balanced across animals/proposed classes.
    proposals = pd.read_csv(source / 'microglia_annotation_template.csv').set_index('cell_uid').pseudo_state
    pool = manifest.copy()
    pool['proposed_state'] = pool.cell_uid.map(proposals)
    pool['candidate_state'] = selected.predicted_state
    pool['candidate_score'] = selected.model_score
    pool['dino_state'] = family_predictions['dino_logistic'].predicted_state
    pool['dino_score'] = family_predictions['dino_logistic'].model_score
    pool['proposal_disagreement'] = pool.proposed_state.ne(pool.candidate_state)
    pool['image_shape_disagreement'] = pool.dino_state.ne(pool.candidate_state)
    pool['review_priority'] = (pool.proposal_disagreement.astype(float) + .5 * pool.image_shape_disagreement.astype(float)
                               + 1. - pool.candidate_score)
    queue = select_review_queue(pool, set(labels.index))
    if len(queue) != 300 or queue.animal_key.nunique() != 15 or not queue.human_state.eq('').all():
        raise ValueError('Expected 300 unlabeled targets across 15 animals.')
    for row in queue.to_dict('records'):
        m.read_cell(row)  # Validate the exact image/mask identities offered for review.
    queue_path = run / 'next_review_queue_300.csv'
    queue.to_csv(queue_path, index=False)
    out_of_fold = pd.read_csv(run / 'validation/nested_selected_out_of_fold.csv')
    out_of_fold['proposed_state'] = out_of_fold.cell_uid.map(proposals)
    out_of_fold['model_score'] = out_of_fold[pcols].max(axis=1)
    errors = out_of_fold.loc[out_of_fold.human_state.ne(out_of_fold.predicted_state)].sort_values(['model_score', 'cell_uid'], ascending=[False, True])
    errors.to_csv(run / 'held_out_errors.csv', index=False)
    render_error_atlas(run, errors, manifest)
    current = report['nested_validation']['morphometric_proposals']['balanced_accuracy']
    candidate = report['nested_validation']['nested_selected']['balanced_accuracy']
    delta = report['nested_validation']['nested_selected']['paired_difference_vs_morphometry_95_interval']
    improved = candidate > current and delta[0] > 0
    assessment = {
        'status': 'audited_candidates', 'scientific_decision': 'candidate_merits_independent_review' if improved else 'retain_existing_figure_classifier',
        'reason': 'The animal-held-out comparison does not establish an improvement over existing morphology proposals.' if not improved else 'Positive paired animal-bootstrap difference; external validation remains necessary.',
        'publication_figures_modified': False, 'original_human_labels_modified': False,
        'nested_balanced_accuracy': audited, 'selected_pipeline_difference_vs_proposals': candidate-current,
        'paired_difference_95_interval': delta,
        'interpretation': 'Internal development validation on a proposed-class-balanced review sample. The historical 40-cell scores are not an independent test. Model scores are uncalibrated. The whole-animal bootstrap conditions on fixed OOF predictions.',
        'checks': {'outer_inner_animal_separation': True, 'every_reviewed_cell_held_out_once': True,
                   'all_metrics_recomputed': True, 'first_outer_animal_four_models_refitted_and_matched': first['held_out_animal'],
                   'all_final_models_use_only_reviewed_training_inputs': True, 'selected_predictions_reproduced': True,
                   'next_review_image_and_mask_hashes_verified': 300},
        'candidate_predictions': receipts, 'next_review_queue': {'path': str(queue_path), 'sha256': m.sha256(queue_path),
                   'cells': len(queue), 'animals': queue.animal_key.nunique(), 'human_labels_added': 0,
                   'per_animal': queue.animal_key.value_counts().to_dict(), 'per_proposed_state': queue.proposed_state.value_counts().to_dict()},
        'training_report_sha256': m.sha256(run / 'transfer_training_report.json'), 'audit_script_sha256': m.sha256(Path(__file__))}
    m.atomic_json(run / 'improvement_assessment.json', assessment)
    shutil.copy2(Path(__file__), run / 'source_code/10_audit_microglia_transfer.py')
    print(json.dumps(assessment, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    audit(parser.parse_args())


if __name__ == '__main__':
    main()
