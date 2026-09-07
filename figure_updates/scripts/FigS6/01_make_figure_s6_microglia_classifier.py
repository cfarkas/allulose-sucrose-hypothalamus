#!/usr/bin/env python3
"""Rebuild Supplementary Figure S6 from frozen, human-referenced predictions.

This renderer recomputes all plotted results; it never trains classifiers or
modifies Figure 5, microscopy, human reviews, or its own frozen inputs.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil

os.environ.setdefault('SOURCE_DATE_EPOCH', '0')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'scripts/setup').is_dir() and (p / 'README.txt').is_file())
DEFAULT_INPUT = PAPER / 'FigS6/raw_data'
STATES = ('Ramified', 'Rod-like', 'Activated', 'Amoeboid')
METHODS = ('morphometric_proposals', 'dino_logistic')
COLORS = ('#277b83', '#c66a2c')
SLUGS = ('A_balanced_accuracy', 'B_class_recall', 'C_morphometric_confusion', 'D_DINOv2_confusion')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def confusion(y, predicted):
    return np.bincount(4 * y + predicted, minlength=16).reshape(4, 4)


def balanced_accuracy(y, predicted):
    matrix = confusion(y, predicted)
    support = matrix.sum(axis=1)
    if np.any(support == 0):
        raise ValueError('Balanced accuracy requires all four reference classes.')
    return float(np.mean(np.diag(matrix) / support))


def validate_design(design, reference):
    mapping = reference.set_index('cell_uid').animal_key.to_dict()
    all_uids = set(mapping)
    outer_seen = []
    for fold in design:
        fit, test = set(fold['training_uids']), set(fold['test_uids'])
        if fit & test or fit | test != all_uids:
            raise ValueError('Invalid outer partition of reviewed cells.')
        fit_animals = {mapping[u] for u in fit}
        test_animals = {mapping[u] for u in test}
        if test_animals != {fold['held_out_animal']} or fit_animals & test_animals:
            raise ValueError('A held-out animal leaked into training.')
        inner_seen = []
        for inner in fold['inner_folds']:
            a, b = set(inner['fit_uids']), set(inner['selection_uids'])
            if a & b or a | b != fit:
                raise ValueError('Invalid inner partition of outer training cells.')
            if {mapping[u] for u in a} & {mapping[u] for u in b}:
                raise ValueError('An animal leaked into model selection.')
            inner_seen.extend(b)
        if Counter(inner_seen) != Counter({u: 1 for u in fit}):
            raise ValueError('Inner validation did not cover each fitting cell once.')
        outer_seen.extend(test)
    if Counter(outer_seen) != Counter({u: 1 for u in all_uids}):
        raise ValueError('Outer validation did not cover each reviewed cell once.')


def load_inputs(root):
    manifest_path = root / 'RAW_DATA_MANIFEST.csv'
    manifest = pd.read_csv(manifest_path)
    required = {'method_metadata.json', 'human_reference.csv', 'morphometric_predictions.csv',
                'dinov2_predictions.csv', 'validation_design.json'}
    if set(manifest.file) != required or manifest.file.duplicated().any():
        raise ValueError('Frozen comparison manifest is incomplete or duplicated.')
    for record in manifest.itertuples(index=False):
        path = root / record.file
        if path.stat().st_size != record.size_bytes or sha256(path) != record.sha256:
            raise ValueError(f'Frozen figure input changed: {record.file}')
    metadata = json.loads((root / 'method_metadata.json').read_text())
    reference = pd.read_csv(root / 'human_reference.csv', keep_default_na=False)
    if reference.cell_uid.duplicated().any() or reference.animal_key.eq('').any():
        raise ValueError('Human reference contains duplicate cells or missing animals.')
    if set(reference.human_state) != set(STATES) or metadata['states'] != list(STATES):
        raise ValueError('Human reference class set differs.')
    if len(reference) != metadata['reviewed_cells_n'] or reference.animal_key.nunique() != metadata['animals_n']:
        raise ValueError('Human reference sample size differs from the frozen protocol.')
    if reference.human_state.value_counts().to_dict() != metadata['human_class_counts']:
        raise ValueError('Human class counts differ from the frozen protocol.')
    design = json.loads((root / 'validation_design.json').read_text())
    validate_design(design, reference)
    frames = {}
    for method in METHODS:
        frame = pd.read_csv(root / metadata['methods'][method]['frozen_predictions'])
        if frame.cell_uid.duplicated().any() or set(frame.cell_uid) != set(reference.cell_uid):
            raise ValueError(f'{method} did not predict the same unique reviewed cells.')
        frame = frame.set_index('cell_uid').loc[reference.cell_uid].reset_index()
        for column in ('cell_uid', 'animal_key', 'human_state'):
            if not frame[column].equals(reference[column]):
                raise ValueError(f'{method}: reference {column} differs.')
        pcols = ['prob_' + s.lower().replace('-', '') for s in STATES]
        probabilities = frame[pcols].to_numpy(float)
        if not np.isfinite(probabilities).all() or np.any(probabilities < 0) or np.any(probabilities > 1):
            raise ValueError(f'{method}: invalid model scores.')
        np.testing.assert_allclose(probabilities.sum(axis=1), 1., atol=1e-7)
        inferred = [STATES[i] for i in probabilities.argmax(axis=1)]
        if frame.predicted_state.tolist() != inferred:
            raise ValueError(f'{method}: stored class disagrees with model score maximum.')
        frames[method] = frame
    return metadata, reference, frames, manifest


def calculate(metadata, reference, frames):
    y = np.asarray([STATES.index(s) for s in reference.human_state], dtype=int)
    groups = reference.animal_key.to_numpy(str)
    predictions = {k: np.asarray([STATES.index(s) for s in frame.predicted_state], dtype=int) for k, frame in frames.items()}
    animals = sorted(set(groups))
    indices = {a: np.flatnonzero(groups == a) for a in animals}
    rng = np.random.default_rng(metadata['bootstrap']['seed'])
    samples = []
    for _ in range(metadata['bootstrap']['replicates']):
        selected = np.concatenate([indices[a] for a in rng.choice(animals, len(animals), replace=True)])
        if np.any(np.bincount(y[selected], minlength=4) == 0):
            continue
        samples.append([balanced_accuracy(y[selected], predictions[k][selected]) for k in METHODS])
    draws = np.asarray(samples)
    intervals = np.percentile(draws, [2.5, 97.5], axis=0)
    result = {'n_cells': len(reference), 'n_animals': len(animals), 'bootstrap_replicates_used': len(draws), 'methods': {}}
    for i, method in enumerate(METHODS):
        p = predictions[method]
        matrix = confusion(y, p)
        values = {'balanced_accuracy': balanced_accuracy(y, p), 'accuracy': float((y == p).mean()),
                  'ci_low': float(intervals[0, i]), 'ci_high': float(intervals[1, i]),
                  'correct_n': int(np.trace(matrix)), 'confusion_matrix': matrix.tolist(),
                  'support': matrix.sum(axis=1).tolist(), 'recall': (np.diag(matrix) / matrix.sum(axis=1)).tolist()}
        expected = metadata['expected_metrics'][method]
        np.testing.assert_allclose([values['balanced_accuracy'], values['accuracy']], [expected['balanced_accuracy'], expected['accuracy']], atol=1e-12, rtol=0)
        np.testing.assert_allclose(intervals[:, i], expected['balanced_accuracy_95_interval'], atol=1e-12, rtol=0)
        result['methods'][method] = values
    delta = np.percentile(draws[:, 0] - draws[:, 1], [2.5, 97.5])
    expected_delta = -np.asarray(metadata['expected_metrics']['dino_logistic']['paired_difference_vs_morphometry_95_interval'])[::-1]
    np.testing.assert_allclose(delta, expected_delta, atol=1e-12, rtol=0)
    result['baseline_minus_dinov2'] = {
        'balanced_accuracy_difference': result['methods'][METHODS[0]]['balanced_accuracy'] - result['methods'][METHODS[1]]['balanced_accuracy'],
        'ci_low': float(delta[0]), 'ci_high': float(delta[1]),
        'interval_conditions_on_fixed_predictions': True}
    return result, draws


def draw(result, spanish=False):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .8, 'savefig.facecolor': 'white'})
    states = ['Ramificada', 'En bastón', 'Activada', 'Ameboide'] if spanish else list(STATES)
    names = ['Propuestas morfométricas', 'DINOv2 + clasificador lineal'] if spanish else ['Morphometric proposals', 'DINOv2 + linear classifier']
    fig = plt.figure(figsize=(11.2, 8.8), facecolor='white')
    grid = fig.add_gridspec(2, 2, left=.145, right=.975, top=.845, bottom=.16,
                           wspace=.58, hspace=.80, height_ratios=[.70, 1.])
    axes = [fig.add_subplot(grid[i, j]) for i in range(2) for j in range(2)]
    for letter, ax in zip('ABCD', axes):
        ax.text(-.18, 1.10, letter, transform=ax.transAxes, fontsize=16, weight='bold', ha='left', va='top')
    title = ('Figura suplementaria S6 | Clasificación de la microglía' if spanish
             else 'Supplementary Figure S6 | Microglial classification')
    subtitle = (f"{result['n_cells']} células revisadas por una persona · {result['n_animals']} animales · cuatro clases morfológicas" if spanish
                else f"{result['n_cells']} human-reviewed cells · {result['n_animals']} animals · four morphology classes")
    fig.text(.5, .967, title, ha='center', va='top', fontsize=16, weight='bold', color='#213943')
    fig.text(.5, .925, subtitle, ha='center', va='top', fontsize=10.5, color='#46585e')
    ax = axes[0]
    for i, method in enumerate(METHODS):
        val = result['methods'][method]
        y = 1 - i
        value, lo, hi = [100 * val[k] for k in ('balanced_accuracy', 'ci_low', 'ci_high')]
        ax.errorbar(value, y, xerr=[[value-lo], [hi-value]], fmt='o' if i == 0 else 'D',
                    markersize=8, color=COLORS[i], linewidth=2.1, capsize=5, zorder=3)
        ax.text(value, y+.22, f'{value:.1f}%', color=COLORS[i], ha='center', weight='bold', fontsize=13)
    ax.set(yticks=[1, 0], yticklabels=['Propuestas\nmorfométricas', 'DINOv2 +\nclasificador lineal'] if spanish
           else ['Morphometric\nproposals', 'DINOv2 +\nlinear classifier'],
           xlim=(0, 100), ylim=(-.48, 1.60), xticks=[0, 25, 50, 75, 100],
           xlabel='Exactitud balanceada (%)' if spanish else 'Balanced accuracy (%)')
    ax.set_title('Concordancia con la revisión humana' if spanish else 'Agreement with human review', loc='left', fontsize=11, weight='bold', pad=13)
    ax.grid(axis='x', color='#e5e9ea', linewidth=.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis='y', length=0, labelsize=9)
    ax.spines['left'].set_visible(False)
    delta = result['baseline_minus_dinov2']
    diff, lo, hi = [100 * delta[k] for k in ('balanced_accuracy_difference', 'ci_low', 'ci_high')]
    text = (f'Diferencia: +{diff:.1f} puntos porcentuales\nIC del 95% pareado: {lo:.1f} a {hi:.1f}' if spanish
            else f'Difference: +{diff:.1f} percentage points\nPaired 95% CI: {lo:.1f} to {hi:.1f}')
    ax.text(0, -.42, text, transform=ax.transAxes, fontsize=8.5, va='top', color='#46585e')
    ax = axes[1]
    for i, state in enumerate(STATES):
        vals = [100 * result['methods'][m]['recall'][i] for m in METHODS]
        ax.plot(vals, [3-i, 3-i], color='#c5ccce', linewidth=2, zorder=1)
    for i, method in enumerate(METHODS):
        ax.scatter(np.asarray(result['methods'][method]['recall'])*100, np.arange(3, -1, -1),
                   marker='o' if i == 0 else 'D', s=48, color=COLORS[i], zorder=3,
                   label=('Morfometría' if spanish else 'Morphometry') if i == 0 else 'DINOv2')
    support = result['methods'][METHODS[0]]['support']
    ax.set(yticks=range(3, -1, -1), yticklabels=[f'{s} (n={n})' for s, n in zip(states, support)],
           xlim=(0, 100), ylim=(-.45, 3.55), xticks=[0, 25, 50, 75, 100],
           xlabel='Sensibilidad por clase (%)' if spanish else 'Recall within human class (%)')
    ax.set_title('Concordancia por clase' if spanish else 'Performance by class', loc='left', fontsize=11, weight='bold', pad=13)
    ax.grid(axis='x', color='#e5e9ea', linewidth=.8)
    ax.set_axisbelow(True)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, labelsize=9)
    ax.legend(loc='upper left', bbox_to_anchor=(-.02, -.32), ncol=2, frameon=False, fontsize=8.5,
              handletextpad=.4, columnspacing=1.)
    for ax, method, name in zip(axes[2:], METHODS, names):
        matrix = np.asarray(result['methods'][method]['confusion_matrix'])
        percent = matrix / matrix.sum(axis=1, keepdims=True) * 100
        ax.imshow(percent, cmap='Blues', vmin=0, vmax=100, interpolation='nearest')
        for y in range(4):
            for x in range(4):
                ax.text(x, y, f'{matrix[y,x]}\n{percent[y,x]:.0f}%', ha='center', va='center', fontsize=10,
                        color='white' if percent[y,x] >= 60 else '#213943')
        ax.set(xticks=range(4), yticks=range(4), xticklabels=states, yticklabels=states,
               xlabel='Clase predicha' if spanish else 'Predicted class',
               ylabel='Clase revisada por una persona' if spanish else 'Human-reviewed class')
        ax.tick_params(axis='x', rotation=25, labelsize=8.5, length=0)
        ax.tick_params(axis='y', labelsize=8.5, length=0)
        ax.set_title(name, fontsize=11, weight='bold', pad=11)
        for spine in ax.spines.values():
            spine.set_visible(False)
    footer = ('Exactitud balanceada: media de las sensibilidades por clase. Intervalos: remuestreo de animales.\nValidación interna de una muestra seleccionada para revisión; las propuestas morfométricas permanecen fijas.' if spanish
              else 'Balanced accuracy: mean recall across classes. Intervals: bootstrap resampling of animals.\nInternal validation of a selected review sample; the morphometric proposals are fixed.')
    fig.text(.5, .025, footer, ha='center', va='bottom', fontsize=8, color='#46585e')
    fig.canvas.draw()
    return fig, axes


def captions(result, metadata, spanish=False):
    base, dino = [result['methods'][k] for k in METHODS]
    delta = result['baseline_minus_dinov2']
    if not spanish:
        title = 'Supplementary Figure S6. Internal benchmark of microglial morphology classification.'
        panels = [
            f"(A) Balanced accuracy, defined as the unweighted mean recall across Ramified, Rod-like, Activated and Amoeboid human-reviewed classes. Existing morphometric proposals achieved {100*base['balanced_accuracy']:.1f}% (95% interval {100*base['ci_low']:.1f}–{100*base['ci_high']:.1f}%), compared with {100*dino['balanced_accuracy']:.1f}% ({100*dino['ci_low']:.1f}–{100*dino['ci_high']:.1f}%) for frozen DINOv2 features and a linear classifier. The paired baseline-minus-DINOv2 difference was {100*delta['balanced_accuracy_difference']:.1f} percentage points (95% interval {100*delta['ci_low']:.1f}–{100*delta['ci_high']:.1f}).",
            '(B) Recall for each human-reviewed class, evaluated on the same cells for both methods. Connecting segments pair the two methods within a class. The reference counts are Ramified n=50, Rod-like n=32, Activated n=50 and Amoeboid n=71.',
            '(C) Confusion matrix for the frozen morphometric proposals. Rows indicate the human-reviewed class and columns the predicted class. Each entry gives a cell count and its percentage of that row; color intensity represents row percentage on the same 0–100% scale as panel D.',
            '(D) Confusion matrix for the DINOv2 classifier using predictions made with the entire target animal excluded from supervised fitting, preprocessing and hyperparameter selection. Counts, row normalization and color scale match panel C.'
        ]
        methods = ('All panels use the same 203 human-reviewed cells from 15 animals in the Figure 5 soma-refined segmentation run. The review campaign selected 300 targets in equal proposed-class quotas; 203 reviews were completed. The morphometric comparator is the pre-existing four-cluster K-means proposal set, using standardized cell/branch morphology and a predefined cluster-to-class mapping. No human class labels fit that comparator; its unsupervised preprocessing and clustering used the full 14,415-cell collection, including images from the evaluated animals. It is evaluated as a fixed comparator rather than as a newly fitted supervised model. '
                   'DINOv2 ViT-S/14 uses externally pretrained, frozen weights. Close target and wider branch-context Iba1 crops provide class-token and target-mask-weighted patch-token features, averaged across four right-angle rotations. Within each outer leave-one-animal-out fold, feature standardization, 32-component whitened PCA, balanced multinomial logistic regression, and selection of regularization strength by three-fold grouped inner validation use only the training animals. All 15 outer folds are pooled for the reported metrics. '
                   'Intervals are percentile 95% intervals from 2,000 paired bootstrap draws of whole animals with replacement (seed 20260907), recomputing the pooled class-balanced metric while retaining the fitted out-of-fold predictions. They describe uncertainty conditional on these predictions and do not include variability from retraining. '
                   f"Ordinary cell-level accuracy, a separate metric, was {100*base['accuracy']:.1f}% ({base['correct_n']}/203) for morphometric proposals and {100*dino['accuracy']:.1f}% ({dino['correct_n']}/203) for DINOv2. "
                   'The selected review sample and prior model development limit inference: this is internal validation against one human review, not independent biological validation or evidence about population prevalence. The comparison does not establish general superiority across microscopy datasets or validate the morphological classes as functional activation states. Original microscopy, human labels and Figure 5 publication outputs were preserved.')
    else:
        title = 'Figura suplementaria S6. Comparación interna de la clasificación morfológica de la microglía.'
        panels = [
            f"(A) Exactitud balanceada, definida como la media no ponderada de la sensibilidad de las clases Ramificada, En bastón, Activada y Ameboide revisadas por una persona. Las propuestas morfométricas alcanzaron {100*base['balanced_accuracy']:.1f}% (intervalo del 95%: {100*base['ci_low']:.1f}–{100*base['ci_high']:.1f}%), frente a {100*dino['balanced_accuracy']:.1f}% ({100*dino['ci_low']:.1f}–{100*dino['ci_high']:.1f}%) para DINOv2 congelado y un clasificador lineal. La diferencia pareada fue de {100*delta['balanced_accuracy_difference']:.1f} puntos porcentuales (intervalo del 95%: {100*delta['ci_low']:.1f}–{100*delta['ci_high']:.1f}).",
            '(B) Sensibilidad por clase revisada, calculada sobre las mismas células para ambos métodos. Las líneas conectan los métodos dentro de cada clase. Recuentos de referencia: Ramificada n=50, En bastón n=32, Activada n=50 y Ameboide n=71.',
            '(C) Matriz de confusión de las propuestas morfométricas fijas. Las filas representan la clase revisada y las columnas la predicción. Cada entrada muestra el recuento celular y su porcentaje dentro de la fila; la escala de color de 0–100% coincide con D.',
            '(D) Matriz de confusión de DINOv2, con el animal completo de cada célula excluido del ajuste supervisado, del preprocesamiento aprendido y de la selección de hiperparámetros. Los recuentos y la normalización coinciden con C.'
        ]
        methods = ('Los paneles incluyen las mismas 203 células revisadas de 15 animales de la segmentación con somas refinados de la Figura 5. La campaña seleccionó 300 objetivos con cuotas iguales según la clase propuesta; se completaron 203 revisiones. El comparador morfométrico corresponde a propuestas previas obtenidas mediante K-means de cuatro grupos, medidas estandarizadas de células y prolongaciones y una asignación predefinida de grupos a clases. No se utilizaron etiquetas humanas para ajustar este comparador; su preprocesamiento no supervisado y agrupamiento incluyeron las 14.415 células, también las imágenes de los animales evaluados. Se evalúa como comparador fijo. '
                   'DINOv2 ViT-S/14 utiliza pesos preentrenados externamente y congelados. Una vista cercana al objetivo y otra del contexto de las prolongaciones Iba1 aportan características de tokens globales y de parches ponderados por la máscara celular, promediadas en cuatro rotaciones de ángulo recto. La estandarización, el PCA blanqueado de 32 componentes, la regresión logística multinomial balanceada y la selección de regularización mediante tres particiones internas por animal usan exclusivamente los animales de entrenamiento de cada una de las 15 particiones externas. '
                   'Los intervalos percentiles del 95% usan 2.000 remuestreos pareados de animales completos con reemplazo (semilla 20260907). Se recalculan las métricas agrupadas manteniendo fijas las predicciones de validación; los intervalos no incluyen la variabilidad del reentrenamiento. '
                   f"La exactitud celular ordinaria, una métrica distinta, fue {100*base['accuracy']:.1f}% ({base['correct_n']}/203) para morfometría y {100*dino['accuracy']:.1f}% ({dino['correct_n']}/203) para DINOv2. "
                   'La selección de la muestra y el desarrollo previo limitan la inferencia: esta es una validación interna frente a una revisión humana, sin validación biológica independiente ni estimación de prevalencia poblacional. No demuestra superioridad general entre conjuntos de microscopía ni que las clases morfológicas definan estados funcionales de activación. Se conservaron las imágenes originales, las etiquetas humanas y los archivos publicados de la Figura 5.')
    return title, panels, methods


def export(args, metadata, reference, frames, manifest, result, draws):
    out = args.output_dir.expanduser().resolve()
    for protected in (PAPER / 'FigS6', args.input_dir.resolve()):
        if out == protected or out in protected.parents or protected in out.parents:
            raise ValueError('Render into a fresh stage outside the figure/input directories.')
    out.mkdir(parents=True, exist_ok=False)
    for name in ('panels', 'legends', 'source_data', 'provenance'):
        (out / name).mkdir()
    english, english_axes = draw(result)
    spanish, spanish_axes = draw(result, spanish=True)
    eng_renderer = english.canvas.get_renderer()
    spa_renderer = spanish.canvas.get_renderer()
    # Identical panel bounds in both languages include every tick/legend artist.
    bounds = []
    for a, b in zip(english_axes, spanish_axes):
        bounds.append(Bbox.union([a.get_tightbbox(eng_renderer).transformed(english.dpi_scale_trans.inverted()),
                                  b.get_tightbbox(spa_renderer).transformed(spanish.dpi_scale_trans.inverted())]).padded(.08))
    for suffix in ('pdf', 'png'):
        english.savefig(out / f'Figure_S6.{suffix}', dpi=args.dpi)
    for fig, axes, language in ((english, english_axes, ''), (spanish, spanish_axes, '_spanish')):
        for index, (slug, bbox) in enumerate(zip(SLUGS, bounds)):
            previous = [a.get_visible() for a in axes]
            text_previous = [t.get_visible() for t in fig.texts]
            try:
                for j, ax in enumerate(axes):
                    ax.set_visible(j == index)
                for artist in fig.texts:
                    artist.set_visible(False)
                for suffix in ('pdf', 'png'):
                    fig.savefig(out / 'panels' / f'Figure_S6_panel_{slug}{language}.{suffix}', bbox_inches=bbox, dpi=args.dpi)
            finally:
                for ax, visible in zip(axes, previous):
                    ax.set_visible(visible)
                for artist, visible in zip(fig.texts, text_previous):
                    artist.set_visible(visible)
    plt.close(english)
    plt.close(spanish)
    for is_spanish, suffix in ((False, ''), (True, '_spanish')):
        title, panel_captions, methods = captions(result, metadata, is_spanish)
        full = '\n\n'.join([title, *panel_captions, methods]) + '\n'
        (out / 'legends' / f'Figure_S6_LEGEND{suffix}.txt').write_text(full)
        if not is_spanish:
            (out / 'Figure_S6_caption.txt').write_text(full)
        for slug, caption in zip(SLUGS, panel_captions):
            (out / 'legends' / f'Figure_S6_panel_{slug}_LEGEND{suffix}.txt').write_text('\n\n'.join([title, caption, methods]) + '\n')
    rows = []
    for method in METHODS:
        metrics = result['methods'][method]
        rows.append({'method': method, **{k: metrics[k] for k in ('balanced_accuracy', 'accuracy', 'ci_low', 'ci_high', 'correct_n')},
                     'cells_n': len(reference), 'animals_n': reference.animal_key.nunique()})
        pd.DataFrame(metrics['confusion_matrix'], index=STATES, columns=STATES).to_csv(out / 'source_data' / f'{method}_confusion.csv')
        pd.DataFrame({'human_state': STATES, 'n': metrics['support'], 'recall': metrics['recall']}).to_csv(out / 'source_data' / f'{method}_class_recall.csv', index=False)
    pd.DataFrame(rows).to_csv(out / 'source_data/balanced_accuracy_comparison.csv', index=False)
    joined = reference.copy()
    for method in METHODS:
        joined[method] = frames[method].predicted_state
    joined.to_csv(out / 'source_data/paired_cell_predictions.csv', index=False)
    pd.crosstab(reference.animal_key, reference.human_state).reindex(columns=STATES).to_csv(out / 'source_data/reference_cells_by_animal.csv')
    pd.DataFrame(draws, columns=METHODS).to_csv(out / 'source_data/animal_bootstrap_draws.csv', index=False)
    write_json(out / 'source_data/comparison_statistics.json', result)
    shutil.copy2(args.input_dir / 'method_metadata.json', out / 'source_data/method_metadata.json')
    source_paths = {str(args.input_dir / row.file): row.sha256 for row in manifest.itertuples(index=False)}
    outputs = {str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write_json(out / 'provenance/FIGURE_S6_BUILD_RECEIPT.json', {
        'figure': 'FigS6', 'created_at': datetime.now(timezone.utc).isoformat(),
        'input_files_sha256': source_paths, 'input_manifest_sha256': sha256(args.input_dir / 'RAW_DATA_MANIFEST.csv'),
        'renderer_sha256': sha256(Path(__file__)), 'dpi': args.dpi,
        'cells_n': len(reference), 'animals_n': reference.animal_key.nunique(),
        'validation': metadata['validation'], 'metric': 'balanced accuracy, not ordinary cell-level accuracy',
        'checks': ['frozen source SHA-256', 'same human-referenced cells for both methods', 'whole-animal outer and inner separation',
                   'each reviewed cell held out exactly once', 'all metrics and bootstrap intervals independently recomputed'],
        'output_files_sha256': outputs, 'existing_figures_modified': False})
    print(json.dumps({'status': 'PASS', 'output': str(out), 'balanced_accuracy': {m: result['methods'][m]['balanced_accuracy'] for m in METHODS},
                      'paired_difference': result['baseline_minus_dinov2']}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    if args.dpi != 600:
        parser.error('Publication master and panel exports require 600 dpi.')
    if not args.validate_only and args.output_dir is None:
        parser.error('--output-dir is required for a fresh render.')
    args.input_dir = args.input_dir.expanduser().resolve()
    metadata, reference, frames, manifest = load_inputs(args.input_dir)
    result, draws = calculate(metadata, reference, frames)
    if args.validate_only:
        print(json.dumps({'status': 'PASS', 'cells_n': len(reference), 'animals_n': reference.animal_key.nunique(),
                          'balanced_accuracy': {m: result['methods'][m]['balanced_accuracy'] for m in METHODS}}, indent=2))
        return
    export(args, metadata, reference, frames, manifest, result, draws)


if __name__ == '__main__':
    main()
