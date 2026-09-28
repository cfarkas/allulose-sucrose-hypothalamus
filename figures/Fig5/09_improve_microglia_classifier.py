#!/usr/bin/env python3
"""Figure 5: frozen DINOv2 features, reviewed masks, and animal-blocked learning.

Run prepare/train/report with paper_apotome_repro; run extract on cellpose_env
for the compatible CUDA stack. All outputs are versioned, and no original
microscopy, human labels, accepted masks, or publication figures are modified.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import inspect
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
from skimage.measure import regionprops
from skimage.morphology import skeletonize

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/06_retrain_microglia_classes.py').is_file())
DEFAULT_REVIEW = PAPER / 'Fig5/hil_review/microglia_masks_20260907_v2'
DEFAULT_ASSETS = PAPER / 'analyses/Fig5/model_assets/dinov2'
STATES = ('Ramified', 'Rod-like', 'Activated', 'Amoeboid')
SEED = 20260907
DINO_WIDTH = 384
EMBEDDING_WIDTH = 4 * DINO_WIDTH
RECIPE = {'version': 'target_context_dino_v1', 'input_size': 224, 'tight_min_px': 24,
          'tight_margin_fraction': .25, 'context_min_px': 64, 'context_scale': 2.,
          'outside_signal_weight': .20, 'target_dilation_px': 2, 'target_blur_sigma': 1.,
          'rotations_degrees': [0, 90, 180, 270], 'encoder': 'dinov2_vits14',
          'features_per_view': ['class_token', 'target_mask_weighted_patch_tokens'],
          'normalization_mean': [.485, .456, .406], 'normalization_std': [.229, .224, .225],
          'pca_components': 32, 'uncertain_score_below': .60}
FAMILIES = ('morphology_logistic', 'morphology_forest', 'dino_logistic', 'hybrid_logistic')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, value):
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_cell(row, verify=True):
    patch = Path(row['patch_path'])
    if verify and sha256(patch) != row['patch_sha256']:
        raise ValueError(f'Patch changed: {row["cell_uid"]}')
    with Image.open(patch) as image:
        rgb = np.array(image.convert('RGB'))
    cell_path = str(row.get('cell_mask_path', ''))
    if cell_path:
        if verify and sha256(cell_path) != row['cell_mask_sha256']:
            raise ValueError(f'Cell mask changed: {row["cell_uid"]}')
        with Image.open(cell_path) as image:
            cell = np.array(image) > 127
    else:
        cell = rgb[..., 1] > 127
    if verify and sha256(row['soma_mask_path']) != row['soma_mask_sha256']:
        raise ValueError(f'Soma mask changed: {row["cell_uid"]}')
    with Image.open(row['soma_mask_path']) as image:
        soma = np.array(image) > 127
    if cell.shape != rgb.shape[:2] or soma.shape != cell.shape:
        raise ValueError('Mask and image dimensions differ.')
    if not cell.any() or not soma.any() or np.any(soma & ~cell):
        raise ValueError('Nonempty cell/soma masks with soma inside cell are required.')
    if not np.array_equal(cell, rgb[..., 1] > 127):
        raise ValueError('Review mask and original CNN mask channel differ.')
    return rgb[..., 0], cell, soma


def square_crop(array, center, side):
    """Crop with background padding; never resize or wrap outside image bounds."""
    side = int(side)
    y0, x0 = [int(math.floor(c - side / 2)) for c in center]
    y1, x1 = y0 + side, x0 + side
    result = np.zeros((side, side), dtype=array.dtype)
    sy0, sx0 = max(0, y0), max(0, x0)
    sy1, sx1 = min(array.shape[0], y1), min(array.shape[1], x1)
    if sy0 < sy1 and sx0 < sx1:
        result[sy0-y0:sy1-y0, sx0-x0:sx1-x0] = array[sy0:sy1, sx0:sx1]
    return result


def native_views(signal, cell, soma):
    yy, xx = np.nonzero(cell)
    center = ((yy.min() + yy.max() + 1) / 2, (xx.min() + xx.max() + 1) / 2)
    extent = int(max(yy.max() - yy.min() + 1, xx.max() - xx.min() + 1))
    tight_side = max(RECIPE['tight_min_px'], int(math.ceil(extent * (1 + 2 * RECIPE['tight_margin_fraction']))))
    context_side = max(RECIPE['context_min_px'], int(math.ceil(extent * RECIPE['context_scale'])))
    envelope = ndi.gaussian_filter(ndi.binary_dilation(cell, iterations=RECIPE['target_dilation_px']).astype(np.float32), RECIPE['target_blur_sigma'])
    envelope[cell] = 1.
    focus = signal.astype(np.float32) * (RECIPE['outside_signal_weight'] + (1 - RECIPE['outside_signal_weight']) * envelope)
    views = [(square_crop(focus, center, tight_side), square_crop(cell, center, tight_side)),
             (square_crop(signal, center, context_side), square_crop(cell, center, context_side))]
    for _, mask in views:
        if int(mask.sum()) != int(cell.sum()):
            raise ValueError('A target branch was clipped by a generated crop.')
    return views


def resized_views(signal, cell, soma):
    output = []
    for view, mask in native_views(signal, cell, soma):
        image = Image.fromarray(np.clip(view, 0, 255).astype(np.uint8)).resize((RECIPE['input_size'],) * 2, Image.Resampling.BICUBIC)
        # No overlay or class-dependent colors enter the pretrained encoder.
        rgb = np.repeat(np.array(image)[..., None], 3, axis=2)
        weights = np.array(Image.fromarray(mask.astype(np.float32)).resize((RECIPE['input_size'],) * 2, Image.Resampling.BILINEAR))
        output.append((rgb, weights))
    return output


def shape_features(signal, cell, soma):
    features = {}
    for name, mask in [('cell', cell), ('soma', soma)]:
        region = regionprops(mask.astype(np.uint8))[0]
        area = float(mask.sum())
        perimeter = float(region.perimeter)
        features.update({
            f'{name}_log_area': np.log1p(area), f'{name}_log_perimeter': np.log1p(perimeter),
            f'{name}_log_major_axis': np.log1p(region.axis_major_length),
            f'{name}_log_minor_axis': np.log1p(region.axis_minor_length),
            f'{name}_aspect_ratio': region.axis_major_length / max(region.axis_minor_length, 1.),
            f'{name}_eccentricity': float(region.eccentricity),
            f'{name}_circularity': 4 * math.pi * area / max(perimeter ** 2, 1.),
            f'{name}_solidity': float(region.solidity), f'{name}_extent': float(region.extent),
        })
        intensity = signal[mask].astype(float) / 255
        features[f'{name}_mean_intensity'] = float(intensity.mean())
        features[f'{name}_intensity_sd'] = float(intensity.std())
    skeleton = skeletonize(cell)
    degree = ndi.convolve(skeleton.astype(np.int8), np.ones((3, 3), np.int8), mode='constant') - skeleton
    length = int(skeleton.sum())
    ends = int(np.sum(skeleton & (degree == 1)))
    branches = int(ndi.label(skeleton & (degree >= 3))[1])
    processes = cell & ~soma
    features.update(soma_fraction=float(soma.sum() / cell.sum()),
                    process_log_area=float(np.log1p(processes.sum())),
                    skeleton_log_length=float(np.log1p(length)),
                    skeleton_length_norm=float(length / np.sqrt(cell.sum())),
                    skeleton_endpoints=float(ends), skeleton_branch_regions=float(branches),
                    endpoints_per_length=float(ends / max(length, 1)),
                    branches_per_length=float(branches / max(length, 1)),
                    mean_process_width=float(processes.sum() / max(np.sum(skeleton & processes), 1)))
    cy, cx = ndi.center_of_mass(soma)
    yy, xx = np.indices(cell.shape)
    distance = np.hypot(yy - cy, xx - cx)
    crossings = []
    for radius in (5, 10, 15, 20, 30):
        count = 0
        for dy, dx in ((0, 1), (1, -1), (1, 0), (1, 1)):
            sy = slice(0, cell.shape[0] - dy)
            sx = slice(max(0, -dx), min(cell.shape[1], cell.shape[1] - dx))
            ty = slice(dy, cell.shape[0])
            tx = slice(max(0, dx), min(cell.shape[1], cell.shape[1] + dx))
            edge = skeleton[sy, sx] & skeleton[ty, tx]
            count += int(np.sum(edge & ((distance[sy, sx] < radius) != (distance[ty, tx] < radius))))
        crossings.append(count)
    features['sholl_max_crossings'] = float(max(crossings))
    features['sholl_mean_crossings'] = float(np.mean(crossings))
    features['process_mean_intensity'] = float(signal[processes].mean() / 255) if processes.any() else 0.
    ring = ndi.binary_dilation(cell, iterations=5) & ~cell
    features['local_intensity_contrast'] = float((signal[cell].mean() - (signal[ring].mean() if ring.any() else 0)) / 255)
    if not np.isfinite(list(features.values())).all():
        raise ValueError('Nonfinite morphological feature.')
    return features


def recipe_hash():
    text = json.dumps(RECIPE, sort_keys=True) + ''.join(inspect.getsource(f) for f in (read_cell, square_crop, native_views, resized_views, shape_features))
    return hashlib.sha256(text.encode()).hexdigest()


def configurations():
    output = []
    for family in FAMILIES:
        if family == 'morphology_forest':
            output.extend({'family': family, 'min_samples_leaf': leaf} for leaf in (2, 5))
        else:
            output.extend({'family': family, 'C': c} for c in (.1, 1., 10.))
    return output


def prepare(args):
    trainer = load_module('fig5_transfer_review_validator', PAPER / 'Fig5/06_retrain_microglia_classes.py')
    review = args.review_dir.resolve()
    metadata, cells, _ = trainer.load_source(review)
    labels_path = review / 'labels_used_for_training.csv'
    labels, lab, groups, tr, va, seed = trainer.validate_labels(labels_path, metadata, cells, 200)
    if labels_path.read_bytes() != (review / 'microglia_reviewed_labels.csv').read_bytes():
        raise ValueError('The frozen training labels differ from the current review.')
    assets = json.loads((args.assets_dir / 'asset_manifest.json').read_text())
    if sha256(assets['weights_path']) != assets['weights_sha256']:
        raise ValueError('Encoder weights changed.')
    out = args.run_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / 'source_code').mkdir()
    shutil.copy2(Path(__file__).resolve(), out / 'source_code/09_improve_microglia_classifier.py')
    shutil.copy2(labels_path, out / 'human_labels_used.csv')
    shutil.copy2(review / 'review_batch.json', out / 'review_batch_snapshot.json')
    lookup = {r['cell_uid']: r for r in metadata['selected_cells']}
    reviewed = set(labels.cell_uid)
    rows, features = [], []
    for i, row in enumerate(cells.to_dict('records')):
        uid = row['cell_uid']
        selected = uid in lookup
        soma_path = review / 'masks' / f'{uid}_soma.png' if selected else Path(row['soma_patch_path'])
        cell_path = review / 'masks' / f'{uid}_cell.png' if selected else None
        item = {k: str(row[k]) for k in ('cell_uid', 'sample', 'section', 'animal_key', 'region', 'patch_path')}
        item.update(patch_sha256=sha256(item['patch_path']), soma_mask_path=str(soma_path),
                    soma_mask_sha256=sha256(soma_path), cell_mask_path=str(cell_path) if cell_path else '',
                    cell_mask_sha256=sha256(cell_path) if cell_path else '', is_reviewed=uid in reviewed)
        signal, cell, soma = read_cell(item)
        features.append({'cell_uid': uid, **shape_features(signal, cell, soma)})
        rows.append(item)
        if (i + 1) % 2000 == 0:
            print(f'Prepared {i + 1}/{len(cells)} cell/mask feature records', flush=True)
    manifest = pd.DataFrame(rows)
    manifest.to_csv(out / 'cells_manifest.csv', index=False)
    shape = pd.DataFrame(features)
    shape.to_csv(out / 'shape_features.csv', index=False)
    legacy = lab[['state']].copy()
    legacy['animal_key'] = groups
    legacy['split'] = 'training'
    legacy.loc[lab.iloc[va].index, 'split'] = 'validation'
    legacy.to_csv(out / 'legacy_development_split.csv')
    protocol = {
        'created_at': datetime.now(timezone.utc).isoformat(), 'review_dir': str(review),
        'source_analysis': metadata['source_analysis'], 'source_measurements_sha256': metadata['source_measurements_sha256'],
        'human_labels_sha256': sha256(out / 'human_labels_used.csv'), 'review_batch_sha256': sha256(out / 'review_batch_snapshot.json'),
        'cells_manifest_sha256': sha256(out / 'cells_manifest.csv'), 'shape_features_sha256': sha256(out / 'shape_features.csv'),
        'cells_n': len(cells), 'reviewed_cells_n': len(labels), 'reviewed_animals_n': len(set(groups)),
        'class_counts': dict(Counter(labels.state)), 'shape_features': shape.columns[1:].tolist(),
        'recipe': RECIPE, 'recipe_sha256': recipe_hash(), 'assets': assets, 'seed': SEED,
        'configurations': configurations(),
        'validation': 'Nested leave-one-animal-out outer validation; 3-fold stratified-group inner selection. Scaling, PCA, class weights and all model choices fit on training animals only.',
        'selection_rule': 'Highest pooled inner balanced accuracy, then lowest inner log loss, then configuration order. Final candidate selected by grouped CV on all human labels, never by the legacy 40-cell result.',
        'legacy_comparison': 'The same historical 163/40 split is a development comparison, not an independent test set.',
        'input_exclusions': ['animal identity', 'region', 'treatment', 'sex', 'genotype', 'coordinates', 'proposed class', 'proposed confidence'],
        'interpretation': 'Internal validation of a proposed-class-balanced human-review sample; no independent final test cohort and no unbiased prevalence claim.'}
    atomic_json(out / 'protocol.json', protocol)
    print(f'Frozen {len(labels)} reviewed cells; {len(cells)} source cells; {len(shape.columns)-1} mask/intensity features.', flush=True)
    render_input_examples(out, manifest, labels)


def verify_run(run):
    protocol = json.loads((run / 'protocol.json').read_text())
    for name, key in [('cells_manifest.csv', 'cells_manifest_sha256'), ('shape_features.csv', 'shape_features_sha256'), ('human_labels_used.csv', 'human_labels_sha256'), ('review_batch_snapshot.json', 'review_batch_sha256')]:
        if sha256(run / name) != protocol[key]:
            raise ValueError(f'Frozen experiment input changed: {name}')
    if recipe_hash() != protocol['recipe_sha256']:
        raise ValueError('Feature extraction recipe changed; use a fresh run directory.')
    return protocol


def render_input_examples(run, manifest, labels):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    chosen = labels.groupby('state', sort=False).head(2).cell_uid.tolist()
    lookup = manifest.set_index('cell_uid')
    fig, axes = plt.subplots(len(chosen), 3, figsize=(8, 2.05 * len(chosen)))
    label_lookup = labels.set_index('cell_uid').state
    for i, uid in enumerate(chosen):
        signal, cell, soma = read_cell(lookup.loc[uid].to_dict() | {'cell_uid': uid})
        views = resized_views(signal, cell, soma)
        axes[i, 0].imshow(signal, cmap='gray', vmin=0, vmax=255)
        axes[i, 0].contour(cell, levels=[.5], colors=['#13b5c6'], linewidths=.6)
        axes[i, 0].set_title(f'Human: {label_lookup[uid]}', fontsize=9)
        for j, (rgb, _) in enumerate(views, 1):
            axes[i, j].imshow(rgb)
            axes[i, j].set_title(['Closer target', 'Wider branch context'][j-1], fontsize=9)
        for axis in axes[i]:
            axis.axis('off')
    fig.suptitle('Figure 5 input crops: source fluorescence and actual encoder inputs', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .98])
    fig.savefig(run / 'input_crop_examples.png', dpi=160)
    fig.savefig(run / 'input_crop_examples.pdf')
    plt.close(fig)


def extract(args):
    run = args.run_dir.resolve()
    protocol = verify_run(run)
    assets = protocol['assets']
    if sha256(assets['weights_path']) != assets['weights_sha256']:
        raise ValueError('Pretrained encoder weights changed.')
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    import torch
    from torch.nn import functional as F
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    sys.path.insert(0, assets['source_dir'])
    from dinov2.hub.backbones import dinov2_vits14
    model = dinov2_vits14(pretrained=False)
    model.load_state_dict(torch.load(assets['weights_path'], map_location='cpu', weights_only=True), strict=True)
    model.requires_grad_(False).eval().to(args.device)
    manifest = pd.read_csv(run / 'cells_manifest.csv', keep_default_na=False)
    path = run / 'dino_features.npy'
    if not path.exists():
        features = np.lib.format.open_memmap(path, mode='w+', dtype='float32', shape=(len(manifest), EMBEDDING_WIDTH))
        features[:] = np.nan
        features.flush()
    else:
        features = np.lib.format.open_memmap(path, mode='r+')
        if features.shape != (len(manifest), EMBEDDING_WIDTH):
            raise ValueError('Cached embedding shape differs from the frozen manifest.')
    wanted = np.ones(len(manifest), bool) if args.all_cells else manifest.is_reviewed.astype(str).str.lower().eq('true').to_numpy()
    pending = np.flatnonzero(wanted & ~np.isfinite(features).all(axis=1))
    mean = torch.tensor(RECIPE['normalization_mean'], device=args.device)[None, :, None, None]
    std = torch.tensor(RECIPE['normalization_std'], device=args.device)[None, :, None, None]
    begin = time.monotonic()
    for start in range(0, len(pending), args.cells_per_batch):
        indices = pending[start:start + args.cells_per_batch]
        images, masks = [], []
        for idx in indices:
            signal, cell, soma = read_cell(manifest.iloc[int(idx)].to_dict())
            for rgb, weights in resized_views(signal, cell, soma):
                for rotation in range(4):
                    images.append(np.rot90(rgb, rotation).copy())
                    masks.append(np.rot90(weights, rotation).copy())
        x = torch.from_numpy(np.stack(images).transpose(0, 3, 1, 2)).to(args.device, dtype=torch.float32) / 255
        weights = torch.from_numpy(np.stack(masks)).to(args.device)
        batches = []
        with torch.inference_mode():
            for offset in range(0, len(x), args.image_batch_size):
                encoded = model.forward_features((x[offset:offset+args.image_batch_size] - mean) / std)
                tokens = encoded['x_norm_patchtokens']
                side = int(math.isqrt(tokens.shape[1]))
                pooled_mask = F.interpolate(weights[offset:offset+args.image_batch_size, None], size=(side, side), mode='area').flatten(1)
                masked = (tokens * pooled_mask[..., None]).sum(1) / pooled_mask.sum(1).clamp_min(1e-6)[..., None]
                batches.append(torch.cat((encoded['x_norm_clstoken'], masked), dim=1).cpu().numpy())
        result = np.concatenate(batches).reshape(len(indices), 2, 4, 2, DINO_WIDTH).mean(axis=2)
        result /= np.maximum(np.linalg.norm(result, axis=-1, keepdims=True), 1e-12)
        if not np.isfinite(result).all():
            raise ValueError('Nonfinite image embedding.')
        features[indices] = result.reshape(len(indices), EMBEDDING_WIDTH)
        features.flush()
        if start == 0 or (start // args.cells_per_batch + 1) % 10 == 0 or start + len(indices) == len(pending):
            complete = int(np.isfinite(features).all(axis=1).sum())
            atomic_json(run / 'extraction_status.json', {'complete_cells': complete, 'total_cells': len(manifest), 'requested_pending_cells': len(pending), 'processed_this_run': start + len(indices), 'seconds': time.monotonic()-begin, 'device': args.device, 'torch': torch.__version__, 'recipe_sha256': recipe_hash(), 'encoder_weights_sha256': assets['weights_sha256']})
            print(f'DINO features {start + len(indices)}/{len(pending)} requested; {complete}/{len(manifest)} cached; {time.monotonic()-begin:.1f}s', flush=True)
    print('Requested pretrained features are complete.', flush=True)


def load_matrix(run, reviewed_only):
    protocol = verify_run(run)
    manifest = pd.read_csv(run / 'cells_manifest.csv', keep_default_na=False)
    shape = pd.read_csv(run / 'shape_features.csv')
    if not manifest.cell_uid.equals(shape.cell_uid):
        raise ValueError('Cell and morphology row order differs.')
    embeddings = np.load(run / 'dino_features.npy', mmap_mode='r')
    ids = np.flatnonzero(manifest.is_reviewed.astype(str).str.lower().eq('true')) if reviewed_only else np.arange(len(manifest))
    selected = np.array(embeddings[ids])
    if not np.isfinite(selected).all():
        raise ValueError('Requested image features are incomplete; finish extraction first.')
    x = np.concatenate((shape.loc[ids, protocol['shape_features']].to_numpy(np.float64), selected), axis=1)
    return protocol, manifest.iloc[ids].reset_index(drop=True), x


def make_model(config, n_shape, n_training):
    from sklearn.compose import ColumnTransformer
    from sklearn.decomposition import PCA
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    family = config['family']
    columns = []
    if family.startswith('morphology') or family == 'hybrid_logistic':
        columns.append(('morphology', StandardScaler(), list(range(n_shape))))
    if family in ('dino_logistic', 'hybrid_logistic'):
        dino = Pipeline([('scale', StandardScaler()),
                         ('pca', PCA(n_components=min(RECIPE['pca_components'], n_training - 1), whiten=True, svd_solver='full'))])
        columns.append(('dino', dino, list(range(n_shape, n_shape + EMBEDDING_WIDTH))))
    transform = ColumnTransformer(columns, remainder='drop')
    if family == 'morphology_forest':
        classifier = RandomForestClassifier(n_estimators=200, min_samples_leaf=config['min_samples_leaf'],
                                            max_features=1., class_weight='balanced', random_state=SEED, n_jobs=4)
    else:
        classifier = LogisticRegression(C=config['C'], class_weight='balanced', solver='lbfgs',
                                        max_iter=2000, tol=1e-6, random_state=SEED)
    return Pipeline([('features', transform), ('classifier', classifier)])


def probabilities(model, x):
    p = model.predict_proba(x)
    if set(model.classes_) != set(range(4)):
        raise ValueError('All four morphology classes must be represented in fitting.')
    return p[:, [list(model.classes_).index(i) for i in range(4)]]


def grouped_inner_splits(y, groups):
    from sklearn.model_selection import StratifiedGroupKFold
    splits = list(StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED).split(np.zeros(len(y)), y, groups))
    for fit, validation in splits:
        if set(groups[fit]) & set(groups[validation]):
            raise ValueError('Animal leakage in model selection.')
        if set(y[fit]) != set(range(4)):
            raise ValueError('Model-selection training fold lacks a human-reviewed class.')
    return splits


def scores(y, p):
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, log_loss
    predicted = p.argmax(axis=1)
    return {'balanced_accuracy': float(balanced_accuracy_score(y, predicted)),
            'accuracy': float(accuracy_score(y, predicted)),
            'log_loss': float(log_loss(y, p, labels=list(range(4)))),
            'classification_report': classification_report(y, predicted, labels=list(range(4)), target_names=STATES, output_dict=True, zero_division=0)}


def choose_config(x, y, groups, n_shape, candidates):
    splits = grouped_inner_splits(y, groups)
    results = []
    for ordinal, config in enumerate(candidates):
        p = np.full((len(y), 4), np.nan)
        for fit, validation in splits:
            model = make_model(config, n_shape, len(fit))
            model.fit(x[fit], y[fit])
            p[validation] = probabilities(model, x[validation])
        result = scores(y, p)
        results.append({'configuration': config, 'ordinal': ordinal,
                        'balanced_accuracy': result['balanced_accuracy'], 'log_loss': result['log_loss']})
    ranked = sorted(results, key=lambda r: (-r['balanced_accuracy'], r['log_loss'], r['ordinal']))
    best_by_family = {family: next(r['configuration'] for r in ranked if r['configuration']['family'] == family) for family in FAMILIES}
    return ranked[0]['configuration'], best_by_family, results


def prediction_frame(manifest, y, p):
    result = manifest[['cell_uid', 'animal_key']].copy()
    result['human_state'] = [STATES[i] for i in y]
    result['predicted_state'] = [STATES[i] for i in p.argmax(axis=1)]
    for i, state in enumerate(STATES):
        result['prob_' + state.lower().replace('-', '')] = p[:, i]
    return result


def animal_bootstrap(y, predicted, reference, groups, repetitions=2000):
    animals = sorted(set(groups))
    indices = {animal: np.flatnonzero(groups == animal) for animal in animals}
    rng = np.random.default_rng(SEED)
    values, differences = [], []
    for _ in range(repetitions):
        selected = np.concatenate([indices[a] for a in rng.choice(animals, len(animals), replace=True)])
        count = np.bincount(y[selected], minlength=4)
        if np.any(count == 0):
            continue
        correct = np.bincount(y[selected][predicted[selected] == y[selected]], minlength=4) / count
        baseline = np.bincount(y[selected][reference[selected] == y[selected]], minlength=4) / count
        values.append(float(correct.mean()))
        differences.append(float((correct - baseline).mean()))
    return {'balanced_accuracy_95_interval': np.percentile(values, [2.5, 97.5]).tolist(),
            'paired_difference_vs_morphometry_95_interval': np.percentile(differences, [2.5, 97.5]).tolist(),
            'bootstrap_unit': 'whole animal', 'replicates': len(values),
            'interpretation': 'Intervals resample the fixed out-of-fold predictions by animal; they do not measure model-training variability.'}


def train(args):
    import joblib
    from sklearn.metrics import confusion_matrix
    from threadpoolctl import threadpool_limits
    run = args.run_dir.resolve()
    protocol, manifest, x = load_matrix(run, reviewed_only=True)
    out = run / 'validation'
    out.mkdir(exist_ok=False)
    models = run / 'models'
    models.mkdir(exist_ok=False)
    labels = pd.read_csv(run / 'human_labels_used.csv').set_index('cell_uid')
    y = np.asarray([STATES.index(s) for s in labels.loc[manifest.cell_uid, 'state']], dtype=int)
    groups = manifest.animal_key.to_numpy(str)
    n_shape = len(protocol['shape_features'])
    candidates = protocol['configurations']
    # Write every outer/inner identity before calculating any performance.
    design = []
    for animal in sorted(set(groups)):
        outer_fit = np.flatnonzero(groups != animal)
        outer_test = np.flatnonzero(groups == animal)
        inner = grouped_inner_splits(y[outer_fit], groups[outer_fit])
        design.append({'held_out_animal': animal,
                       'training_uids': manifest.iloc[outer_fit].cell_uid.tolist(),
                       'test_uids': manifest.iloc[outer_test].cell_uid.tolist(),
                       'inner_folds': [{'fit_uids': manifest.iloc[outer_fit[fit]].cell_uid.tolist(),
                                        'selection_uids': manifest.iloc[outer_fit[validation]].cell_uid.tolist()} for fit, validation in inner]})
    atomic_json(out / 'validation_design.json', design)
    predicted = {name: np.full((len(y), 4), np.nan) for name in (*FAMILIES, 'nested_selected')}
    choices = []
    begin = time.monotonic()
    with threadpool_limits(limits=4):
        for fold, item in enumerate(design, 1):
            animal = item['held_out_animal']
            fit, test = np.flatnonzero(groups != animal), np.flatnonzero(groups == animal)
            best, best_by_family, selection = choose_config(x[fit], y[fit], groups[fit], n_shape, candidates)
            for family, config in best_by_family.items():
                model = make_model(config, n_shape, len(fit)).fit(x[fit], y[fit])
                predicted[family][test] = probabilities(model, x[test])
            predicted['nested_selected'][test] = predicted[best['family']][test]
            choices.append({'held_out_animal': animal, 'selected': best, 'inner_results': selection})
            atomic_json(out / 'model_selection_audit.json', choices)
            print(f'Outer animal {fold}/{len(design)}: {animal}; selected {best}; {time.monotonic()-begin:.1f}s', flush=True)
        source = Path(protocol['source_analysis'])
        template = pd.read_csv(source / 'microglia_annotation_template.csv').set_index('cell_uid')
        reference = np.asarray([STATES.index(s) for s in template.loc[manifest.cell_uid, 'pseudo_state']], dtype=int)
        predicted['morphometric_proposals'] = np.eye(4)[reference]
        metrics = {}
        for name, p in predicted.items():
            if not np.isfinite(p).all():
                raise ValueError(f'Incomplete out-of-fold predictions for {name}')
            metrics[name] = scores(y, p)
            metrics[name].update(animal_bootstrap(y, p.argmax(axis=1), reference, groups))
            prediction_frame(manifest, y, p).to_csv(out / f'{name}_out_of_fold.csv', index=False)
            pd.DataFrame(confusion_matrix(y, p.argmax(axis=1), labels=list(range(4))), index=STATES, columns=STATES).to_csv(out / f'{name}_confusion.csv')
        # Historical 40-cell comparison uses only the historical 163 training cells
        # for scaling, parameter/model choice, and fitting. It never picks the final model.
        legacy = pd.read_csv(run / 'legacy_development_split.csv').set_index('cell_uid')
        fit = np.flatnonzero(legacy.loc[manifest.cell_uid, 'split'].eq('training'))
        test = np.flatnonzero(legacy.loc[manifest.cell_uid, 'split'].eq('validation'))
        legacy_best, legacy_by_family, legacy_selection = choose_config(x[fit], y[fit], groups[fit], n_shape, candidates)
        legacy_metrics = {}
        legacy_prob = {}
        for family, config in legacy_by_family.items():
            model = make_model(config, n_shape, len(fit)).fit(x[fit], y[fit])
            p = probabilities(model, x[test])
            legacy_prob[family] = p
            legacy_metrics[family] = scores(y[test], p)
            prediction_frame(manifest.iloc[test].reset_index(drop=True), y[test], p).to_csv(out / f'{family}_legacy_40_cells.csv', index=False)
        legacy_metrics['nested_selected'] = scores(y[test], legacy_prob[legacy_best['family']])
        legacy_metrics['morphometric_proposals'] = scores(y[test], np.eye(4)[reference[test]])
        old_path = Path(protocol['review_dir']) / 'retrained_bncal_20260907_v3/per_microglia_classifier_predictions.csv'
        old = pd.read_csv(old_path).set_index('cell_uid')
        probability_columns = ['prob_' + s.lower().replace('-', '') for s in STATES]
        legacy_metrics['previous_tiny_cnn'] = scores(y[test], old.loc[manifest.iloc[test].cell_uid, probability_columns].to_numpy())
        best, by_family, selection = choose_config(x, y, groups, n_shape, candidates)
        for family, config in by_family.items():
            model = make_model(config, n_shape, len(y)).fit(x, y)
            bundle = {'model': model, 'configuration': config, 'states': STATES,
                      'shape_features': protocol['shape_features'], 'recipe_sha256': protocol['recipe_sha256'],
                      'encoder_weights_sha256': protocol['assets']['weights_sha256'],
                      'human_labels_sha256': protocol['human_labels_sha256'],
                      'training_uids': manifest.cell_uid.tolist(), 'training_animals': sorted(set(groups))}
            joblib.dump(bundle, models / f'{family}.joblib')
        shutil.copy2(models / f'{best["family"]}.joblib', models / 'selected_classifier.joblib')
    score_table = pd.DataFrame([{'method': name, 'balanced_accuracy': data['balanced_accuracy'], 'accuracy': data['accuracy'],
                                'ci_low': data['balanced_accuracy_95_interval'][0], 'ci_high': data['balanced_accuracy_95_interval'][1]} for name, data in metrics.items()])
    score_table.to_csv(out / 'nested_validation_summary.csv', index=False)
    pd.DataFrame([{'method': name, 'balanced_accuracy': data['balanced_accuracy'], 'accuracy': data['accuracy']} for name, data in legacy_metrics.items()]).to_csv(out / 'legacy_validation_summary.csv', index=False)
    report = {'status': 'trained_candidate', 'manual_labels_n': len(y), 'animals_n': len(set(groups)),
              'encoder': 'Frozen DINOv2 ViT-S/14, pretrained on external images; no microscopy labels or images update the encoder.',
              'input_views': RECIPE, 'selected_final_configuration': best, 'final_inner_selection': selection,
              'nested_validation': metrics, 'legacy_development_validation': legacy_metrics,
              'legacy_training_selected_configuration': legacy_best, 'legacy_inner_selection': legacy_selection,
              'primary_validation_method': protocol['validation'], 'interpretation': protocol['interpretation'],
              'selected_model_path': str(models / 'selected_classifier.joblib'), 'selected_model_sha256': sha256(models / 'selected_classifier.joblib'),
              'human_labels_sha256': protocol['human_labels_sha256'], 'protocol_sha256': sha256(run / 'protocol.json'),
              'validation_design_sha256': sha256(out / 'validation_design.json'), 'training_script_sha256': sha256(Path(__file__).resolve()),
              'elapsed_seconds': time.monotonic()-begin, 'publication_status': 'Candidate; requires scientific review before promotion.'}
    shutil.copy2(Path(__file__).resolve(), run / 'source_code/09_improve_microglia_classifier_trained.py')
    atomic_json(run / 'transfer_training_report.json', report)
    print(score_table.to_string(index=False), flush=True)
    print('Final candidate configuration:', best, flush=True)
    render_validation_report(run, report)


def render_validation_report(run, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = {'morphometric_proposals': 'Existing morphology proposals', 'morphology_logistic': 'Supervised shape: logistic',
             'morphology_forest': 'Supervised shape: forest', 'dino_logistic': 'DINO image features',
             'hybrid_logistic': 'DINO + shape', 'nested_selected': 'Model selected within each fold'}
    order = ['morphometric_proposals', *FAMILIES, 'nested_selected']
    summary = report['nested_validation']
    manual_n, animals_n = report['manual_labels_n'], report['animals_n']
    selected = summary['nested_selected']
    baseline = summary['morphometric_proposals']
    improvement = (selected['balanced_accuracy'] > baseline['balanced_accuracy']
                   and selected['paired_difference_vs_morphometry_95_interval'][0] > 0)
    decision = ('The candidate has a positive paired animal-bootstrap difference from existing morphology proposals and merits independent review.'
                if improvement else 'The tested candidates do not establish an improvement over existing morphology proposals in nested validation.')
    decision += ' The new models and full-cell predictions remain development candidates; this experiment does not modify Figure 5 publication outputs.'
    legacy = report['legacy_development_validation']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw={'width_ratios': [1.4, 1]})
    values = np.asarray([summary[name]['balanced_accuracy'] for name in order])
    intervals = np.asarray([summary[name]['balanced_accuracy_95_interval'] for name in order])
    axes[0].barh(np.arange(len(order)), values, color=['#a0a0a0', '#77aab7', '#5194a5', '#4278a5', '#315c85', '#235448'])
    axes[0].errorbar(values, np.arange(len(order)), xerr=np.maximum(0, np.vstack((values - intervals[:, 0], intervals[:, 1]-values))), fmt='none', color='#222222', capsize=3)
    axes[0].set_yticks(np.arange(len(order)), [names[name] for name in order])
    axes[0].invert_yaxis()
    axes[0].set(xlim=(0, 1), xlabel='Balanced accuracy (animal bootstrap 95% interval)', title=f'Nested validation across {animals_n} animals')
    for i, value in enumerate(values):
        axes[0].text(.99, i, f'{value:.1%}', va='center', ha='right', fontsize=9)
    cm = pd.read_csv(run / 'validation/nested_selected_confusion.csv', index_col=0).loc[list(STATES), list(STATES)].to_numpy()
    rates = cm / cm.sum(axis=1, keepdims=True)
    axes[1].imshow(rates, cmap='Blues', vmin=0, vmax=1)
    for i in range(4):
        for j in range(4):
            axes[1].text(j, i, f'{cm[i,j]}\n{rates[i,j]:.0%}', ha='center', va='center', fontsize=9, color='white' if rates[i,j] > .55 else '#142536')
    axes[1].set(xticks=range(4), yticks=range(4), xticklabels=STATES, yticklabels=STATES,
                xlabel='Predicted class', ylabel='Human-reviewed class', title='Animal-held-out errors')
    axes[1].tick_params(axis='x', rotation=25)
    fig.suptitle(f'Figure 5 classifier comparison: {manual_n} confirmed human reviews', fontsize=13)
    fig.text(.5, .025, 'Learned models exclude whole animals from fitting and parameter selection.\nReview sample selected by proposed class; independent validation is still required.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .10, 1, .94])
    fig.savefig(run / 'classifier_comparison.pdf')
    fig.savefig(run / 'classifier_comparison.png', dpi=300)
    plt.close(fig)
    lines = ['# Figure 5 classifier comparison', '',
             decision, '',
             f'Frozen pretrained DINOv2 features from a closer target crop and a wider branch-context crop are compared with mask morphology. Four rotations are averaged; masks identify the target. Only the {manual_n} human labels supervise learning.', '',
             f'Primary evaluation: nested leave-one-animal-out validation, with all normalization, PCA and model selection inside the training animals. The {animals_n} held-out animals provide {manual_n} out-of-fold predictions. This is internal validation, not a new independent test cohort.', '',
             '| Method | Balanced accuracy | Accuracy |', '|---|---:|---:|']
    for name in order:
        lines.append(f"| {names[name]} | {summary[name]['balanced_accuracy']:.1%} | {summary[name]['accuracy']:.1%} |")
    lines += ['', f"Final model selected using grouped CV on all reviewed cells: `{report['selected_final_configuration']}`.",
              f'The final fitted model uses all {manual_n} labels; its own fitted predictions are not validation results.', '',
              '[Comparison figure](classifier_comparison.pdf) · [Training report](transfer_training_report.json) · [Frozen protocol](protocol.json) · [Input crops](input_crop_examples.pdf)', '',
              'The historical 40-cell comparison is reported separately in validation/legacy_validation_summary.csv because those cells have already guided earlier model development. No prior 40-cell score was used to select this final model.', '',
              f"On those same 40 cells: previous CNN {legacy['previous_tiny_cnn']['balanced_accuracy']:.1%}, DINO image features {legacy['dino_logistic']['balanced_accuracy']:.1%}, and the selected candidate {legacy['nested_selected']['balanced_accuracy']:.1%} balanced accuracy. The broader animal-held-out results above determine the interpretation.", '',
              f'The audit writes held_out_errors.csv, held_out_error_examples.pdf, next_review_queue_300.csv and improvement_assessment.json. The new review queue contains no human labels, excludes the {manual_n} reviewed cells, balances 20 targets per animal, and prioritizes uncertainty and disagreement. It is an acquisition sample, not an independent test set.', '',
              'Implementation and reproducible commands: ../../../../Fig5/README.txt. The archived source_code directory preserves the exact training implementation. Official pretrained encoder source: https://github.com/facebookresearch/dinov2 (commit and weight hashes are frozen in protocol.json).', '',
              'Raw microscopy, accepted anatomy, previous masks and human review CSVs are preserved. Figure outputs derived from these proposals remain candidates for scientific review.', '']
    (run / 'README.md').write_text('\n'.join(lines))


def predict(args):
    import joblib
    from threadpoolctl import threadpool_limits
    run = args.run_dir.resolve()
    protocol, manifest, x = load_matrix(run, reviewed_only=False)
    report = json.loads((run / 'transfer_training_report.json').read_text())
    if sha256(report['selected_model_path']) != report['selected_model_sha256']:
        raise ValueError('Selected classifier changed after fitting.')
    bundle = joblib.load(report['selected_model_path'])
    if bundle['recipe_sha256'] != protocol['recipe_sha256'] or bundle['human_labels_sha256'] != protocol['human_labels_sha256']:
        raise ValueError('Model preprocessing or label provenance differs.')
    with threadpool_limits(limits=4):
        p = probabilities(bundle['model'], x)
    result = manifest[['cell_uid']].copy()
    result['microglia_state'] = [STATES[i] for i in p.argmax(axis=1)]
    family = bundle['configuration']['family']
    result['classifier_source'] = 'human_supervised_' + family
    result['classifier_confidence'] = p.max(axis=1)
    result['classifier_uncertain'] = p.max(axis=1) < RECIPE['uncertain_score_below']
    for i, state in enumerate(STATES):
        result['prob_' + state.lower().replace('-', '')] = p[:, i]
    result.to_csv(run / 'per_microglia_classifier_predictions.csv', index=False)
    atomic_json(run / 'prediction_receipt.json', {'cells_n': len(result), 'selected_model_sha256': report['selected_model_sha256'],
                'predictions_sha256': sha256(run / 'per_microglia_classifier_predictions.csv'),
                'dino_features_sha256': sha256(run / 'dino_features.npy'), 'family': family,
                'manual_labels_n': protocol['reviewed_cells_n'], 'uncertain_cells_n': int(result.classifier_uncertain.sum()),
                'uncertain_score_threshold': RECIPE['uncertain_score_below'],
                'class_counts': result.microglia_state.value_counts().to_dict(),
                'interpretation': 'All source cells classified by the refitted candidate, including training cells. These counts are model outputs, not validated population prevalence.'})
    print(f'Classified {len(result)} cells with {family}.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'extract', 'train', 'predict'])
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--review-dir', type=Path, default=DEFAULT_REVIEW)
    parser.add_argument('--assets-dir', type=Path, default=DEFAULT_ASSETS)
    parser.add_argument('--all-cells', action='store_true')
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--cells-per-batch', type=int, default=8)
    parser.add_argument('--image-batch-size', type=int, default=32)
    args = parser.parse_args()
    if min(args.cells_per_batch, args.image_batch_size) < 1:
        parser.error('Batch sizes must be positive.')
    {'prepare': prepare, 'extract': extract, 'train': train, 'predict': predict}[args.stage](args)


if __name__ == '__main__':
    main()
