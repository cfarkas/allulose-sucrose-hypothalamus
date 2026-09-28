#!/usr/bin/env python3
"""Reconstruct Nancy's missing sucrose cartoon from immutable native masks.

Uses exactly the accepted Figure 3 native overlap classifier and display
transform. A ventricular boundary is drawn only from the saved human review.
Run Fig3/02b_review_sucrose_ventricle.py to delineate the native registered DAPI.
"""
from pathlib import Path
import importlib.util
import json
import sys
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root",type=Path,default=Path.cwd())
parser.add_argument("--output-dir",type=Path)
args=parser.parse_args()
ROOT=args.root.resolve()
OUT=args.output_dir or ROOT/"analyses/Fig3/results/complete_conditions_20260908/sucrose_panels"
OUT.mkdir(parents=True, exist_ok=True)

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m

m = module('npy_original_renderer', ROOT / 'Fig3/03_make_figure_3_cfos_npy.py')
c = module('npy_original_cartoons', ROOT / 'Fig3/02_make_cfos_npy_cartoons.py')
spec = next(s for s in m.MICROSCOPY_SOURCES if s['condition'] == 'Sucrose') if hasattr(m, 'MICROSCOPY_SOURCES') else None
if spec is None:
    for value in vars(m).values():
        if isinstance(value, tuple) and value and isinstance(value[0], dict) and 'condition' in value[0]:
            spec = next((s for s in value if s.get('condition') == 'Sucrose' and 'mask_paths' in s), None)
            if spec: break
assert spec is not None
base = ROOT / 'Fig3/raw/legacy_experiment_2025_07_28'
labels = {}
sources = []
for marker, name in spec['mask_paths'].items():
    path = base / name
    labels[marker], _ = m._load_cellpose_labels(path, name)
    sources.append({'path':str(path.relative_to(ROOT)), 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
cfos_ids, npy_ids, dual_ids, observed = m._analysis_positive_sets(labels['dapi'], labels['cfos'], labels['npy'])
for key in ('dapi_nuclei', 'cfos_positive_nuclei', 'npy_positive_nuclei', 'dual_positive_nuclei'):
    assert observed[key] == spec['expected_analysis_counts'][key], (key, observed)
display = m._warp_labels_to_display(labels['dapi'], spec)
hil = module('sucrose_hil_validation', ROOT / 'Fig3/02b_review_sucrose_ventricle.py')
accepted = hil.load_accepted(ROOT, ROOT / hil.DEFAULT_REVIEW, display.shape)
if accepted is None:
    raise RuntimeError('Delineate and save the sucrose ventricle in HIL before rendering: http://127.0.0.1:34248/')
ventricle_mask, ventricle_receipt = accepted
# Same DAPI-supported display envelope used for Water and Allulose.
# It is a presentation aid, never an analysis ROI or denominator.
tissue = c._tissue_envelope(display)
tissue_outline = np.asarray(c.binary_fill_holes(tissue), dtype=bool)
assert c._contour_path_counts(tissue_outline) == {'external': 1, 'tree': 1}
visible = set(np.unique(display)) - {0}
counts = {'NPY':len(visible & set(npy_ids)), 'c-FOS':len(visible & set(cfos_ids)), 'dual':len(visible & set(dual_ids))}
colors = {'NPY':'#00e85e', 'c-FOS':'#ff37d4', 'dual':'#f2a900'}
for language in ('en', 'es'):
    fig, axes = plt.subplots(1,3, figsize=(7.44,2.55), facecolor='white')
    fig.subplots_adjust(left=.025,right=.988,bottom=.12,top=.735,wspace=.045)
    title = 'Sucrose — marker-positive DAPI nuclei' if language == 'en' else 'Sacarosa — núcleos DAPI positivos para marcador'
    fig.text(.5,.955,title,fontsize=12.5,fontweight='bold',ha='center',va='top')
    titles = ['NPY-positive nuclei','c-FOS-positive nuclei','Combined status'] if language == 'en' else ['Núcleos NPY positivos','Núcleos c-FOS positivos','Estado combinado']
    layers = [[(npy_ids,colors['NPY'])],[(cfos_ids,colors['c-FOS'])],[(set(npy_ids)-set(dual_ids),colors['NPY']),(set(cfos_ids)-set(dual_ids),colors['c-FOS']),(dual_ids,colors['dual'])]]
    for ax,ttl,ls in zip(axes,titles,layers):
        c._draw_background(ax, tissue, ventricle_mask)
        rgba = np.zeros((*display.shape,4),dtype=np.float32)
        rgba[:,:,:3] = matplotlib.colors.to_rgb(c.NEGATIVE_NUCLEUS_COLOR)
        rgba[:,:,3] = (display > 0).astype(np.float32) * .25
        ax.imshow(rgba,interpolation='nearest')
        for ids,col in ls:
            layer = np.zeros((*display.shape,4),dtype=np.float32)
            layer[:,:,:3] = matplotlib.colors.to_rgb(col)
            layer[:,:,3] = np.isin(display, list(ids))*.94
            ax.imshow(layer,interpolation='nearest')
        # The human contour annotates the image; underlying cell calls remain intact.
        ax.contour(ventricle_mask.astype(float), levels=[.5], colors='black', linewidths=1.25)
        ax.set_title(ttl,fontsize=10.5,fontweight='bold',color='#172033',pad=3)
        c._matched_add_scalebar(ax,display.shape,spec['um_per_px'])
        ax.set_axis_off()
    names = ['Other DAPI nuclei','NPY-positive nucleus','c-FOS-positive nucleus','Dual-positive nucleus'] if language == 'en' else ['Otros núcleos DAPI','Núcleo NPY positivo','Núcleo c-FOS positivo','Núcleo doble positivo']
    fig.legend([Patch(facecolor=x) for x in ['#d1e0f5',*colors.values()]],names,loc='lower center',bbox_to_anchor=(.5,.002),ncol=4,frameon=False,fontsize=7,handlelength=1.4,columnspacing=1.4)
    stem=OUT/f'Sucrose_NPY_cartoon_{language}'
    fig.savefig(stem.with_suffix('.pdf'),dpi=600)
    fig.savefig(stem.with_suffix('.png'),dpi=600)
    plt.close(fig)
receipt={'native_counts':observed,'display_counts':counts,'sources':sources,'display_transform':spec['display_to_raw'],'pixel_size_um':spec['um_per_px'],'ventricle_hil_receipt':str((ROOT / hil.DEFAULT_REVIEW / 'sucrose_ventricle_receipt.json').relative_to(ROOT)), 'ventricle_hil_receipt_sha256':hil.sha(ROOT / hil.DEFAULT_REVIEW / 'sucrose_ventricle_receipt.json'), 'reviewer':ventricle_receipt['reviewer'], 'tissue_display_envelope':{'method':'Original Figure 3 _tissue_envelope and _draw_background; DAPI density; display only', 'area_px':int(tissue.sum()), 'outline_area_px':int(tissue_outline.sum()), 'outline_sha256':hashlib.sha256(np.ascontiguousarray(tissue_outline, dtype=np.uint8).tobytes()).hexdigest(), 'contour_paths':c._contour_path_counts(tissue_outline)}, 'note':'Native DAPI labels with the human-reviewed black ventricular boundary and DAPI-supported external tissue outline. Display only; no cell classification or animal statistic changed.'}
(OUT/'receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
print(json.dumps(receipt,ensure_ascii=False))
