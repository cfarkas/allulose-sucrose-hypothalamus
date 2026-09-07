#!/usr/bin/env python3
"""Create versioned Iba1-driven cell/soma mask proposals and reviewer inputs.

Seeds come from locally prominent, thick Iba1 signal. DAPI is used only to
check association after a soma candidate exists. These are machine proposals,
not accepted human masks; the existing raw data and accepted masks are read only.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import logging
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu, sobel
from skimage.measure import regionprops
from skimage.morphology import binary_closing, disk, remove_small_holes, remove_small_objects
from skimage.segmentation import watershed
from PIL import Image

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'Fig5/01_analyze_gfap_iba1_microglia.py').is_file())
PARAMETERS = dict(background_width_px=31, smooth_sigma=.6, soma_opening_width_px=3,
                  soma_score_sigma=.7, minimum_peak_spacing_px=5, minimum_soma_score=.22,
                  minimum_peak_intensity=.32, minimum_local_contrast=.12,
                  seed_noise_multiplier=3., max_nucleus_distance_px=5,
                  process_radius_px=32, soma_radius_px=7, min_cell_area_px=8)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def annulus(signal, y, x):
    y0,y1=max(0,y-14),min(signal.shape[0],y+15)
    x0,x1=max(0,x-14),min(signal.shape[1],x+15)
    yy,xx=np.indices((y1-y0,x1-x0));rr=np.hypot(yy+y0-y,xx+x0-x)
    values=signal[y0:y1,x0:x1][(rr>=8)&(rr<=14)]
    background=float(np.median(values)) if values.size else 0.
    noise=1.4826*float(np.median(abs(values-background))) if values.size else 0.
    return background,noise


def segment(signal, dapi):
    signal=np.asarray(signal,dtype=np.float32)
    if signal.ndim!=2 or signal.shape!=dapi.shape:
        raise ValueError('Iba1 and DAPI must have the same 2-D shape.')
    if not np.isfinite(signal).all():raise ValueError('Nonfinite Iba1 signal.')
    corr=np.maximum(0,signal-ndi.grey_opening(signal,size=(31,31)))
    corr=np.clip(corr/max(float(np.percentile(corr,99.8)),1),0,1)
    empty=np.zeros(signal.shape,np.int32)
    if not np.any(corr>0) or not np.any(dapi>0):
        return empty,empty.copy(),corr,[]
    smooth=ndi.gaussian_filter(corr,.6)
    threshold=float(threshold_otsu(corr[corr>0]))*.65
    fg=remove_small_objects(remove_small_holes(binary_closing(smooth>threshold,disk(1)),area_threshold=12),min_size=6)
    score=ndi.gaussian_filter(ndi.grey_opening(corr,size=(3,3)),.7)
    peaks=peak_local_max(score,min_distance=5,threshold_abs=max(threshold*1.5,.22),exclude_border=5)
    distance,nearest=ndi.distance_transform_edt(dapi==0,return_indices=True)
    markers=empty.copy();taken=set();seeds=[]
    for y,x in sorted(peaks,key=lambda p:-score[tuple(p)]):
        if not fg[y,x] or distance[y,x]>5:continue
        background,noise=annulus(smooth,y,x)
        if smooth[y,x]<max(.32,background+max(.12,3*noise)):continue
        nucleus=int(dapi[tuple(nearest[:,y,x])])
        if nucleus in taken:continue
        taken.add(nucleus);identity=len(seeds)+1;markers[y,x]=identity
        seeds.append(dict(label=identity,y=int(y),x=int(x),dapi_nucleus_label=nucleus,
                          local_background=background,local_noise=noise,peak=float(smooth[y,x])))
    if not seeds:return empty,empty.copy(),corr,[]
    distance,nearest=ndi.distance_transform_edt(markers==0,return_indices=True)
    near_marker=markers[nearest[0],nearest[1]]
    levels=np.zeros(len(seeds)+1)
    for seed in seeds:
        levels[seed['label']]=max(threshold,seed['local_background']+1.5*seed['local_noise'],seed['peak']*.15)
    fg &= (distance<=32)&(smooth>levels[near_marker])
    fg=remove_small_holes(binary_closing(fg,disk(1)),area_threshold=12)
    labels=watershed(sobel(smooth),markers,mask=fg).astype(np.int32)
    somas=empty.copy()
    for rp in regionprops(labels):
        if rp.area<8:
            labels[labels==rp.label]=0
            continue
        seed=seeds[rp.label-1];sy,sx=seed['y'],seed['x'];sl=rp.slice
        yy,xx=np.indices(rp.image.shape);radius=np.hypot(yy+sl[0].start-sy,xx+sl[1].start-sx)
        core=rp.image&(smooth[sl]>max(threshold,seed['peak']*.38))&(radius<7)
        core=ndi.binary_fill_holes(binary_closing(core,disk(1)))&rp.image
        components,_=ndi.label(core);value=int(components[sy-sl[0].start,sx-sl[1].start])
        core=components==value if value else np.zeros_like(core)
        if not core.any():
            labels[labels==rp.label]=0
            continue
        somas[sl][core]=rp.label
    return labels,somas,corr,seeds


def crop(array, center, size=96):
    y0,x0=int(round(center[0]))-size//2,int(round(center[1]))-size//2
    y1,x1=y0+size,x0+size
    patch=np.zeros((size,size),array.dtype)
    sy0,sx0=max(0,y0),max(0,x0);sy1,sx1=min(array.shape[0],y1),min(array.shape[1],x1)
    patch[sy0-y0:sy1-y0,sx0-x0:sx1-x0]=array[sy0:sy1,sx0:sx1]
    return patch,(y0,x0)


def analyzer():
    path=PAPER/'Fig5/01_analyze_gfap_iba1_microglia.py'
    spec=importlib.util.spec_from_file_location('fig5_refinement_analyzer',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    module.import_science_stack()
    return module


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--limit',type=int,default=0)
    cli=ap.parse_args();out=cli.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    module=analyzer();args=module.build_argparser_gfap_iba1().parse_args([])
    args.manual_region_dir=str(module.DEFAULT_MANUAL_REGION_ROOT)
    args.require_manual_regions=True;args.save_reconstructed_tifs=True
    args.prepare_labels_only=True;args.save_microglia_patches=True;args.save_qc=False
    logger=module.setup_logger(out,False)
    paths=sorted((PAPER/'Fig5/raw_data').glob('*/S*/Image_*_Iba1.tif'))
    if cli.limit:paths=paths[:cli.limit]
    elif len(paths)!=61:raise ValueError(f'Expected 61 sections, found {len(paths)}')
    all_rows=[];section_audit=[];all_regions=[];section_summaries=[]
    for index,path in enumerate(paths,1):
        started=time.monotonic();sample=path.parent.parent.name;section=path.parent.name
        dapi_path=PAPER/'Fig5/hil_review/dapi_final_20260830_v1'/sample/section/f'{sample}_{section}_DAPI_seg.npy'
        dapi=module.ensure_label_image(module.load_seg_mask(dapi_path))
        signal=module._gi_read_signal(path,'Iba1')
        labels,somas,corr,seeds=segment(signal,dapi)
        target=out/'mask_proposals'/sample/section;target.mkdir(parents=True)
        mask_path=target/f'{sample}_{section}_Iba1_seg.npy'
        np.save(mask_path,{'masks':labels,'source':'Iba1_signal_soma_seeded','human_accepted':False,'parameters':PARAMETERS},allow_pickle=True)
        tifffile.imwrite(target/'cell_labels.tif',labels.astype(np.uint32),compression='deflate')
        tifffile.imwrite(target/'soma_labels.tif',somas.astype(np.uint32),compression='deflate')
        images={name:path.with_name(path.name.replace('_Iba1.',f'_{name}.')) for name in ['Iba1','DAPI','GFAP']}
        sec=module.GISection(sample=sample,sample_dir=path.parent.parent,section_index=int(section[1:]),section=section,
                            section_dir=path.parent,images=images,segs={'DAPI':dapi_path,'Iba1':mask_path})
        rows,regions,summary=module.process_gfap_iba1_section(sec,args,out,logger)
        if not rows or summary.get('status')!='ok':raise ValueError(f'Section extraction failed: {sample}/{section}: {summary.get("status")}')
        soma_props={r.label:r for r in regionprops(somas)}
        for row in rows:
            identity=int(row['iba1_label']);body=soma_props.get(identity)
            row['segmentation_version']='iba1_signal_soma_v1'
            row['mask_human_accepted']=False
            if body is None:
                row['cell_qc_pass']=False
                row['cell_qc_rejection_reasons'] += ';missing_soma'
                continue
            center=body.centroid
            uid=row['cell_uid'].replace('__iba1_','__iba1sig_')
            row['cell_uid']=uid
            row.update(soma_centroid_y=float(center[0]),soma_centroid_x=float(center[1]),soma_area_px=float(body.area),
                       soma_fraction=float(body.area/max(row['area_px'],1)),process_fraction=float(1-body.area/max(row['area_px'],1)),
                       source_height=signal.shape[0],source_width=signal.shape[1])
            if row['cell_qc_pass']:
                patch=module._gi_extract_patch(corr,labels==identity,center,96)
                patch_path=out/'microglia_patches'/sample/section/f'{uid}.png'
                module._gi_save_patch(patch_path,patch);row['patch_path']=str(patch_path)
                soma_patch,origin=crop(somas==identity,center)
                soma_path=out/'soma_patches'/sample/section/f'{uid}.png';soma_path.parent.mkdir(parents=True,exist_ok=True)
                Image.fromarray(soma_patch.astype(np.uint8)*255).save(soma_path)
                row['soma_patch_path']=str(soma_path);row['patch_origin_y'],row['patch_origin_x']=origin
                dapi_patch,_=crop(dapi==int(row['dapi_nucleus_label']),center)
                dapi_out=out/'nucleus_patches'/sample/section/f'{uid}.png';dapi_out.parent.mkdir(parents=True,exist_ok=True)
                Image.fromarray(dapi_patch.astype(np.uint8)*255).save(dapi_out);row['nucleus_patch_path']=str(dapi_out)
            all_rows.append(row)
        all_regions.extend(regions);section_summaries.append(summary)
        receipt={'sample':sample,'section':section,'status':'PROPOSED_FOR_HUMAN_REVIEW',
                 'algorithm':'Iba1 fluorescence soma seeds, DAPI association check, signal-constrained watershed',
                 'parameters':PARAMETERS,'source_iba1':str(path),'source_iba1_sha256':digest(path),
                 'source_dapi_labels':str(dapi_path),'source_dapi_labels_sha256':digest(dapi_path),
                 'cell_masks_sha256':digest(target/'cell_labels.tif'),'soma_masks_sha256':digest(target/'soma_labels.tif'),
                 'candidates':len(regionprops(labels)),'qc_accepted':sum(bool(r['cell_qc_pass']) for r in rows),
                 'seconds':time.monotonic()-started,'raw_data_modified':False}
        (target/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        section_audit.append(receipt)
        print(f'[{index}/{len(paths)}] {sample}/{section}: {receipt["candidates"]} soma-seeded cells, {receipt["qc_accepted"]} pass QC ({receipt["seconds"]:.1f}s)',flush=True)
    candidates=pd.DataFrame(all_rows);candidates.to_csv(out/'microglia_object_qc.csv',index=False)
    good=candidates.loc[candidates.cell_qc_pass].copy()
    classified,meta=module.classify_microglia(good,args,out,logger)
    classified.to_csv(out/'per_microglia_cell_measurements.csv',index=False)
    pd.DataFrame(all_regions).to_csv(out/'per_section_region_summary.csv',index=False)
    pd.DataFrame(section_summaries).to_csv(out/'per_section_summary.csv',index=False)
    pd.DataFrame(section_audit).to_csv(out/'mask_refinement_audit.csv',index=False)
    report={'status':'MASK_PROPOSALS_READY_FOR_REVIEW','sections':len(paths),'cells':len(classified),
            'source_analyzer_sha256':digest(PAPER/'Fig5/01_analyze_gfap_iba1_microglia.py'),
            'refinement_script_sha256':digest(Path(__file__)),'parameters':PARAMETERS,
            'manual_labels_n':0,'cnn_used':False,'raw_data_modified':False,
            'previous_accepted_masks_modified':False,'class_counts':classified.microglia_state.value_counts().to_dict()}
    (out/'mask_refinement_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
