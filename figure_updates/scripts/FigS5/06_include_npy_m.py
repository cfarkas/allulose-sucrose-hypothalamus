#!/usr/bin/env python3
"""Apply the author's NPY-M inclusion and channel clarification to copied manifests.

Source exports, microscopy and native masks are never edited. The overlay
includes the eight WT acquisitions plus NPY-M (male NPY-transgenic Water).
"""
from pathlib import Path
import argparse,json,hashlib
import pandas as pd

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path.cwd());ap.add_argument('--output-dir',type=Path);args=ap.parse_args()
    root=args.root.resolve();source=root/'FigS5/channel_tiffs';out=args.output_dir or root/'FigS5/provenance/inclusion_NPY_M_20260908';out.mkdir(parents=True,exist_ok=True)
    table=pd.read_csv(source/'sample_manifest.csv');channels=pd.read_csv(source/'channel_map.csv');acq='2026-07-10_NPY-M_Water'
    selected=table['include_main_analysis'].astype(str).str.lower().isin(['true','1']) | table.acquisition_id.eq(acq)
    table=table.loc[selected].copy();assert len(table)==9 and table.acquisition_id.nunique()==9
    row=table.acquisition_id.eq(acq);assert row.sum()==1
    table.loc[row,'include_main_analysis']=True;table.loc[row,'exclusion_reason']=''
    table.loc[row,'assay_hint']='cFOS_ACTH_CLIP';table.loc[row,'channel_map_status']='CONFIRMED_BY_AUTHOR_20260908'
    table.loc[row,'provisional_marker_map']=json.dumps({'DAPI':'DAPI','AF488':'NPY','AF546':'ACTH_CLIP','AF633':'cFOS'})
    table.loc[row,'notes']='Male; author explicitly requested inclusion and confirmed DAPI/NPY-GFP488/ACTH-CLIP546/cFOS633 mapping on 8 September 2026. Catalogue ACTH/CLIP F-3 sc-373878, 1:100.'
    channels=channels[channels.acquisition_id.isin(table.acquisition_id)].copy()
    for idx,r in channels[channels.acquisition_id.eq(acq)].iterrows():
        # Fluor labels are recorded in the exporter; marker identities are author-confirmed.
        fluor=r['fluor'];marker={'DAPI':'DAPI','AF488':'NPY','AF546':'ACTH_CLIP','AF633':'cFOS'}[fluor]
        channels.loc[idx,'marker']=marker
        for col in ['mapping_status','channel_map_status','status']:
            if col in channels:channels.loc[idx,col]='CONFIRMED_BY_AUTHOR_20260908'
        if 'mapping_rule' in channels:channels.loc[idx,'mapping_rule']='AUTHOR_CONFIRMED_20260908'
        if 'review_notes' in channels:channels.loc[idx,'review_notes']='Author confirmed 488 NPY-GFP, 546 ACTH/CLIP and 633 c-FOS; include NPY-M in the Water group.'
        if 'mapping_source' in channels:channels.loc[idx,'mapping_source']='Author confirmation; NPY-M source filename POMC_TRITC_CFOS_647'
        if 'use_for_analysis' in channels:channels.loc[idx,'use_for_analysis']=marker!='NPY'
        if 'notes' in channels:channels.loc[idx,'notes']='Author-confirmed four-channel mapping; NPY fluorescence is not the c-FOS channel in this animal.'
    table.to_csv(out/'sample_manifest.csv',index=False);channels.to_csv(out/'channel_map.csv',index=False)
    receipt={'schema':'figs5_NPY_M_author_inclusion_v1','date':'2026-09-08','author_instruction':'Include NPY-M','author_confirmed_channel_assignment':{'DAPI':'DAPI','AF488':'NPY-GFP','AF546':'ACTH/CLIP','AF633':'c-FOS'},'genotype':'NPY-Tg','sex':'male','condition':'Water','cohort_counts':table.condition.value_counts().to_dict(),'previous_exclusion':'NON_WT_COHORT_REQUIRES_SEPARATE_ANALYSIS','source_sha256':{str(p.relative_to(root)):sha(p) for p in [source/'sample_manifest.csv',source/'channel_map.csv']},'overlay_sha256':{p.name:sha(p) for p in [out/'sample_manifest.csv',out/'channel_map.csv']},'raw_data_unchanged':True}
    (out/'INCLUSION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'cohort':receipt['cohort_counts'],'overlays':str(out)}))
if __name__=='__main__':main()
