#!/usr/bin/env python3
"""Install verified spatial replacements while retaining nonspatial panel files."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import pandas as pd
from PIL import Image, ImageChops


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--analysis',type=Path,required=True)
    args=parser.parse_args()
    root=args.root.resolve(); analysis=args.analysis.resolve()
    checks=[]
    for figure,kept in [('Fig3','ABCDEF'),('Fig4','ABCDEFGH')]:
        stage=analysis/'assembled'/figure; current=root/figure
        for source in sorted((stage/'panels').glob('*.png')):
            match=re.search(r'(?:Panel|panel)_([A-Z])(?:_|\.)',source.name)
            if match and match.group(1) in kept:
                target=current/'panels'/source.name
                if not target.exists():raise FileNotFoundError(target)
                with Image.open(source) as a,Image.open(target) as b:
                    same=a.size==b.size and ImageChops.difference(a.convert('RGB'),b.convert('RGB')).getbbox() is None
                checks.append(dict(panel=str(target.relative_to(root)),identical_pixels=same,
                                   original_sha256=digest(target),staged_sha256=digest(source)))
        if any(not item['identical_pixels'] for item in checks):
            (analysis/'nonspatial_panel_comparison.json').write_text(json.dumps(checks,indent=2)+'\n')
            raise RuntimeError('An unchanged panel differs from its previous raster; inspect comparison before installing')
    # Comparison of all unchanged panels completes before any canonical write.
    for figure,cohort,kept in [('Fig3','NPY','ABCDEF'),('Fig4','POMC','ABCDEFGH')]:
        stage=analysis/'assembled'/figure; current=root/figure
        for source in stage.iterdir():
            if source.is_file() and (source.suffix in ('.pdf','.png') or 'caption' in source.name):
                shutil.copy2(source,current/source.name)
        for folder in ('panels','legends','provenance'):
            (current/folder).mkdir(exist_ok=True)
            for source in (stage/folder).glob('*'):
                if not source.is_file():continue
                match=re.search(r'(?:Panel|panel)_([A-Z])(?:_|\.)',source.name)
                if folder=='panels' and match and match.group(1) in kept:continue
                shutil.copy2(source,current/folder/source.name)
        archive=analysis/'archive'/figure
        for source in list((current/'source_data').glob('*'))+list((current/'provenance').glob('*')):
            if source.is_file() and ('radial_occurrence' in source.name or 'spatial_exact_' in source.name or 'spatial_ring_verification' in source.name):
                destination=archive/source.parent.name/source.name
                destination.parent.mkdir(parents=True,exist_ok=True)
                if destination.exists() and digest(destination)!=digest(source):raise RuntimeError('Archive collision: '+str(destination))
                shutil.move(str(source),str(destination))
        if figure=='Fig4':
            # Keep the previously published master-caption aliases current.
            for name in ('Figure_ARC_ME_multipanel_LEGEND.txt','Figure_ARC_ME_multipanel_caption.txt'):
                alias=current/'legends'/name
                if alias.exists():
                    previous=archive/'legends'/name
                    previous.parent.mkdir(parents=True,exist_ok=True)
                    if not previous.exists():shutil.copy2(alias,previous)
                    shutil.copy2(current/'legends'/'Figure_4_LEGEND.txt',alias)
        prefix='Figure_3_cFos_NPY' if figure=='Fig3' else 'Figure_4'
        for name in ('animal_metrics','unit_metrics','condition_summary','treatment_contrasts','stratum_metrics','cluster_counts_by_condition','contour_manifest','bandwidth_sensitivity_condition_summary'):
            path=analysis/(name+'.csv')
            if path.exists():
                frame=pd.read_csv(path)
                frame.loc[frame.cohort.eq(cohort)].to_csv(current/'source_data'/f'{prefix}_spatial_shape_{name}.csv',index=False)
        members=pd.read_csv(analysis/'cluster_membership.csv.gz')
        members.loc[members.cohort.eq(cohort)].to_csv(current/'source_data'/f'{prefix}_spatial_shape_cluster_membership.csv.gz',index=False,compression='gzip')
        for name in ('analysis_metadata.json','validation.json'):
            shutil.copy2(analysis/name,current/'provenance'/('spatial_shape_'+name))
        if (analysis/'README.md').exists():shutil.copy2(analysis/'README.md',current/'source_data'/f'{prefix}_spatial_shape_README.md')
    for item in checks:
        item['original_file_preserved']=digest(root/item['panel'])==item['original_sha256']
    if not all(item['original_file_preserved'] for item in checks):raise AssertionError('Unchanged panel bytes modified')
    (analysis/'nonspatial_panel_comparison.json').write_text(json.dumps(checks,indent=2)+'\n')
    record={'status':'PASS','unchanged_panel_rasters':len(checks),'unchanged_panel_files_preserved':True,
            'main_figures':{str(p.relative_to(root)):digest(p) for p in [root/'Fig3/Figure_3_cFos_NPY.pdf',root/'Fig3/Figure_3_cFos_NPY.png',root/'Fig4/Figure_4.pdf',root/'Fig4/Figure_4.png']}}
    (analysis/'publication_receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
