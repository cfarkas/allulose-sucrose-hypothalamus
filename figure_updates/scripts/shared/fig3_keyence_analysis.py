#!/usr/bin/env python3
"""Extend Figure3 with three cohort-specific Keyence animals; retain source masks."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,importlib.util,json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from skimage.measure import regionprops_table
import shape_spatial_analysis as shape
PAPER=Path(__file__).resolve().parents[2]
SOURCE=PAPER/'thesis_and_manuscript/Scientific_Reports_Segura_et_al_2026/spatial_reanalysis_20260919'
DEFAULT=PAPER/'analyses/fig3_keyence_extension_20260920'
CONDITIONS=('Water','Sucrose','Allulose')

def load_legacy():
 spec=importlib.util.spec_from_file_location('fig3_accepted_association',PAPER/'Fig3/01_analyze_cfos_npy.py')
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m

def load_fixed_keyence_inputs(out):
 """Verify frozen inputs without consulting the mutable Figure 4 result tree."""
 out=Path(out).resolve();manifest_path=out/'source_records/fixed_keyence_inputs.json'
 if not manifest_path.is_file():
  raise FileNotFoundError('Required frozen Keyence input manifest is missing: '+str(manifest_path))
 frozen=json.loads(manifest_path.read_text())
 if frozen.get('schema')!='fig3-keyence-fixed-inputs-v1':raise ValueError('Unsupported fixed Keyence input schema')
 for key,name in [('source_section_manifest','keyence_accepted_canonical_sections.csv'),('source_mask_audit','nonpomc_masks_before_vs_after_size_qc.json')]:
  path=out/'source_records'/name
  if shape.sha256(path)!=frozen[key]['sha256']:raise ValueError('Frozen Keyence provenance record changed: '+str(path))
 selected=pd.read_csv(out/'source_records/keyence_accepted_canonical_sections.csv')
 audit=json.loads((out/'source_records/nonpomc_masks_before_vs_after_size_qc.json').read_text())
 audited={(r['sample'],r['section'],r['channel']):r['current_sha256'] for r in audit['files']}
 entries={(r['sample'],r['section']):r for r in frozen['sections']}
 expected=set(zip(selected['sample'],selected['section']))
 if len(entries)!=len(frozen['sections']) or len(expected)!=len(selected) or set(entries)!=expected:
  raise ValueError('Frozen Keyence section roster differs from accepted canonical sections')
 fixed_root=(out/'inputs/fixed_keyence_rois').resolve()
 def checked_path(record,original):
  path=(out/record['relative_path']).resolve()
  if not path.is_relative_to(fixed_root):raise ValueError('Frozen Keyence input is outside fixed input directory: '+str(path))
  if record['original_source_path']!=str(original):raise ValueError('Frozen Keyence source attribution changed: '+str(original))
  if not path.is_file() or shape.sha256(path)!=record['sha256']:
   raise ValueError('Frozen Keyence file is missing or its SHA256 changed: '+str(path))
  return str(path)
 rows=[]
 for row in selected.to_dict('records'):
  entry=entries[(row['sample'],row['section'])]
  if set(entry['masks'])!={'DAPI','NPY','cFOS'}:raise ValueError('Frozen Keyence section does not contain exactly three selected channels')
  row['_fixed_mask_paths']={}
  for channel in ('DAPI','NPY','cFOS'):
   record=entry['masks'][channel]
   if record['sha256']!=audited[(row['sample'],row['section'],channel)]:
    raise ValueError('Frozen Keyence mask differs from the pre-extension audit: '+row['sample']+' '+row['section']+' '+channel)
   row['_fixed_mask_paths'][channel]=checked_path(record,row[channel+'_labels'])
  row['_fixed_registration_path']=checked_path(entry['registration'],row['current_registration'])
  calibration=entry['calibration'];native=frozen['native_calibration_source']['frozen_records_by_sample'][row['sample']]
  sx=float(calibration['pixel_x_um']);sy=float(calibration['pixel_y_um'])
  source_sx=float(np.median([float(r['pixel_x_um']) for r in native]))
  if not np.isfinite([sx,sy]).all() or sx<=0 or sy<=0 or sx!=source_sx or sy!=source_sx:
   raise ValueError('Invalid frozen Keyence physical calibration: '+row['sample'])
  if calibration['source_sample']!=row['sample'] or calibration['source_record_count']!=len(native):
   raise ValueError('Frozen Keyence calibration source mapping changed: '+row['sample'])
  row['pixel_x_um']=sx;row['pixel_y_um']=sy;row['_fixed_input_manifest']=str(manifest_path)
  rows.append(row)
 return pd.DataFrame(rows)

def keyence_section(task):
 row,out=task;out=Path(out);old=load_legacy()
 paths={k:Path(row[k+'_labels']) for k in ('DAPI','NPY','cFOS')}
 read_paths={k:Path(row.get('_fixed_mask_paths',paths)[k]) for k in paths}
 registration_path=Path(row.get('_fixed_registration_path',row['current_registration']))
 masks={k:tifffile.imread(p) for k,p in read_paths.items()}
 dapi,npy,cfos=[masks[k] for k in ('DAPI','NPY','cFOS')]
 assert dapi.ndim==2 and dapi.shape==npy.shape==cfos.shape
 fields=pd.DataFrame(regionprops_table(dapi,properties=['label','area','centroid'])).rename(columns={'label':'nucleus_label','area':'area_px','centroid-0':'y_px','centroid-1':'x_px'})
 positives=old.marker_positive_nuclei(dapi,npy,min_overlap_px=old.NPY_MIN_OVERLAP_PX,min_fraction=old.NPY_DAPI_OVERLAP_FRACTION)
 fos,accepted_fos=old.cfos_positive_nuclei(dapi,cfos)
 sx=float(row['pixel_x_um']);sy=float(row['pixel_y_um'])
 reg=json.loads(registration_path.read_text());transform=reg['transform'];seam=transform.get('seam_global_x_px',transform.get('left_retained_width_px'))
 if seam is None:raise ValueError('No traceable source-field seam: '+row['sample']+row['section'])
 side=np.ones(dapi.shape,np.uint8);side[:,int(seam):]=2
 # Match existing Fig3 DAPI support; no new anatomical ROI is invented.
 tissue=ndimage.distance_transform_edt(dapi==0,sampling=(sy,sx))<=15
 roi=np.ones(dapi.shape,np.uint8)
 animal='Keyence_'+row['sample'];aid=animal+'__'+row['section']
 geometry=out/'geometry'/(aid+'.npz');np.savez_compressed(geometry,tissue=tissue,side=side,roi=roi)
 fields['cohort']='NPY';fields['marker']='NPY';fields['animal_id']=animal;fields['condition']=row['condition'];fields['cage']='Keyence_'+row['cage'];fields['genotype']='NPY-GFP';fields['sex']=row['sex'];fields['acquisition_id']=aid;fields['section']=row['section'];fields['region']='FIELD';fields['region_code']=1
 fields['cell_id']='NPY|'+aid+'|N'+fields.nucleus_label.astype(str)
 fields['pixel_x_um']=sx;fields['pixel_y_um']=sy;fields['x_um']=fields.x_px*sx;fields['y_um']=fields.y_px*sy
 ix=np.rint(fields.x_px).astype(int).clip(0,dapi.shape[1]-1);iy=np.rint(fields.y_px).astype(int).clip(0,dapi.shape[0]-1)
 fields['hemifield']=np.where(side[iy,ix]==1,'LEFT_ACQUISITION','RIGHT_ACQUISITION');fields['in_tissue']=tissue[iy,ix];fields['is_cfos']=fields.nucleus_label.isin(fos);fields['accepted_legacy']=fields.nucleus_label.isin(positives);fields['is_npy_legacy']=fields.accepted_legacy
 fields['marker_channel_present']=True;fields['include_acquisition']=True;fields['npy_channel_informative']=True;fields['platform']='Keyence';fields['biological_cohort']='two_bottle_2025'
 for marker,k in [('dapi','DAPI'),('marker','NPY'),('cfos','cFOS')]:fields['source_'+marker+'_mask']=str(paths[k])
 manifest=dict(cohort='NPY',animal_id=animal,condition=row['condition'],cage='Keyence_'+row['cage'],acquisition_id=aid,section=row['section'],genotype='NPY-GFP',sex=row['sex'],pixel_x_um=sx,pixel_y_um=sy,geometry_path=str(geometry),include_acquisition=True,dapi_mask=str(paths['DAPI']),marker_mask=str(paths['NPY']),cfos_mask=str(paths['cFOS']),platform='Keyence',biological_cohort='two_bottle_2025',registration_status=reg.get('status'),region_source='whole sampled field with 15 micrometer DAPI support; separate left/right acquisition geometry',source_record=str(row['current_registration']))
 if '_fixed_input_manifest' in row:
  manifest.update(fixed_dapi_mask=str(read_paths['DAPI']),fixed_marker_mask=str(read_paths['NPY']),fixed_cfos_mask=str(read_paths['cFOS']),fixed_registration_record=str(registration_path),fixed_input_manifest=row['_fixed_input_manifest'])
 check=dict(animal_id=animal,section=row['section'],dapi_nuclei=len(fields),npy_positive_nuclei=len(positives),cfos_positive_nuclei=len(fos),double_positive_nuclei=len(positives&fos),accepted_cfos_rois=accepted_fos,registration_status=reg.get('status'),human_region_review_accepted=row['human_region_review_accepted'],source_mask_hashes={k:shape.sha256(p) for k,p in read_paths.items()})
 if '_fixed_input_manifest' in row:check.update(fixed_input_manifest=row['_fixed_input_manifest'],registration_sha256=shape.sha256(registration_path))
 print('Prepared',aid, 'nuclei',len(fields),'NPY',len(positives),'double',len(positives&fos),flush=True)
 return fields,manifest,check

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=DEFAULT);parser.add_argument('--workers',type=int,default=24);args=parser.parse_args();out=args.output.resolve()
 for folder in ('inputs','geometry','contours','logs'): (out/folder).mkdir(parents=True,exist_ok=True)
 old=load_legacy();legacy=pd.read_csv(SOURCE/'data/NPY_cells.csv.gz',low_memory=False);man=pd.read_csv(SOURCE/'audit/NPY_acquisition_manifest.csv')
 for col in ['is_cfos','accepted_legacy','marker_channel_present','include_acquisition','in_tissue']:legacy[col]=shape.bool_column(legacy,col)
 legacy['platform']=np.where(legacy.condition.eq('Water'),'Leica_SP8','Zeiss_legacy');legacy['biological_cohort']=np.where(legacy.condition.eq('Water'),'legacy_Water','single_bottle_2025')
 legacy.loc[legacy.animal_id.eq('WATER_NPY3'),'biological_cohort']='Water_2026'
 meta=legacy.groupby('animal_id')[['platform','biological_cohort']].first()
 man=man.merge(meta,left_on='animal_id',right_index=True,validate='one_to_one')
 selected=load_fixed_keyence_inputs(out)
 assert len(selected)==15 and selected.groupby('sample').size().eq(5).all()
 assert selected.human_region_review_accepted.all() and selected.canonical_series.all()
 assert set(selected['sample'])=={'FR6-1','FR6-2','FR722'}
 with ProcessPoolExecutor(max_workers=min(args.workers,15)) as pool:new=list(pool.map(keyence_section,[(r,str(out)) for r in selected.to_dict('records')]))
 cells=pd.concat([legacy]+[r[0] for r in new],ignore_index=True,sort=False)
 manifest=pd.concat([man,pd.DataFrame([r[1] for r in new])],ignore_index=True,sort=False)
 assert not cells.cell_id.duplicated().any()
 assert cells.groupby('animal_id').condition.nunique().eq(1).all()
 roster=cells.groupby('animal_id',sort=True).agg(condition=('condition','first'),platform=('platform','first'),biological_cohort=('biological_cohort','first'),cage=('cage','first'),sections=('acquisition_id','nunique')).reset_index()
 assert roster.groupby('condition').size().to_dict()=={'Allulose':5,'Sucrose':4,'Water':3}
 roster.to_csv(out/'animal_roster.csv',index=False);cells.to_csv(out/'inputs/NPY_cells.csv.gz',index=False,compression='gzip');manifest.to_csv(out/'inputs/NPY_acquisition_manifest.csv',index=False)
 shape.write_json(out/'inputs/keyence_association_audit.json',dict(rule_source=str(PAPER/'Fig3/01_analyze_cfos_npy.py'),rule_sha256=shape.sha256(PAPER/'Fig3/01_analyze_cfos_npy.py'),NPY_overlap_pixels=int(old.NPY_MIN_OVERLAP_PX),NPY_overlap_fraction=float(old.NPY_DAPI_OVERLAP_FRACTION),cFOS_min_ROI_area_px=int(old.CFOS_MIN_ROI_AREA_PX),cFOS_overlap_fraction=float(old.CFOS_DAPI_OVERLAP_FRACTION),sections=[r[2] for r in new]))
 counts=[]
 for animal,g in cells.groupby('animal_id',sort=True):
  nd=len(g);nn=int(g.accepted_legacy.sum());nf=int(g.is_cfos.sum());nb=int((g.accepted_legacy&g.is_cfos).sum())
  counts.append(dict(animal=animal,animal_id=animal,condition=g.condition.iloc[0],sample=animal,platform=g.platform.iloc[0],biological_cohort=g.biological_cohort.iloc[0],sections=g.acquisition_id.nunique(),dapi_nuclei=nd,npy_positive_nuclei=nn,cfos_positive_nuclei=nf,double_positive_nuclei=nb,cfos_over_dapi=nf/nd,double_over_npy=nb/nn,npy_over_dapi=nn/nd))
 counts=pd.DataFrame(counts)
 baseline=pd.read_csv(out/'archive/Fig3/source_data/Figure_3_cFos_NPY_animal_plot_values.csv').set_index('animal')
 for r in counts[counts.animal.isin(baseline.index)].itertuples():
  for name in ('cfos_over_dapi','double_over_npy'):assert np.isclose(getattr(r,name),baseline.loc[r.animal,name],rtol=0,atol=1e-12),(r.animal,name)
 counts.to_csv(out/'animal_plot_values.csv',index=False)
 stats=[]
 for endpoint in ['cfos_over_dapi','double_over_npy','npy_over_dapi']:
  groups={condition:counts.loc[counts.condition.eq(condition),endpoint].to_numpy(float) for condition in CONDITIONS};ns=dict(n_water=len(groups['Water']),n_sucrose=len(groups['Sucrose']),n_allulose=len(groups['Allulose']))
  tests=[('one_way_ANOVA','','',old.one_way_anova(list(groups.values())))]
  tests.extend(('exact_two_sided_Mann_Whitney',a,b,old.exact_mwu(groups[a],groups[b])) for a,b in combinations(CONDITIONS,2))
  for test,a,b,result in tests:stats.append(dict(endpoint=endpoint,test=test,group_a=a,group_b=b,statistic=result['statistic'],p_value=result['p_value_raw'],p_adjustment='none; nominal inherited abundance comparison',**ns,enumerated_labelings=result['enumerated_labelings'],extreme_labelings=result['extreme_labelings']))
 pd.DataFrame(stats).to_csv(out/'abundance_statistics.csv',index=False)
 config=dict(k=6,edge_cap_um=75.0,bandwidth_um=40.0,bandwidth_sensitivity=True,grid_um=5.0,contour_mass=.8,minimum_marker_anchors=5)
 tasks=[dict(cells=cells[cells.acquisition_id.eq(r['acquisition_id'])].copy(),manifest=r,config=config,output=str(out)) for r in manifest.to_dict('records')]
 with ProcessPoolExecutor(max_workers=min(args.workers,len(tasks))) as pool:processed=list(pool.map(shape.process_acquisition,tasks))
 strata=pd.DataFrame([row for result in processed for row in result['strata']]);members=pd.DataFrame([row for result in processed for row in result['memberships']],columns=shape.MEMBERSHIP_COLUMNS);contours=pd.DataFrame([row for result in processed for row in result['contours']])
 animals=shape.aggregate_animals(strata);units=shape.unit_metrics(animals)
 # Copy cohort metadata onto the already aggregated animal points.
 animals=animals.merge(roster[['animal_id','platform','biological_cohort']],on='animal_id',validate='many_to_one')
 units=units.merge(roster[['animal_id','platform','biological_cohort']],left_on='unit_id',right_on='animal_id',validate='many_to_one').drop(columns='animal_id')
 summaries=shape.condition_summary(animals,units);cluster_counts=shape.cluster_counts(animals)
 sensitivity=[]
 varied=pd.DataFrame([row for result in processed for row in result['sensitivity']])
 for bandwidth,g in varied.groupby('bandwidth_um'):
  aa=shape.aggregate_animals(g);uu=shape.unit_metrics(aa);ss=shape.condition_summary(aa,uu);ss['bandwidth_um']=bandwidth;sensitivity.append(ss)
 for name,frame in [('animal_metrics',animals),('unit_metrics',units),('condition_summary',summaries),('cluster_counts_by_condition',cluster_counts),('stratum_metrics',strata),('contour_manifest',contours),('bandwidth_sensitivity_condition_summary',pd.concat(sensitivity))]:frame.to_csv(out/(name+'.csv'),index=False)
 members.to_csv(out/'cluster_membership.csv.gz',index=False,compression='gzip')
 # Essential numerical and identity checks, independent of p-value results.
 old_metrics=pd.read_csv(PAPER/'analyses/spatial_shape_20260920/animal_metrics.csv');old_metrics=old_metrics[old_metrics.cohort.eq('NPY')].set_index(['animal_id','endpoint'])
 for r in animals[animals.animal_id.isin(old_metrics.index.get_level_values(0))].itertuples():
  for name in ['N','m','core_enrichment','cluster_fraction','cluster_expected_fraction','cluster_excess']:
   assert np.isclose(getattr(r,name),old_metrics.loc[(r.animal_id,r.endpoint),name],equal_nan=True,atol=1e-10),('legacy changed',r.animal_id,r.endpoint,name)
 assert len(animals)==24 and len(roster)==12
 assert strata.graph_boundary_pass.all() and strata.graph_unique_pass.all()
 assert (strata.clustered_cells<=strata.m).all() and (strata.cluster_expected_cells<=strata.m+1e-8).all()
 assert animals.loc[animals.m.eq(0),'core_enrichment'].isna().all() and animals.loc[animals.m.eq(0),'cluster_fraction'].isna().all()
 assert (members.cluster_size>=2).all() and not members.duplicated(['endpoint','cell_id']).any()
 found=members.groupby(['endpoint','animal_id']).size().to_dict()
 for r in animals.itertuples():assert int(r.clustered_cells)==found.get((r.endpoint,r.animal_id),0)
 assert contours.loc[contours.available,'achieved_mass'].ge(.8-1e-12).all()
 shape.write_json(out/'analysis_validation.json',dict(status='PASS',original_nine_abundance_and_spatial_values_unchanged=True,animals=12,conditions={'Water':3,'Sucrose':4,'Allulose':5},keyence_animals=3,keyence_canonical_sections=15,total_acquisitions=len(manifest),accepted_association_rules='existing Figure3 rules applied to Keyence masks',graphs_cross_source_hemifields=False,missing_conditional_features_not_imputed=True,original_images_and_masks_modified=False))
 shape.write_json(out/'analysis_metadata.json',dict(parameters=config,analysis_script_sha256=shape.sha256(Path(__file__)),shared_shape_script_sha256=shape.sha256(Path(shape.__file__)),inputs=[dict(path=str(p),sha256=shape.sha256(p)) for p in [SOURCE/'data/NPY_cells.csv.gz',SOURCE/'audit/NPY_acquisition_manifest.csv',out/'source_records/keyence_accepted_canonical_sections.csv',out/'source_records/fixed_keyence_inputs.json']],keyence_input_policy='Frozen 45 accepted masks, 15 stitch-registration records and native physical calibration; SHA256 verified against fixed manifest and pre-extension mask audit; original source paths retained only as provenance',interpretation='Exploratory across-cohort comparison; conditional within-platform Allulose-Sucrose inference reported separately',Keyence_animal_weighting='sum count numerators and denominators across five canonical sections; one animal point; no cross-section/hemifield graph edges',source_registration_status_counts=selected.registration_status.value_counts().to_dict(),biological_identity_basis='source_records/identity_evidence.json'))
 print('Analysis complete; 3 Water, 4 Sucrose, 5 Allulose; original nine values preserved.',flush=True)

if __name__=='__main__':main()
