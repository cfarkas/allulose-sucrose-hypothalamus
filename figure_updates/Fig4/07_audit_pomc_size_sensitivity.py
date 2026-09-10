"""Reassign nuclei at three fixed POMC size thresholds; retain original data.

The default 1.0 threshold must reproduce every main-analysis POMC assignment.
The other thresholds characterize sensitivity, not selection by significance.
"""
from pathlib import Path
import importlib.util,sys,json,numpy as np,pandas as pd,tifffile
from scipy.stats import f_oneway,mannwhitneyu
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/shared/pomc_size_qc.py').is_file());sys.path.insert(0,str(ROOT/'scripts/shared'));from pomc_size_qc import filter_pomc_by_local_dapi
spec=importlib.util.spec_from_file_location('pomc_analysis',ROOT/'Fig4/02_analyze_pomc_cfos.py');m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);m.import_science_stack()
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--analysis-dir',type=Path,default=ROOT/'analyses/Fig4/results/human_final_run');args=parser.parse_args();r=args.analysis_dir.resolve();cells=pd.read_csv(r/'per_cell_reconstructed_measurements.csv',low_memory=False);keys=['sample','section_index','tile'];cells['tile']=cells.tile.fillna('');ss=pd.read_csv(r/'per_reconstructed_section_summary.csv');ss['tile']=ss.tile.fillna('');records=[]
for row in ss.itertuples():
 folder=r/'reconstructed_sections'/row.animal_id/(f'{row.sample}_S{row.section_index:02d}'+(f'_T{row.tile}' if row.tile else ''))
 d=tifffile.imread(folder/'reconstructed_DAPI_labels.tif');p=tifffile.imread(folder/'reconstructed_POMC_labels.tif');n=tifffile.imread(folder/'reconstructed_NPY_labels.tif');region=tifffile.imread(folder/'reconstructed_ARC_ME_accepted_mask.tif')
 c=cells[(cells['sample']==row.sample)&(cells.section_index==row.section_index)&(cells.tile==row.tile)].copy()
 npy=m.nuclei_positive_by_expanded_label(d,n,round(7.5/.755),5)
 for factor in [.5,1.,1.5]:
  filtered,audit=filter_pomc_by_local_dapi(p,d,region,area_factor=factor)
  pos=m.nuclei_positive_by_expanded_label(d,filtered,round(7.5/.755),5);pos,_=m.apply_strict_pomc_vs_npy(d,npy,pos,n,filtered,.60)
  predicted=c.nucleus_label.isin(pos)
  if factor==1.:assert np.array_equal(predicted.to_numpy(),c.is_pomc.to_numpy()),row.sample
  for reg in ['ARC','ME']:
   keep=c.region.eq(reg);yes=keep&predicted
   if keep.any():records.append({'factor':factor,'animal_id':row.animal_id,'sample':row.sample,'section':row.section_index,'region':reg,'condition':c.cond.iloc[0],'cage':c.cage.iloc[0],'pomc':int(yes.sum()),'double':int((yes&c.is_cfos).sum())})

a=pd.DataFrame(records);a.to_csv(r/'pomc_size_sensitivity_section_counts.csv',index=False);a=a[~a.condition.str.contains('negative',case=False,na=False)]
an=a.groupby(['factor','animal_id','region','condition','cage'],dropna=False)[['pomc','double']].sum().reset_index();an['fraction']=an['double']/an.pomc.replace(0,np.nan);an.to_csv(r/'pomc_size_sensitivity_animal_values.csv',index=False)
st=[]
for (factor,reg),g in an.groupby(['factor','region']):
 gs={cond:g.loc[g.condition.eq(cond),'fraction'].dropna().to_numpy() for cond in ['Water','Sucrose','Allulose']}
 st.append({'factor':factor,'region':reg,'test':'ANOVA','comparison':'global','p_nominal':f_oneway(*gs.values()).pvalue})
 for one,two in [('Water','Sucrose'),('Water','Allulose'),('Sucrose','Allulose')]:
  st.append({'factor':factor,'region':reg,'test':'exact Mann-Whitney','comparison':one+' vs '+two,'p_nominal':mannwhitneyu(gs[one],gs[two],alternative='two-sided',method='exact').pvalue,'mean_a':np.mean(gs[one]),'mean_b':np.mean(gs[two])})
s=pd.DataFrame(st);s.to_csv(r/'pomc_size_sensitivity_statistics.csv',index=False);print(s[s.region=='ARC'].to_string(index=False))
