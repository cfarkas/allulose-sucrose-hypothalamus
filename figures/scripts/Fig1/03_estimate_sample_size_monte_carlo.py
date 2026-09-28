#!/usr/bin/env python3
"""Reconstruct prospective sample-size scenarios from single-bottle cage endpoints.

This dated reconstruction is not the historical planning simulation. All cages
are retained with original treatment at the day-6 measurement. Mice sharing a
cage and repeated observations are never counted as independent consumption units.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import f,ncf

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rejection_rates(rng,means,sds,n,reps,test):
    rejected=0
    for start in range(0,reps,2000):
        count=min(2000,reps-start)
        samples=rng.normal(means[None,:,None],sds[None,:,None],size=(count,3,n))
        m=samples.mean(axis=2);v=samples.var(axis=2,ddof=1)
        if test=='ordinary_ANOVA':
            stat=(n*((m-m.mean(axis=1,keepdims=True))**2).sum(axis=1)/2)/v.mean(axis=1)
            p=f.sf(stat,2,3*n-3)
        else:
            weights=n/v;total=weights.sum(axis=1,keepdims=True);weighted=(weights*m).sum(axis=1)/total[:,0]
            correction=(((1-weights/total)**2)/(n-1)).sum(axis=1)/8
            stat=(weights*(m-weighted[:,None])**2).sum(axis=1)/2/(1+2*correction)
            p=f.sf(stat,2,1/(3*correction))
        rejected+=int((p<.05).sum())
    power=rejected/reps
    return power,np.sqrt(power*(1-power)/reps),rejected

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path.cwd());ap.add_argument('--replicates',type=int,default=100000);args=ap.parse_args();root=args.root.resolve();out=root/'analyses/Fig1/results/sample_size_monte_carlo_20260908';out.mkdir(parents=True,exist_ok=True)
    pilot_dir=root/'Fig1/source_data/sample_size_pilot'
    if not pilot_dir.is_dir():pilot_dir=root/'analyses/Fig1/results'
    paths=[pilot_dir/'consumption_primary_values.csv',pilot_dir/'body_weight_cage_endpoint_values.csv']
    a=pd.read_csv(paths[0]);a=a.loc[a.study_day.eq(6)].copy();b=pd.read_csv(paths[1])
    inputs=[('bottle_volume_removed_day6',a,'cumulative_cage_volume_removed_ml','mL/cage'),('body_weight_change_day6',b,'cage_mean_body_weight_change_pct','percent cage-mean change')]
    rng=np.random.default_rng(20260908);curves=[];summaries=[];pilot=[];thresholds=[]
    for endpoint,df,value,unit in inputs:
        assert df.E.nunique()==12 and len(df)==12
        groups=[df.loc[df.Treatment.eq(g),value].dropna().to_numpy() for g in ['Water','Sucrose','Allulose']]
        means=np.array([g.mean() for g in groups]);sds=np.array([g.std(ddof=1) for g in groups]);pooled=np.sqrt(sum(((g-g.mean())**2).sum() for g in groups)/(sum(map(len,groups))-3));effect=np.std(means)/pooled
        for condition,g in zip(['Water','Sucrose','Allulose'],groups):summaries.append({'endpoint':endpoint,'condition':condition,'n_cages':len(g),'mean':g.mean(),'sd':g.std(ddof=1),'unit':unit})
        for _,r in df.iterrows():pilot.append({'endpoint':endpoint,'cage':r.E,'condition':r.Treatment,'value':r[value],'unit':unit})
        for scenario,test,sim_sds in [('common_pooled_variance','ordinary_ANOVA',np.repeat(pooled,3)),('group_specific_variance','Welch_ANOVA',sds)]:
            ns=list(range(3,41));minimum=None
            for n in ns:
                power,se,hits=rejection_rates(rng,means,sim_sds,n,args.replicates,test)
                null,se_null,hits_null=rejection_rates(rng,np.repeat(means.mean(),3),sim_sds,n,20000,test)
                rec={'endpoint':endpoint,'scenario':scenario,'test':test,'n_cages_per_group':n,'total_cages':3*n,'simulations':args.replicates,'rejections':hits,'power':power,'mc_se':se,'null_rejection_rate':null,'null_simulations':20000,'equal_allocation_effect_f':effect,'noncentral_F_power':float(ncf.sf(f.isf(.05,2,3*n-3),2,3*n-3,3*n*effect**2)) if test=='ordinary_ANOVA' else None}
                curves.append(rec)
                if power>=.80 and minimum is None:minimum=rec.copy()
                if minimum and n>=minimum['n_cages_per_group']+2:break
            thresholds.append(minimum or {'endpoint':endpoint,'scenario':scenario,'minimum':'above 40'})
            print(json.dumps(thresholds[-1]),flush=True)
    pd.DataFrame(pilot).to_csv(out/'pilot_cage_values.csv',index=False);pd.DataFrame(summaries).to_csv(out/'pilot_group_parameters.csv',index=False);pd.DataFrame(curves).to_csv(out/'power_curves.csv',index=False);pd.DataFrame(thresholds).to_csv(out/'power_thresholds.csv',index=False)
    receipt={'date':'2026-09-08','purpose':'Author-requested reconstruction from experiment 1; not the original historical planning simulation','seed':20260908,'rng':'NumPy default_rng PCG64','alpha':.05,'target_power':.8,'allocation':'equal numbers of independent cages in three groups','pilot_counts':{'Water':3,'Sucrose':5,'Allulose':4},'distribution':'normal with pilot-estimated means and either pooled or group-specific SD','test_scope':'omnibus group difference, not pairwise contrasts or neuronal/glial outcomes','limitations':'Small pilot groups give uncertain effects and variances; these estimates do not establish adequacy of the completed study, and raw pilot effects may overestimate reproducible differences. Bottle disappearance includes possible spillage. No attrition or multiplicity inflation included.','welfare_transfer':'FR5-4 moved after day-6 measurement; original treatment at measurement retained','source_sha256':{str(p.relative_to(root)):digest(p) for p in paths},'script_sha256':digest(Path(__file__)),'thresholds':thresholds}
    (out/'MONTE_CARLO_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
