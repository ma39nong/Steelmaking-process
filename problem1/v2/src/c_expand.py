"""Additional C development candidates decided from controlled DEVELOPMENT diagnostics.
No independent completion/interpolation test scores are read.
"""
import numpy as np,pandas as pd,json
from common import *
from models import fit_c
from detector import diagnostics

def main():
 if (P/'config/FINAL_SELECTION_LOCKED.json').exists():raise RuntimeError('Final selection already locked')
 lock=json.loads((P/'config/LOCKED.json').read_text());rows=[];candidates=[]
 for kt,kw in [(24,48),(48,96),(96,192)]:
  for delta in [2.5,20,1e8]:
   ratio=kt*kw/(24*48);candidates.append(dict(kt=kt,kw=kw,alpha_t=.001*ratio,alpha_w=.001*ratio,delta=delta))
 for h,stage,y,t,w in panels('development'):
  for seed in CFG['development_seeds']:
   for scenario,hidden in hides(y.shape,seed+stage):
    if scenario not in ['random_0.05','time_8','wave_16']:continue
    visible=y.copy();visible[hidden]=np.nan;cfg=lock['detector'];mask=diagnostics(visible,t,w,**cfg)[0] if lock['detector_method']=='auxiliary' else detect(visible,t,w,window=5,threshold=cfg['threshold'])[0]
    for ci,p in enumerate(candidates):
     params=dict(irls=30,cg_iter=250,cg_tol=1e-6,irls_tol=1e-4,objective_tol=1e-6)|p;start=__import__('time').perf_counter();m=fit_c(visible,t,w,mask,params);rows.append(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model='C',candidate=ci,seconds=__import__('time').perf_counter()-start,converged=m.diagnostics['converged'],**score(y,m.grid,hidden,w)))
  pd.DataFrame(rows).to_csv(P/'results/c_expanded_development.csv',index=False);print('C expand',h,stage,flush=True)
 df=pd.DataFrame(rows);a=df.groupby(['candidate','heat_id','scenario']).nrmse.mean().groupby('candidate').mean();ci=int(a.idxmin());dump(P/'config/C_EXPANDED_LOCKED.json',dict(candidate=ci,params=candidates[ci],development_macro_nrmse=float(a.loc[ci]),all_candidate_scores=a.to_dict(),full_heat_scaling='Maintain same coefficient density per native coordinate count; alpha scaled to preserve effective coefficient-difference penalty per sample. Not validated full-heat capacity extrapolation.'))
 options={k:v['development_macro_nrmse'] for k,v in lock['models'].items()};options['C']=float(a.loc[ci]);dump(P/'config/FINAL_SELECTION_LOCKED.json',dict(model=min(options,key=options.get),development_macro_nrmse=options,B_lock_sha256=digest(P/'config/LOCKED.json'),C_lock_sha256=digest(P/'config/C_EXPANDED_LOCKED.json'),selection_basis='Identical development panels/scenarios/seeds and masks; expanded C candidates driven by development capacity diagnostic only; completion test unopened'))
if __name__=='__main__':main()
