import numpy as np,pandas as pd,json
from common import *
from models import fit_c
from detector import diagnostics

def main(split):
 lock=json.loads((P/'config/LOCKED.json').read_text());rows=[]
 if split=='test':cl=json.loads((P/'config/C_LOCKED.json').read_text())
 candidates=[dict(kt=96,kw=192,alpha_t=0,alpha_w=1e-6,delta=1e8),dict(kt=96,kw=192,alpha_t=1e-6,alpha_w=1e-6,delta=20)]
 for h,stage,y,t,w in panels(split):
  for seed in CFG[split+'_seeds']:
   for scenario,hidden in hides(y.shape,seed+stage):
    if split=='development' and scenario not in ['random_0.05','time_8','wave_16']:continue
    visible=y.copy();visible[hidden]=np.nan;cfg=lock['detector'];mask=diagnostics(visible,t,w,**cfg)[0] if lock['detector_method']=='auxiliary' else detect(visible,t,w,window=5,threshold=cfg['threshold'])[0]
    for ci,p in (list(enumerate(candidates)) if split=='development' else [(cl['candidate'],cl['params'])]):
     params=dict(irls=20,cg_iter=250,cg_tol=1e-6,irls_tol=1e-4,objective_tol=1e-6)|p;start=__import__('time').perf_counter();m=fit_c(visible,t,w,mask,params)
     rows.append(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model='C',candidate=ci,seconds=__import__('time').perf_counter()-start,converged=m.diagnostics['converged'],**score(y,m.grid,hidden,w)))
  pd.DataFrame(rows).to_csv(P/f'results/c_{split}.csv',index=False);print('C',split,h,stage,flush=True)
 if split=='development':
  df=pd.DataFrame(rows);a=df.groupby(['candidate','heat_id','scenario']).nrmse.mean().groupby('candidate').mean();ci=int(a.idxmin());dump(P/'config/C_LOCKED.json',dict(candidate=ci,params=candidates[ci],development_macro_nrmse=float(a.loc[ci])))
  options={k:v['development_macro_nrmse'] for k,v in lock['models'].items()};options['C']=float(a.loc[ci]);dump(P/'config/FINAL_LOCKED.json',dict(model=min(options,key=options.get),development_macro_nrmse=options,B_lock_sha256=digest(P/'config/LOCKED.json'),C_lock_sha256=digest(P/'config/C_LOCKED.json'),selection_basis='Only identical disjoint development panels/masks/scenarios/seeds; test unopened'))
if __name__=='__main__':main(__import__('sys').argv[1])
