"""Known-reference artificial censoring; cap held fixed per quantile experimental condition."""
import numpy as np,pandas as pd,json
from common import *
from detector import diagnostics
from solver import fit
from convex_censor import fit_convex

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());rows=[];hist=[]
 # gamma sensitivity is predeclared, never used to alter locked choice based on test scores.
 for h,stage,y,t,w in panels('test'):
  for seed in CFG['test_seeds']:
   for quantile in CFG['clip_quantiles']:
    cap=float(np.quantile(y[(y>0)&(y<64000)],quantile));rng=np.random.default_rng(seed+stage)
    for nr,nw in CFG['clip_sizes']:
     i=int(rng.integers(4,len(t)-nr-3));j=int(rng.integers(4,len(w)-nw-3));region=np.zeros(y.shape,bool);region[i:i+nr,j:j+nw]=True
     sat=region&(y>=cap)&(y<64000);visible=y.copy();visible[sat]=cap
     # Hide saturated codes before detecting any other anomalies; explicit sat mask is sent to solver.
     dv=visible.copy();dv[sat]=np.nan;cfg=lock['detector'];mask=diagnostics(dv,t,w,**cfg)[0] if lock['detector_method']=='auxiliary' else detect(dv,t,w,window=5,threshold=cfg['threshold'])[0]
     params=lock['models']['B4']['params'].copy()
     for method,gamma,hard in [('mask',0,False),('soft_0.5',.5,False),('soft_4',4,False),('hard_ADMM',0,True),('convex_nuclear_bound',0,False)]:
      cp=params|dict(gamma=gamma,max_iter=200,tol=1e-4,cap=cap,hard=hard)
      r=fit_convex(visible,t,w,mask,sat,cap) if method=='convex_nuclear_bound' else fit(visible,t,w,mask,'B4',**cp);pred=r['grid'];n=int(sat.sum())
      row=dict(heat_id=h,stage=stage,seed=seed,quantile=quantile,cap=cap,time_size=nr,wave_size=nw,method=method,n=n,rmse=float(np.sqrt(np.mean((pred[sat]-y[sat])**2))) if n else np.nan,lower_bound_violation=float(np.mean(pred[sat]<cap-1e-3)) if n else np.nan,shortfall=float(np.mean(np.maximum(cap-pred[sat],0))) if n else np.nan,seconds=r['seconds'],converged=r['converged'])
      rows.append(row)
      if h in ['6C02057','6C02058','6C02107'] and stage==1 and seed==CFG['test_seeds'][0]:hist.extend(dict(heat_id=h,quantile=quantile,time_size=nr,method=method,**v) for v in r['history'])
  pd.DataFrame(rows).to_csv(P/'results/artificial_censoring.csv',index=False);print('censor',h,stage,flush=True)
 pd.DataFrame(hist).to_csv(P/'results/censor_solver_curves.csv',index=False)
if __name__=='__main__':main()
