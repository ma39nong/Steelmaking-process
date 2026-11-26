import json,time
import numpy as np,pandas as pd
from common import *
from detector import diagnostics
from solver import fit
from interpolation import Surface

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());final=json.loads((P/'config/FINAL_SELECTION_LOCKED.json').read_text());inter=json.loads((P/'config/INTERPOLATION_LOCKED.json').read_text());rows=[]
 assert final['model']=='B2'
 for h in HEATS:
  cache=P/f'results/{h}_full_interpolation.json'
  if cache.exists():rows.extend(json.loads(cache.read_text()));continue
  y,t,w=load(h);local=[]
  for axis,stride in [('time',2),('time',3),('time',4),('wave',2),('wave',3)]:
   ti=np.arange(len(t));wi=np.arange(len(w))
   if axis=='time':ti=np.unique(np.r_[ti[::stride],len(t)-1])
   else:wi=np.unique(np.r_[wi[::stride],len(w)-1])
   tt=t[ti];ww=w[wi];visible=y[np.ix_(ti,wi)];dc=lock['detector'];mask=diagnostics(visible,tt,ww,**dc)[0] if lock['detector_method']=='auxiliary' else detect(visible,tt,ww,window=5,threshold=dc['threshold'])[0];clean=fill_training(np.where(mask,visible,np.nan),tt,ww)
   hidden=np.ones(y.shape,bool);hidden[np.ix_(ti,wi)]=False
   for method in CFG['interpolation_methods']+['locked_completion_plus_locked_interpolation']:
    start=time.perf_counter()
    if method.startswith('locked_'):r=fit(visible,tt,ww,mask,final['model'],**lock['models'][final['model']]['params'],max_iter=CFG['solver_max_iter'],tol=CFG['solver_tol']);grid=r['grid'];interp=inter['method']
    else:grid=clean;interp=method
    model=Surface(tt,ww,grid,interp);pred=model.evaluate(t,w);s=score(y,pred,hidden,w)
    it=np.clip(np.searchsorted(tt,t,side='right')-1,0,len(tt)-2);iw=np.clip(np.searchsorted(ww,w,side='right')-1,0,len(ww)-2)
    # Accumulate bounds without stacking four full matrices.
    lo=grid[np.ix_(it,iw)].copy();hi=lo.copy()
    for a,b in [(0,1),(1,0),(1,1)]:c=grid[np.ix_(it+a,iw+b)];np.minimum(lo,c,out=lo);np.maximum(hi,c,out=hi)
    over=np.maximum(pred-hi,lo-pred);s['overshoot_fraction']=float(np.mean(over[hidden]>1e-6));s['max_overshoot']=float(np.maximum(over[hidden],0).max());local.append(dict(heat_id=h,scenario=f'{axis}_stride_{stride}',method=method,seconds=time.perf_counter()-start,**s))
   print('full continuous',h,axis,stride,flush=True)
  dump(cache,local);rows.extend(local);pd.DataFrame(rows).to_csv(P/'results/full_interpolation.csv',index=False)
 dump(P/'results/FULL_INTERPOLATION_COMPLETE.json',dict(rows=len(rows),config_sha256=digest(P/'config/full_interpolation.json')))
if __name__=='__main__':main()
