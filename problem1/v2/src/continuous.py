import json
import numpy as np,pandas as pd
from common import *
from detector import diagnostics
from interpolation import Surface

def main(split):
 lock=json.loads((P/'config/LOCKED.json').read_text());rows=[]
 if split=='test' and not (P/'config/INTERPOLATION_LOCKED.json').exists():raise RuntimeError('Interpolation must be selected on development only')
 for h,stage,y,t,w in panels(split):
  for axis,stride in [('time',k) for k in CFG['time_decimation']]+[('wave',k) for k in CFG['wave_decimation']]:
   ti=np.arange(len(t));wi=np.arange(len(w))
   if axis=='time':ti=np.unique(np.r_[ti[::stride],len(t)-1])
   else:wi=np.unique(np.r_[wi[::stride],len(w)-1])
   # Only retained coordinates ever reach preprocessing and evaluate model construction.
   visible=y[np.ix_(ti,wi)].copy();tt=t[ti];ww=w[wi];cfg=lock['detector']
   mask=diagnostics(visible,tt,ww,**cfg)[0] if lock['detector_method']=='auxiliary' else detect(visible,tt,ww,window=5,threshold=cfg['threshold'])[0]
   clean=fill_training(np.where(mask,visible,np.nan),tt,ww);hidden=np.ones(y.shape,bool);hidden[np.ix_(ti,wi)]=False
   for method in CFG['interpolation_methods']:
    model=Surface(tt,ww,clean,method);pred=model.evaluate(t,w);s=score(y,pred,hidden,w)
    # Local bracketing rectangle bounds use retained reconstructed observations only.
    it=np.clip(np.searchsorted(tt,t,side='right')-1,0,len(tt)-2);iw=np.clip(np.searchsorted(ww,w,side='right')-1,0,len(ww)-2)
    corners=np.stack([clean[np.ix_(it+a,iw+b)] for a in [0,1] for b in [0,1]]);lo=corners.min(0);hi=corners.max(0)
    over=np.maximum(pred-hi,lo-pred);s['overshoot_fraction']=float(np.mean(over[hidden]>1e-6));s['max_overshoot']=float(np.maximum(over[hidden],0).max())
    rows.append(dict(heat_id=h,stage=stage,scenario=f'{axis}_stride_{stride}',method=method,dt_seconds=.5*stride if axis=='time' else .5,**s))
  print('continuous',split,h,stage,flush=True)
 df=pd.DataFrame(rows);df.to_csv(P/f'results/interpolation_{split}.csv',index=False)
 heat=df.groupby(['heat_id','scenario','method']).mean(numeric_only=True);heat.to_csv(P/f'results/interpolation_{split}_by_heat.csv');heat.groupby(['scenario','method']).agg(['mean','std']).to_csv(P/f'results/interpolation_{split}_macro.csv')
 if split=='development':
  a=df.groupby(['method','heat_id','scenario']).nrmse.mean().groupby('method').mean();dump(P/'config/INTERPOLATION_LOCKED.json',dict(method=a.idxmin(),development_macro_nrmse=a.to_dict(),source_sha256=digest(P/'src/interpolation.py'),selection_basis='disjoint development panels only'))
if __name__=='__main__':main(__import__('sys').argv[1])
