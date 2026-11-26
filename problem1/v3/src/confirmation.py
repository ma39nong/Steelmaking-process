"""Finite candidate gate on gap 0, then frozen confirmation on gap 4."""
import sys,json,time
import numpy as np,pandas as pd
from ctx import *
from preprocess import apply
from model3 import project_record

def candidates(y,t,w,pp,h):
 mask=pp['valid'];par=json.loads((V1/'results/final/model_selection.json').read_text())['per_heat'][h]['models']['B']['params'];v1=fit_b(y,t,w,mask,par)
 r=legacy.fit(y,t,w,mask,'B2',rank=4,eta=.003,max_iter=100,tol=1e-4);L=r['latent'];q=legacy.fit(y,t,w,mask,'B1',rank=4,eta=0,max_iter=100,tol=1e-4)
 grids={'V1':v1.grid,'V2':r['grid']};details={}
 for name,latent,guard,sat in [('V3-A',L,False,False),('V3-B',q['latent'],False,False),('V3-C',L,True,True)]:
  grids[name],_,details[name]=project_record(y,mask,latent,t,w,guard,sat)
 assert np.array_equal(grids['V2'][~mask&~pp['sat']],grids['V3-A'][~mask&~pp['sat']])
 costs={name:r['seconds'] for name in ['V2','V3-A','V3-C']};costs['V3-B']=q['seconds'];costs['V1']=np.nan
 return grids,details,costs

def physical(y,grid,pp,t,w):
 valid=pp['valid'];sat=pp['sat'];_,_,d=project_record(y,valid,grid,t,w,False,True);q=~(valid|sat)
 return dict(valid_max_error=float(np.max(abs(grid[valid]-y[valid]))) if valid.any() else 0,negative_count=int((grid<0).sum()),saturated_record_max_error=float(np.max(abs(grid[sat]-65535))) if sat.any() else 0,unsupported_missing_count=int((d['unsupported']&q).sum()))

def main(split):
 assert (P/'results/PREPROCESSING_ACCEPTED.json').exists();policy=json.loads((P/'config/final_preprocessing_policy.json').read_text())['policy'];stagewanted=0 if split=='gate' else 1
 if split!='gate':assert (P/'config/FINAL_MODEL_LOCKED.json').exists()
 rows=[];continuous=[];equiv=[]
 for h,stage,y,t,w,i,j in panels('model'):
  if stage!=stagewanted:continue
  for seed in CFG['confirmation_seeds']:
   for scene,hidden in hides(y.shape,seed):
    visible=y.copy();visible[hidden]=np.nan;pp=apply(visible,t,w,policy);grids,details,costs=candidates(visible,t,w,pp,h)
    equiv.append(dict(heat_id=h,seed=seed,scenario=scene,max_missing_difference=float(np.max(abs(grids['V2'][hidden]-grids['V3-A'][hidden])))))
    for name,grid in grids.items():
     m=Surface(t,w,grid,'pchip' if name=='V1' else 'linear');pred=m.evaluate(t,w)
     rows.append(dict(heat_id=h,seed=seed,scenario=scene,model=name,seconds=costs[name],**physical(visible,grid,pp,t,w),**score(y,pred,hidden,w)))
    base=fill_training(np.where(pp['valid']|pp['sat'],visible,np.nan),t,w);rows.append(dict(heat_id=h,seed=seed,scenario=scene,model='linear_baseline',seconds=0,**score(y,base,hidden,w)))
  for axis in ['time','wave']:
   ti=np.arange(len(t));wi=np.arange(len(w))
   if axis=='time':ti=np.unique(np.r_[ti[::2],len(t)-1])
   else:wi=np.unique(np.r_[wi[::2],len(w)-1])
   tt=t[ti];ww=w[wi];visible=y[np.ix_(ti,wi)];pp=apply(visible,tt,ww,policy);grids,_,costs=candidates(visible,tt,ww,pp,h);hidden=np.ones(y.shape,bool);hidden[np.ix_(ti,wi)]=False
   for name,grid in grids.items():
    m=Surface(tt,ww,grid,'pchip' if name=='V1' else 'linear');pred=m.evaluate(t,w);continuous.append(dict(heat_id=h,scenario=axis+'_stride_2',model=name,seconds=costs[name],**score(y,pred,hidden,w)))
  print('candidate',split,h,flush=True)
 pd.DataFrame(rows).to_csv(P/f'results/model_{split}.csv',index=False);pd.DataFrame(continuous).to_csv(P/f'results/continuous_{split}.csv',index=False);pd.DataFrame(equiv).to_csv(P/f'results/projection_equivalence_{split}.csv',index=False)
 if split=='gate':
  d=pd.DataFrame(rows);scores=d.groupby(['heat_id','model','scenario']).nrmse.mean().groupby('model').mean();accepted=[];gate=[]
  for name in ['V1','V2','V3-A','V3-B','V3-C']:
   g=d[d.model==name];physical_ok=bool((g.valid_max_error<=1e-8).all() and (g.negative_count==0).all() and (g.saturated_record_max_error<=1e-8).all() and (g.unsupported_missing_count==0).all())
   # Recorded-censor preservation is a required method-level property across all eight known saturated heats.
   censor_supported=name=='V3-C';ok=physical_ok and censor_supported;gate.append(dict(model=name,physical_samples_pass=physical_ok,recorded_censor_supported=censor_supported,eligible=ok,macro_nrmse=float(scores[name])))
   if ok:accepted.append(name)
  if not accepted:raise RuntimeError('No physically eligible candidate; do not tune confirmation')
  chosen=min(accepted,key=lambda x:scores[x]);dump(P/'config/FINAL_MODEL_LOCKED.json',dict(model=chosen,rank=4,eta=.003,preprocessing=policy,interpolation='linear',gates=gate,selection_basis='Predeclared physical priorities on new gap0; gap4 confirmation unopened',protocol_sha256=sha(P/'config/protocol.json'),policy_sha256=sha(P/'config/final_preprocessing_policy.json')))
if __name__=='__main__':main(sys.argv[1])
