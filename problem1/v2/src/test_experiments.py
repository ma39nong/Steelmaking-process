"""Locked test evaluation. This file never tunes parameters."""
import json
import numpy as np,pandas as pd
from common import *
from detector import diagnostics,inject,detection_score
from solver import fit

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());assert lock['protocol_sha256']==digest(P/'config/protocol.json')
 for name in ['solver.py','detector.py','common.py']:assert lock['source_sha256'][name]==digest(P/'src'/name),name
 dc=lock['detector'];inj=[];rows=[];hist=[]
 for h,stage,y,t,w in panels('test'):
  cache=P/f'results/test_{h}_{stage}.json'
  if cache.exists():
   a=json.loads(cache.read_text());rows.extend(a['rows']);inj.extend(a['injection']);hist.extend(a['solver']);continue
  local=[];di=[];dh=[]
  for seed in CFG['test_seeds']:
   ix=CFG['test_seeds'].index(seed);strength=[.5,1.5][ix];length=[2,6][ix]
   for kind in ['spike','drop','short_block','row','column','saturation']:
    rng=np.random.default_rng(seed+stage+int(h[-3:]));z,truth=inject(y,kind,rng,strength,length)
    for method in ['MAD','auxiliary']:
     if method=='MAD':base=~detect(y,t,w,window=5,threshold=dc['threshold'])[0];pred=~detect(z,t,w,window=5,threshold=dc['threshold'])[0]
     else:base=~diagnostics(y,t,w,**dc)[0];pred=~diagnostics(z,t,w,**dc)[0]
     di.append(dict(heat_id=h,stage=stage,seed=seed,kind=kind,method=method,strength=strength,length=length,**detection_score(pred,truth,base)))
   for scenario,hidden in hides(y.shape,seed+stage):
    visible=y.copy();visible[hidden]=np.nan
    mask=diagnostics(visible,t,w,**dc)[0] if lock['detector_method']=='auxiliary' else detect(visible,t,w,window=5,threshold=dc['threshold'])[0]
    for variant in ['B0','B1','B2','B3','B4','B5']:
     params=lock['models'][variant]['params'];r=fit(visible,t,w,mask,variant,**params,max_iter=CFG['solver_max_iter'],tol=CFG['solver_tol'])
     s=score(y,r['grid'],hidden,w);s['valid_observation_rmse_not_prediction']=float(np.sqrt(np.mean((r['grid'][mask]-y[mask])**2)))
     local.append(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model=variant,seconds=r['seconds'],converged=r['converged'],iterations=len(r['history']),**s))
     if stage==1 and seed==CFG['test_seeds'][0] and scenario=='time_16':dh.extend(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model=variant,**x) for x in r['history'])
  dump(cache,dict(rows=local,injection=di,solver=dh));rows.extend(local);inj.extend(di);hist.extend(dh)
  pd.DataFrame(rows).to_csv(P/'results/test_completion.csv',index=False);print('test',h,stage,flush=True)
 df=pd.DataFrame(rows);df.to_csv(P/'results/test_completion.csv',index=False);pd.DataFrame(inj).to_csv(P/'results/detector_test.csv',index=False)
 pd.DataFrame(hist).to_csv(P/'results/test_solver_curves.csv',index=False)
 heat=df.groupby(['heat_id','model','scenario']).mean(numeric_only=True);heat.to_csv(P/'results/test_by_heat.csv');heat.groupby(['model','scenario']).agg(['mean','std']).to_csv(P/'results/test_heat_macro.csv')
 pd.DataFrame(inj).groupby(['heat_id','method','kind']).mean(numeric_only=True).groupby(['method','kind']).agg(['mean','std']).to_csv(P/'results/detector_test_heat_macro.csv')
 dump(P/'results/TEST_COMPLETE.json',dict(n=len(df),lock_sha256=digest(P/'config/LOCKED.json'),protocol_sha256=digest(P/'config/protocol.json')))
if __name__=='__main__':main()
