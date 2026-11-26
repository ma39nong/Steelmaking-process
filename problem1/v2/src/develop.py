"""Development only. Writes immutable LOCKED.json before any independent testing."""
import json,time,itertools
import numpy as np,pandas as pd
from common import *
from detector import diagnostics,inject,detection_score
from solver import fit

def main():
 if (P/'config/LOCKED.json').exists():raise RuntimeError('Selection already locked')
 inj=[]
 for h,stage,y,t,w in panels('development'):
  for ci,cfg in enumerate(CFG['detector_candidates']):
   for seed in CFG['development_seeds']:
    for kind in ['spike','drop','short_block','row','column','saturation']:
     rng=np.random.default_rng(seed+stage+int(h[-3:]));strength=[.5,1.5][CFG['development_seeds'].index(seed)];length=[1,4][CFG['development_seeds'].index(seed)]
     z,truth=inject(y,kind,rng,strength,length)
     for method in ['MAD','auxiliary']:
      if method=='MAD':base=~detect(y,t,w,window=5,threshold=cfg['threshold'])[0];pred=~detect(z,t,w,window=5,threshold=cfg['threshold'])[0]
      else:base=~diagnostics(y,t,w,**cfg)[0];pred=~diagnostics(z,t,w,**cfg)[0]
      inj.append(dict(heat_id=h,stage=stage,seed=seed,kind=kind,candidate=ci,method=method,strength=strength,length=length,**detection_score(pred,truth,base)))
  print('detector dev',h,stage,flush=True)
 df=pd.DataFrame(inj);df.to_csv(P/'results/detector_development.csv',index=False)
 agg=df.groupby(['candidate','method'])[['f1','false_positive_rate']].mean();ok=agg[agg.false_positive_rate<=CFG['detector_max_incremental_false_positive']];best=ok.f1.idxmax() if len(ok) else agg.false_positive_rate.idxmin();dc=dict(CFG['detector_candidates'][int(best[0])]);method=best[1]
 rows=[];candidates=[]
 for rank,eta in itertools.product(CFG['rank_candidates'],CFG['eta_candidates']):candidates.append(dict(rank=rank,eta=eta,gamma=.5,beta=.03))
 # Additional saturation and temporal weights are developed independently below; keep finite grid explicit.
 candidates.extend([dict(rank=8,eta=.003,gamma=4,beta=.03),dict(rank=8,eta=.003,gamma=.5,beta=.3)])
 for h,stage,y,t,w in panels('development'):
  for seed in CFG['development_seeds']:
   # Three geometries weighted equally for model selection, with widths also included in test grid.
   chosen=[(s,q) for s,q in hides(y.shape,seed+stage) if s in ['random_0.05','time_8','wave_16']]
   for scenario,hidden in chosen:
    visible=y.copy();visible[hidden]=np.nan
    mask=(diagnostics(visible,t,w,**dc)[0] if method=='auxiliary' else detect(visible,t,w,window=5,threshold=dc['threshold'])[0])
    for variant in ['B0','B1','B2','B3','B4','B5']:
     cs=[(0,candidates[0])] if variant=='B0' else list(enumerate(candidates[:4] if variant in ['B1','B2','B3'] else candidates))
     for ci,params in cs:
      r=fit(visible,t,w,mask,variant,**params,max_iter=CFG['solver_max_iter'],tol=CFG['solver_tol']);rows.append(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model=variant,candidate=ci,seconds=r['seconds'],converged=r['converged'],**score(y,r['grid'],hidden,w)))
  pd.DataFrame(rows).to_csv(P/'results/development_completion.csv',index=False);print('completion dev',h,stage,flush=True)
 df=pd.DataFrame(rows);selection={};table=[]
 for model,g in df.groupby('model'):
  a=g.groupby(['candidate','heat_id','scenario']).nrmse.mean().groupby('candidate').mean();ci=int(a.idxmin());selection[model]=dict(params=candidates[ci],candidate=ci,development_macro_nrmse=float(a.loc[ci]));table.append((float(a.loc[ci]),model))
 chosen=min(table)[1]
 lock=dict(protocol_sha256=digest(P/'config/protocol.json'),detector=dc,detector_method=method,models=selection,selected_model=chosen,selection_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),source_sha256={p.name:digest(p) for p in (P/'src').glob('*.py')},selection_basis='development only; independent test never loaded')
 dump(P/'config/LOCKED.json',lock);print('LOCKED',chosen,selection,flush=True)
if __name__=='__main__':main()
