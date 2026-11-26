import numpy as np,pandas as pd,json
from common import *
from detector import diagnostics
from solver import fit
from models import fit_b

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());selection=json.loads((P/'config/FINAL_SELECTION_LOCKED.json').read_text());v1=json.loads((V1/'results/final/model_selection.json').read_text());cfg=json.loads((P/'config/full_confirmation.json').read_text());rows=[]
 assert selection['model']=='B2'
 for h in HEATS:
  cache=P/f'results/{h}_full_confirmation.json'
  if cache.exists():rows.extend(json.loads(cache.read_text()));continue
  y,t,w=load(h);n,m=CFG['panel_shape'];starts=np.linspace(0,len(t)-n,6).astype(int);testrows=np.zeros(len(t),bool)
  for k in [1,3,5]:testrows[starts[k]:starts[k]+n]=True
  local=[]
  for seed in cfg['seeds']:
   rng=np.random.default_rng(seed);masks=[]
   for f in CFG['random_fractions']:masks.append((f'random_{f:g}',(rng.random(y.shape)<f)&testrows[:,None]))
   for width in CFG['time_widths']:
    hidden=np.zeros(y.shape,bool)
    for k in [1,3,5]:i=int(starts[k]+rng.integers(6,n-width-5));hidden[i:i+width]=True
    masks.append((f'time_{width}',hidden))
   for width in CFG['wave_widths']:
    hidden=np.zeros(y.shape,bool)
    for fraction in [.15,.5,.85]:j=int(fraction*(len(w)-width));hidden[testrows,j:j+width]=True
    masks.append((f'wave_{width}',hidden))
   for scenario,hidden in masks:
    visible=y.copy();visible[hidden]=np.nan;dc=lock['detector'];mask=diagnostics(visible,t,w,**dc)[0] if lock['detector_method']=='auxiliary' else detect(visible,t,w,window=5,threshold=dc['threshold'])[0]
    for name in ['B0',selection['model']]+(['V1_B3_matched'] if seed==cfg['matched_v1_seed'] else []):
     start=__import__('time').perf_counter()
     if name=='V1_B3_matched':r=fit_b(visible,t,w,mask,v1['per_heat'][h]['models']['B']['params']);grid=r.grid;converged=r.diagnostics['converged']
     else:r=fit(visible,t,w,mask,name,**lock['models'][name]['params'],max_iter=CFG['solver_max_iter'],tol=CFG['solver_tol']);grid=r['grid'];converged=r['converged']
     local.append(dict(heat_id=h,model=name,seed=seed,scenario=scenario,seconds=__import__('time').perf_counter()-start,converged=converged,**score(y,grid,hidden,w)))
    print('full confirmation',h,seed,scenario,flush=True)
  dump(cache,local);rows.extend(local);pd.DataFrame(rows).to_csv(P/'results/full_confirmation.csv',index=False)
 dump(P/'results/FULL_CONFIRMATION_COMPLETE.json',dict(rows=len(rows),config_sha256=digest(P/'config/full_confirmation.json')))
if __name__=='__main__':main()
