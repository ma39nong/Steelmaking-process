"""Finite native and held-panel iteration studies; fixed stopping and stability thresholds."""
import sys,json,time
import numpy as np,pandas as pd
from ctx import *
from preprocess import apply
from model3 import solve,project_record

def main(h):
 assert (P/'results/PREPROCESSING_ACCEPTED.json').exists();policy=json.loads((P/'config/final_preprocessing_policy.json').read_text())['policy'];y,t,w=load(h);pp=apply(y,t,w,policy)
 r=solve(y,t,w,pp['valid'],weights=pp['weights'],max_iter=1200,checkpoints=[300,600,1200],wall_seconds=CFG['per_heat_wall_seconds'])
 np.savez_compressed(P/f'data/{h}_native_fit.npz',latent=r['latent']);dump(P/f'results/{h}_native_solver.json',{k:v for k,v in r.items() if k not in ['latent','checkpoints']})
 records=[];previous=None
 for budget,L in sorted(r['checkpoints'].items()):
  D,_,_=project_record(y,pp['valid'],L,t,w,True,True);sat=pp['sat'];missing=~pp['valid']&~sat
  records.append(dict(heat_id=h,scope='native',budget=budget,relative_change=r['history'][budget-1]['relative_change'],objective=r['history'][budget-1]['objective'],strict_converged=r['history'][budget-1]['relative_change']<1e-4,latent_sat_change_rmse=float(np.sqrt(np.mean((L[sat]-previous[0][sat])**2))) if previous is not None and sat.any() else 0,record_output_change_rmse=float(np.sqrt(np.mean((D-previous[1])**2))) if previous is not None else np.nan,record_peak_change_max=float(np.max(abs(D.max(1)-previous[1].max(1)))) if previous is not None else np.nan))
  previous=(L,D)
 # Held confirmation target is chosen by protocol before any scores. No full-native warm start is allowed here.
 for hh,stage,yy,tt,ww,i,j in panels('model'):
  if hh!=h or stage!=1:continue
  hidden=dict(hides(yy.shape,CFG['confirmation_seeds'][0]))['random_0.05'];visible=yy.copy();visible[hidden]=np.nan;vp=apply(visible,tt,ww,policy)
  s=solve(visible,tt,ww,vp['valid'],weights=vp['weights'],max_iter=1200,checkpoints=[300,600,1200],wall_seconds=120);previous=None
  for budget,L in sorted(s['checkpoints'].items()):
   D,_,_=project_record(visible,vp['valid'],L,tt,ww,True,True);delta=D[hidden]-previous[hidden] if previous is not None else np.array([np.nan]);records.append(dict(heat_id=h,scope='held_panel',budget=budget,relative_change=s['history'][budget-1]['relative_change'],objective=s['history'][budget-1]['objective'],strict_converged=s['history'][budget-1]['relative_change']<1e-4,prediction_change_rms=float(np.sqrt(np.mean(delta**2))),prediction_change_q95=float(np.quantile(abs(delta),.95)),held_rmse=score(yy,D,hidden,ww)['rmse'],record_peak_change_max=float(np.max(abs(D.max(1)-previous.max(1)))) if previous is not None else np.nan));previous=D
  dump(P/f'results/{h}_held_solver.json',{k:v for k,v in s.items() if k not in ['latent','checkpoints']})
 pd.DataFrame(records).to_csv(P/f'results/{h}_stability.csv',index=False);print('stability',h,r['converged'],len(r['history']),r['stop_reason'],flush=True)
if __name__=='__main__':main(sys.argv[1])
