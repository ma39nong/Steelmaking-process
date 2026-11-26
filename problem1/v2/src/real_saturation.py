import numpy as np,pandas as pd,json
from common import *
from detector import diagnostics
from solver import fit
from convex_censor import fit_convex

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());rows=[];local=[]
 for h in HEATS:
  y,t,w=load(h);sat=y==65535
  if not sat.any():continue
  dc=lock['detector'];mask=diagnostics(y,t,w,**dc)[0] if lock['detector_method']=='auxiliary' else detect(y,t,w,window=5,threshold=dc['threshold'])[0];params=lock['models']['B4']['params'].copy()
  for name,gamma,hard in [('mask',0,False),('soft_0.5',.5,False),('soft_4',4,False),('hard_ADMM',0,True)]:
   r=fit(y,t,w,mask,'B4',**(params|dict(gamma=gamma,max_iter=300,tol=1e-4,hard=hard)));pred=r['grid'];rows.append(dict(heat_id=h,method=name,count=int(sat.sum()),lower_violation=float(np.mean(pred[sat]<65535-1e-3)),shortfall=float(np.mean(np.maximum(65535-pred[sat],0))),mean_reconstruction=float(pred[sat].mean()),seconds=r['seconds'],converged=r['converged'],actual_error_identifiable=False));dump(P/f'results/{h}_{name}_real_sat_solver.json',dict(history=r['history'],converged=r['converged']))
   np.savez_compressed(P/f'data/{h}_{name}_real_sat.npz',time=t,wavelength=w,saturated_prediction=pred[sat],saturated=sat)
  if h in ['6C02057','6C02058','6C02107']:
   # Window centered at the highest saturation-density row/column; selection uses no unknown latent height.
   counts=sat.sum(1);i=max(0,min(int(counts.argmax())-48,len(t)-96));j=max(0,min(int(sat[i:i+96].sum(0).argmax())-96,len(w)-192));yy=y[i:i+96,j:j+192];mm=mask[i:i+96,j:j+192];ss=sat[i:i+96,j:j+192]
   r=fit_convex(yy,t[i:i+96],w[j:j+192],mm,ss,65535.)
   local.append(dict(heat_id=h,method='convex_nuclear_bound_local',row_start=i,column_start=j,local_saturation_count=int(ss.sum()),whole_heat_sat_count=int(sat.sum()),local_coverage=float(ss.sum()/sat.sum()),lower_violation=float(np.mean(r['grid'][ss]<65535-1e-3)),seconds=r['seconds'],converged=r['converged'],actual_error_identifiable=False))
   np.savez_compressed(P/f'data/{h}_local_convex_sat.npz',time=t[i:i+96],wavelength=w[j:j+192],grid=r['grid'],saturated_unidentifiable=ss)
  pd.DataFrame(rows).to_csv(P/'results/real_saturation_comparison.csv',index=False);pd.DataFrame(local).to_csv(P/'results/focus_local_convex_saturation.csv',index=False);print('real saturation',h,flush=True)
if __name__=='__main__':main()
