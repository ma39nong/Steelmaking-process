"""All eight saturated heats: fixed local nuclear-constraint tiles, full top-code coverage.
Local separable objectives, not a globally solved nuclear-norm matrix recovery.
"""
import numpy as np,pandas as pd,json
from ctx import *
from preprocess import apply

def main():
 policy=json.loads((P/'config/final_preprocessing_policy.json').read_text())['policy'];cp=CFG['conditional_saturation'];nr,nw=cp['core_shape'];hr,hw=cp['halo'];rows=[]
 for h in HEATS:
  y,t,w=load(h);sat=y==65535
  if not sat.any():continue
  pp=apply(y,t,w,policy);estimate=y.copy();covered=np.zeros(y.shape,bool);converged_mask=np.zeros(y.shape,bool)
  for i in range(0,len(t),nr):
   for j in range(0,len(w),nw):
    core=sat[i:i+nr,j:j+nw]
    if not core.any():continue
    a=max(0,i-hr);b=min(len(t),i+nr+hr);c=max(0,j-hw);d=min(len(w),j+nw+hw);yy=y[a:b,c:d];ss=sat[a:b,c:d];mm=pp['valid'][a:b,c:d]
    r=convex.fit_convex(yy,t[a:b],w[c:d],mm,ss,65535.,tau=cp['tau'],rho=cp['rho'],max_iter=cp['max_iter'],tol=cp['tol']);pred=r['latent'];sub=estimate[i:i+nr,j:j+nw];crop=pred[i-a:i-a+core.shape[0],j-c:j-c+core.shape[1]];sub[core]=crop[core];covered[i:i+nr,j:j+nw]|=core
    converged_mask[i:i+nr,j:j+nw]|=core&r['converged'];last=r['history'][-1]
    rows.append(dict(heat_id=h,core_i=i,core_j=j,saturated_count=int(core.sum()),iterations=len(r['history']),converged=r['converged'],primal_residual=last['primal_residual'],dual_residual=last['dual_residual'],seconds=r['seconds'],lower_violations=int((crop[core]<65535).sum())))
  assert np.array_equal(covered,sat);assert (estimate[sat]>=65535).all()
  np.savez_compressed(P/f'data/{h}_conditional_saturation.npz',estimate=estimate,saturated=sat,covered=covered,solver_converged_at_sat=converged_mask,true_upper_bound_known=False)
  pd.DataFrame(rows).to_csv(P/'results/conditional_saturation_tiles.csv',index=False);print('conditional',h,int(sat.sum()),flush=True)
if __name__=='__main__':main()
