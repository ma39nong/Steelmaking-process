"""Reuse identical V3 local censored objective; refit tiles whose visible mask changed."""
from core import *
from ctx import convex

def main():
 rows=[];cp=json.loads((V3/'config/protocol.json').read_text())['conditional_saturation'];nr,nw=cp['core_shape'];hr,hw=cp['halo']
 for h in HEATS:
  y,t,w=load(h);sat=y==65535
  if not sat.any():continue
  m=np.load(P/f'data/{h}_masks.npz');repair=m['repair'];valid=m['valid'];old=np.load(V3/f'data/{h}_conditional_saturation.npz');est=old['estimate'].copy();conv=old['solver_converged_at_sat'].copy()
  for i in range(0,len(t),nr):
   for j in range(0,len(w),nw):
    core=sat[i:i+nr,j:j+nw]
    if not core.any():continue
    a=max(0,i-hr);b=min(len(t),i+nr+hr);c=max(0,j-hw);d=min(len(w),j+nw+hw);changed=bool(repair[a:b,c:d].any())
    if changed:
     r=convex.fit_convex(y[a:b,c:d],t[a:b],w[c:d],valid[a:b,c:d],sat[a:b,c:d],65535,tau=cp['tau'],rho=cp['rho'],max_iter=cp['max_iter'],tol=cp['tol']);sub=est[i:i+nr,j:j+nw];crop=r['latent'][i-a:i-a+core.shape[0],j-c:j-c+core.shape[1]];sub[core]=crop[core];conv[i:i+nr,j:j+nw][core]=r['converged'];its=len(r['history']);converged=r['converged']
    else:its=0;converged=bool(conv[i:i+nr,j:j+nw][core].all())
    rows.append(dict(heat_id=h,core_i=i,core_j=j,sat_count=int(core.sum()),refitted=changed,reused_identical_V3_objective=not changed,iterations=its,converged=converged,lower_violations=int((est[i:i+nr,j:j+nw][core]<65535).sum())))
  np.savez_compressed(P/f'data/{h}_conditional.npz',estimate=est,solver_converged_at_sat=conv,sat=sat,true_upper_bound_known=False);pd.DataFrame(rows).to_csv(P/'results/conditional_tiles.csv',index=False);print('conditional',h,flush=True)
if __name__=='__main__':main()
