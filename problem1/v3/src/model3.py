"""Fixed-rank weighted robust completion plus an explicit recorded-scale box projection."""
import time
import numpy as np
from ctx import legacy,fill_training,Surface,CFG
U=65535.

def solve(y,t,w,mask,eta=.003,rank=4,weights=None,max_iter=300,checkpoints=(),wall_seconds=360):
 start=time.perf_counter();z=np.where(mask,y/U,0);W=mask.astype(float) if weights is None else np.where(mask,weights,0);L=legacy.project(fill_training(np.where(mask,y/U,np.nan),t,w),rank,y.size<100000)
 def fg(a):
  r=np.where(mask,a-z,0)
  if eta:return float(np.sum(W*np.where(abs(r)<=eta,.5*r*r,eta*abs(r)-.5*eta**2))),W*np.clip(r,-eta,eta)
  return float(.5*np.sum(W*r*r)),W*r
 obj,g=fg(L);history=[];saved={};step=1.;converged=False;stop_reason='iteration_budget'
 for k in range(max_iter):
  accepted=False
  for bt in range(12):
   new=legacy.project(L-step*g,rank,y.size<100000);val,grad=fg(new)
   if val<=obj+1e-12:accepted=True;break
   step*=.5
  change=float(np.linalg.norm(new-L)/max(np.linalg.norm(L),1e-12)) if accepted else np.nan
  if accepted:L=new;obj=val;g=grad
  history.append(dict(iteration=k+1,objective=obj,relative_change=change,step=step,accepted=accepted,seconds=time.perf_counter()-start));converged=bool(accepted and change<1e-4)
  if k+1 in checkpoints:saved[k+1]=L.copy()*U
  if not accepted:stop_reason='backtracking_stalled';break
  if converged and not checkpoints:stop_reason='relative_change';break
  if time.perf_counter()-start>wall_seconds:stop_reason='wall_budget';break
  step=min(1,step*1.1)
 saved[len(history)]=L*U
 return dict(latent=L*U,history=history,checkpoints=saved,converged=converged,stop_reason=stop_reason,seconds=time.perf_counter()-start)

def support_bounds(y,observed):
 """Nearest recorded axial supports, only visible values. No hidden references used."""
 n,m=y.shape;lo=np.full(y.shape,np.inf);hi=np.full(y.shape,-np.inf);count=np.zeros(y.shape,np.uint8)
 for axis in [0,1]:
  a=y if axis==0 else y.T;q=observed if axis==0 else observed.T;l=lo if axis==0 else lo.T;h=hi if axis==0 else hi.T;c=count if axis==0 else count.T
  for j in range(a.shape[1]):
   ix=np.flatnonzero(q[:,j])
   if not len(ix):continue
   p=np.searchsorted(ix,np.arange(a.shape[0]));left=ix[np.maximum(p-1,0)];right=ix[np.minimum(p,len(ix)-1)]
   for ref in [left,right]:v=a[ref,j];l[:,j]=np.minimum(l[:,j],v);h[:,j]=np.maximum(h[:,j],v);c[:,j]+=1
 weak=count==0
 if not observed.any():raise ValueError('No visible support')
 lo[weak]=np.min(y[observed]);hi[weak]=np.max(y[observed]);lo=np.maximum(lo,0)
 return lo,hi,weak

def project_record(y,valid,latent,t,w,guard=False,retain_saturation=False):
 sat=np.isfinite(y)&(y==U);observed=valid|(sat if retain_saturation else False);D=np.where(observed,y,latent).copy();quality=np.zeros(y.shape,np.uint16);quality[sat]|=2;quality[~valid&~sat]|=4
 lo,hi,weak=support_bounds(y,observed);unsupported=~observed&((D<lo-1e-8)|(D>hi+1e-8));negative=D<0
 if guard:
  D[~observed]=np.clip(D[~observed],lo[~observed],hi[~observed]);quality[unsupported]|=8;quality[negative]|=16;quality[weak]|=32
 return D,quality,dict(lo=lo,hi=hi,unsupported=unsupported,negative=negative,observed=observed)
