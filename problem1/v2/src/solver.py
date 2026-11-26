"""Rank constrained robust completion: monotone projected gradient on eliminated sparse loss.
No hard lower-bound clipping is applied to the plotted reconstruction.
"""
import time
import numpy as np
from scipy.linalg import svd
from scipy.sparse import diags
from sklearn.utils.extmath import randomized_svd
from common import fill_training
SCALE=65535.
def project(a,rank,exact=True):
 if exact:
  u,s,v=svd(a,full_matrices=False,check_finite=False,lapack_driver='gesdd');return (u[:,:rank]*s[:rank])@v[:rank]
 u,s,v=randomized_svd(a,n_components=rank,n_iter=5,random_state=137);return (u*s)@v

def time_weights(z,mask):
 # Only visible pairwise changes; missing changes receive the median visible rate.
 both=mask[1:]&mask[:-1];diff=np.where(both,abs(z[1:]-z[:-1]),np.nan)
 rate=np.nanmedian(diff,axis=1);s=np.nanmedian(rate);rate=np.nan_to_num(rate,nan=s if np.isfinite(s) else 0)
 rate=(rate[:-1]+rate[1:])*.5;scale=max(float(s) if np.isfinite(s) else 0,.001)
 return 1/(1+(rate/scale)**2)

def loss_grad(L,z,mask,eta,sat,gamma,beta,weights):
 r=np.where(mask,L-z,0)
 if eta:
  clipped=np.clip(r,-eta,eta);val=np.sum(np.where(abs(r)<=eta,.5*r*r,eta*abs(r)-.5*eta*eta));grad=clipped
 else:val=.5*np.sum(r*r);grad=r
 if gamma:
  short=np.where(sat,np.minimum(L-1,0),0);val+=.5*gamma*np.sum(short*short);grad=grad+gamma*short
 if beta:
  d=L[:-2]-2*L[1:-1]+L[2:];v=weights[:,None]*d;val+=beta*np.sum(d*v);g=np.zeros_like(L);g[:-2]+=v;g[1:-1]-=2*v;g[2:]+=v;grad=grad+2*beta*g
 return float(val),grad

def fit(y,t,w,mask,variant='B3',rank=8,eta=.003,gamma=.5,beta=.03,max_iter=100,tol=1e-4,exact=None,callback=None,cap=65535.,hard=False):
 start=time.perf_counter();z=np.where(mask,y/SCALE,0);init=fill_training(np.where(mask,y/SCALE,np.nan),t,w)
 if variant=='B0':return dict(grid=init*SCALE,latent=init*SCALE,history=[],converged=True,seconds=time.perf_counter()-start)
 eta=0 if variant=='B1' else eta;gamma=gamma if variant in ['B4','B5'] else 0;beta=beta if variant=='B5' else 0
 sat=np.isfinite(y)&(y==cap)&~mask
 # cap normalization transforms lower bound to one, retaining original recorded-unit output.
 scale=cap
 z=np.where(mask,y/scale,0);init=init*SCALE/scale
 weights=time_weights(z,mask);exact=(y.size<100000) if exact is None else exact
 L=project(init,rank,exact);obj,grad=loss_grad(L,z,mask,eta,sat,gamma,beta,weights)
 hist=[];step=1/(1+gamma+32*beta);base_step=step;converged=False
 # hard alternative uses augmented Lagrangian; returned L is rank constrained, feasibility measured.
 dual=np.zeros_like(L);X=L.copy();rho=2.;outer_res=np.inf
 for k in range(max_iter):
  if hard:
   X=L+dual;X[sat]=np.maximum(X[sat],1);aug=rho*(L-X+dual)
  else:aug=0
  accepted=False
  for back in range(12):
   new=project(L-step*(grad+aug),rank,exact);newobj,newgrad=loss_grad(new,z,mask,eta,sat,gamma,beta,weights)
   oldmerit=obj+(.5*rho*np.sum((L-X+dual)**2) if hard else 0)
   merit=newobj+(.5*rho*np.sum((new-X+dual)**2) if hard else 0)
   if merit<=oldmerit+1e-12:accepted=True;break
   step*=.5
  change=float(np.linalg.norm(new-L)/max(np.linalg.norm(L),1e-12)) if accepted else 0.
  if accepted:L=new;obj=newobj;grad=newgrad
  if hard:
   dual+=L-X;outer_res=float(np.linalg.norm(L-X)/max(np.linalg.norm(L),1e-12))
  row=dict(iteration=k+1,objective=obj,relative_change=change,step=step,accepted=accepted,seconds=time.perf_counter()-start,lower_bound_violation=float(np.mean(L[sat]<1-1e-6)) if sat.any() else 0.,primal_residual=outer_res if hard else 0.)
  if callback:row.update(callback(L*scale))
  hist.append(row)
  if accepted and change<tol and (not hard or outer_res<tol):converged=True;break
  if not accepted:break
  step=min(base_step,step*1.1)
 latent=L*scale;grid=np.where(mask,y,latent) if variant in ['B3','B4','B5'] else latent
 return dict(grid=grid,latent=latent,history=hist,converged=converged,seconds=time.perf_counter()-start,stalled=not hist[-1]['accepted'],solver='exact SVD monotone projected gradient' if exact else 'randomized SVD with objective backtracking',weights=weights)
