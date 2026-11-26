"""Convex nuclear-norm alternative with a true lower-bound constraint.
Two-block ADMM; X update is exact entrywise constrained quadratic, Z update SVT.
Feasible X is an optimization variable updated every iteration, never a plotting correction.
"""
import time
import numpy as np
from scipy.linalg import svd
from common import fill_training

def fit_convex(y,t,w,mask,sat,cap,tau=.01,rho=1,max_iter=300,tol=1e-4):
 z=np.where(mask,y/cap,0);X=fill_training(np.where(mask,y/cap,np.nan),t,w);Z=X.copy();U=np.zeros_like(X);start=time.perf_counter();history=[];converged=False
 for k in range(max_iter):
  X=(mask*z+rho*(Z-U))/(mask+rho);X[sat]=np.maximum(X[sat],1)
  old=Z.copy();u,s,v=svd(X+U,full_matrices=False,check_finite=False);shr=np.maximum(s-tau/rho,0);Z=(u*shr)@v;U+=X-Z
  primal=float(np.linalg.norm(X-Z)/max(np.linalg.norm(X),1e-12));dual=float(rho*np.linalg.norm(Z-old)/max(np.linalg.norm(X),1e-12));obj=.5*np.sum((np.where(mask,X-z,0))**2)+tau*np.linalg.svd(X,compute_uv=False).sum()
  history.append(dict(iteration=k+1,objective=float(obj),primal_residual=primal,dual_residual=dual,seconds=time.perf_counter()-start))
  if max(primal,dual)<tol:converged=True;break
 return dict(grid=np.where(mask,y,X*cap),latent=X*cap,history=history,seconds=time.perf_counter()-start,converged=converged)
