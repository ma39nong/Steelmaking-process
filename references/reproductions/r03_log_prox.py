"""Reproduce R03 preprint (22)-(23) with objective-checked nonnegative radial candidates.
This is the proximal core only, not full L3S3TV. Derivation follows the 1/2 scaling
in (22), rather than the inconsistent 1/4 expression printed in the proof of (25).
"""
import numpy as np

def log_radius_prox(a,alpha):
 a=np.asarray(a,float)
 if np.any(a<0) or alpha<0:raise ValueError('Nonnegative radius and alpha required')
 disc=(a+1)**2-4*alpha
 candidate=np.maximum((a-1+np.sqrt(np.maximum(disc,0)))/2,0)
 objective=.5*(candidate-a)**2+alpha*np.log1p(candidate)
 return np.where((disc>=0)&(objective<=.5*a*a),candidate,0)

def column_log_prox(y,alpha):
 norm=np.linalg.norm(y,axis=0);r=log_radius_prox(norm,alpha)
 return y*(r/np.maximum(norm,1e-300))[None,:]

def singular_log_prox(y,alpha):
 u,s,v=np.linalg.svd(y,full_matrices=False)
 return (u*log_radius_prox(s,alpha))@v
