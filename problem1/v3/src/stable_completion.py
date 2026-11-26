"""Predeclared finite-iteration safeguard for future missing native cells.
Does not modify locked 100-step model-confirmation results. No truth enters fallback.
"""
import numpy as np
from ctx import fill_training,CFG
from model3 import solve,project_record

def complete_with_stability(y,t,w,valid,weights=None):
 r=solve(y,t,w,valid,weights=weights,max_iter=1200,checkpoints=[600,1200],wall_seconds=360)
 end=project_record(y,valid,r['latent'],t,w,True,True)[0]
 missing=~valid&~(y==65535);fallback=np.zeros(y.shape,bool)
 previous=r['checkpoints'].get(600)
 if previous is None:
  fallback=missing.copy() # budget exhaustion before comparison: unverified predictions
 else:
  a=project_record(y,valid,previous,t,w,True,True)[0];fallback=missing&(abs(end-a)>CFG['stability_output_rms_tolerance_counts'])
 if fallback.any():
  base=fill_training(np.where(valid|(y==65535),y,np.nan),t,w)
  end[fallback]=base[fallback]
 return dict(grid=end,numerical_fallback=fallback,latent_result=r)
