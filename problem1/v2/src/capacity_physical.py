"""Development-only correction for changing knot spacing.
Approximate integrated squared-derivative scaling (not exact B-spline Gram integration).
Never used to alter locked models.
"""
import numpy as np,pandas as pd
from common import *
from models import fit_c
rows=[]
for h,stage,y,t,w in panels('development'):
 hidden=dict(hides(y.shape,137+stage))['random_0.05'];visible=y.copy();visible[hidden]=np.nan;mask=detect(visible,t,w,window=5,threshold=9)[0]
 for kt,kw in [(24,48),(48,96),(96,192)]:
  cell_ratio=kt*kw/(24*48);rt=(kt-3)/(24-3);rw=(kw-3)/(48-3)
  params=dict(kt=kt,kw=kw,alpha_t=.001*cell_ratio*rt**3/rw,alpha_w=.001*cell_ratio*rw**3/rt,delta=2.5,irls=30,cg_iter=250,cg_tol=1e-6,irls_tol=1e-4,objective_tol=1e-6)
  m=fit_c(visible,t,w,mask,params);rows.append(dict(heat_id=h,stage=stage,kt=kt,kw=kw,params=json.dumps(params),converged=m.diagnostics['converged'],**score(y,m.grid,hidden,w)))
 pd.DataFrame(rows).to_csv(P/'results/capacity_physical_scaling_development.csv',index=False)
