"""Controlled C diagnostics on development only; source v1 fit_c kept unchanged."""
import numpy as np,pandas as pd
from common import *
from models import fit_c

def main():
 rows=[]
 for h,stage,y,t,w in panels('development'):
  # All heats, all stages; same visible split for each one-factor setting.
  hidden=dict(hides(y.shape,137+stage))['random_0.05'];visible=y.copy();visible[hidden]=np.nan
  mask=detect(visible,t,w,window=5,threshold=9)[0]
  base=dict(kt=24,kw=48,alpha_t=.001,alpha_w=.001,delta=2.5,irls=30,cg_iter=250,cg_tol=1e-6,irls_tol=1e-4,objective_tol=1e-6)
  settings=[('base',{}),('capacity_48_96',dict(kt=48,kw=96)),('capacity_96_192',dict(kt=96,kw=192)),('time_penalty_zero',dict(alpha_t=0)),('time_penalty_0.03',dict(alpha_t=.03)),('huber_delta_20',dict(delta=20)),('quadratic',dict(delta=1e8)),('capacity_full_quadratic',dict(kt=96,kw=192,delta=1e8,alpha_t=0,alpha_w=1e-6))]
  for name,change in settings:
   params=base|change
   # Keep effective penalties fixed for capacity-only contrasts.
   if name.startswith('capacity_') and name!='capacity_full_quadratic':
    ratio=(params['kt']*params['kw'])/(base['kt']*base['kw']);params['alpha_t']*=ratio;params['alpha_w']*=ratio
   start=__import__('time').perf_counter();m=fit_c(visible,t,w,mask,params);rows.append(dict(heat_id=h,stage=stage,setting=name,params=json.dumps(params),scale_intensity=m.diagnostics['scale_intensity'],seconds=__import__('time').perf_counter()-start,converged=m.diagnostics['converged'],**score(y,m.grid,hidden,w)))
  pd.DataFrame(rows).to_csv(P/'results/capacity_development.csv',index=False);print('capacity',h,stage,flush=True)
if __name__=='__main__':main()
