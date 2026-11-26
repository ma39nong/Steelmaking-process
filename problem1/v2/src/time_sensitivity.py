"""Development-only follow-up: no change to any locked model or test choice."""
import numpy as np,pandas as pd
from common import *
import solver

def main():
 rows=[];original=solver.time_weights
 for h,stage,y,t,w in panels('development'):
  for scenario,hidden in hides(y.shape,137+stage):
   if scenario not in ['random_0.05','time_8','wave_16']:continue
   visible=y.copy();visible[hidden]=np.nan;mask=detect(visible,t,w,window=5,threshold=6)[0]
   for weighting in ['uniform','adaptive']:
    solver.time_weights=original if weighting=='adaptive' else lambda z,m:np.ones(len(z)-2)
    for beta in [0,.0001,.001,.01,.03,.3]:
     r=solver.fit(visible,t,w,mask,'B5',rank=4,eta=.003,gamma=.5,beta=beta,max_iter=100,tol=1e-4);rows.append(dict(heat_id=h,stage=stage,scenario=scenario,weighting=weighting,beta=beta,seconds=r['seconds'],converged=r['converged'],**score(y,r['grid'],hidden,w)))
   solver.time_weights=original
  pd.DataFrame(rows).to_csv(P/'results/time_penalty_development_sensitivity.csv',index=False);print('time sensitivity',h,stage,flush=True)
if __name__=='__main__':main()
