import pandas as pd,numpy as np,json
from common import *
from detector import diagnostics
from models import fit_c

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());c=json.loads((P/'config/C_EXPANDED_LOCKED.json').read_text());rows=[]
 assert (P/'config/FINAL_SELECTION_LOCKED.json').exists()
 for h,stage,y,t,w in panels('test'):
  for seed in CFG['test_seeds']:
   for scenario,hidden in hides(y.shape,seed+stage):
    visible=y.copy();visible[hidden]=np.nan;dc=lock['detector'];mask=diagnostics(visible,t,w,**dc)[0] if lock['detector_method']=='auxiliary' else detect(visible,t,w,window=5,threshold=dc['threshold'])[0]
    params=dict(irls=30,cg_iter=250,cg_tol=1e-6,irls_tol=1e-4,objective_tol=1e-6)|c['params'];start=__import__('time').perf_counter();m=fit_c(visible,t,w,mask,params);s=score(y,m.grid,hidden,w);s['valid_observation_rmse_not_prediction']=float(np.sqrt(np.mean((m.grid[mask]-y[mask])**2)));rows.append(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model='C',seconds=__import__('time').perf_counter()-start,converged=m.diagnostics['converged'],**s))
  pd.DataFrame(rows).to_csv(P/'results/c_expanded_test.csv',index=False);print('C expanded test',h,stage,flush=True)
if __name__=='__main__':main()
