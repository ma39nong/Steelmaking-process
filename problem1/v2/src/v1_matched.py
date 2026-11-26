import pandas as pd,numpy as np,json
from common import *
from detector import diagnostics
from models import fit_b

def main():
 lock=json.loads((P/'config/LOCKED.json').read_text());v1=json.loads((V1/'results/final/model_selection.json').read_text());rows=[]
 for h,stage,y,t,w in panels('test'):
  for seed in CFG['test_seeds']:
   for scenario,hidden in hides(y.shape,seed+stage):
    visible=y.copy();visible[hidden]=np.nan;dc=lock['detector'];mask=diagnostics(visible,t,w,**dc)[0] if lock['detector_method']=='auxiliary' else detect(visible,t,w,window=5,threshold=dc['threshold'])[0]
    start=__import__('time').perf_counter();m=fit_b(visible,t,w,mask,v1['per_heat'][h]['models']['B']['params']);rows.append(dict(heat_id=h,stage=stage,seed=seed,scenario=scenario,model='V1_B3_matched',seconds=__import__('time').perf_counter()-start,converged=m.diagnostics['converged'],**score(y,m.grid,hidden,w)))
  pd.DataFrame(rows).to_csv(P/'results/v1_matched_test.csv',index=False);print('V1 matched',h,stage,flush=True)
if __name__=='__main__':main()
