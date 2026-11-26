"""Fresh-process memory/time benchmark for each heat/family; peak RSS includes interpreter/input."""
import sys,json,time,resource,subprocess,os
import numpy as np,pandas as pd
from audit import P,HEATS
from models import fit

def child(h,f):
 d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];mask=np.load(P/f'data/masks/{h}.npz')['valid'];selection=json.loads((P/'results/final/model_selection.json').read_text());params=selection['per_heat'][h]['models'][f]['params'];before=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
 start=time.perf_counter();model=fit(f,y,t,w,mask,params);pred=model.grid;elapsed=time.perf_counter()-start
 print(json.dumps(dict(heat_id=h,model=f,fit_and_grid_seconds=elapsed,isolated_peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,pre_fit_peak_rss_mb=before,solver_converged=model.diagnostics['converged'])))
if __name__=='__main__':
 if len(sys.argv)==3:child(sys.argv[1],sys.argv[2])
 else:
  rows=[]
  for h in HEATS:
   for f in ['A','B','C']:
    out=subprocess.run([sys.executable,__file__,h,f],capture_output=True,text=True,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'));rows.append(json.loads(out.stdout.strip().splitlines()[-1]));print(h,f,flush=True)
  df=pd.DataFrame(rows);df.to_csv(P/'results/model_comparison/resources_isolated.csv',index=False)
  path=P/'results/model_comparison/model_comparison.csv'
  if path.exists():
   scores=pd.read_csv(path).drop(columns=['isolated_peak_rss_mb','pre_fit_peak_rss_mb','fit_and_grid_seconds'],errors='ignore').merge(df.drop(columns=['solver_converged']),on=['heat_id','model'],validate='many_to_one');scores.to_csv(path,index=False)
