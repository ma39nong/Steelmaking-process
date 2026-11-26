"""Re-profile only changed family B in fresh processes; retain unchanged A/C benchmarks."""
import subprocess,sys,json,os
import pandas as pd
from audit import P,HEATS
from pathlib import Path
path=P/'results/model_comparison/resources_isolated.csv';old=pd.read_csv(path);rows=[]
for h in HEATS:
 out=subprocess.run([sys.executable,str(Path(__file__).with_name('benchmark.py')),h,'B'],capture_output=True,text=True,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'));rows.append(json.loads(out.stdout.strip().splitlines()[-1]));print(h,flush=True)
new=pd.concat([old[old.model!='B'],pd.DataFrame(rows)],ignore_index=True).sort_values(['heat_id','model']);new.to_csv(path,index=False)
score=P/'results/model_comparison/model_comparison.csv';d=pd.read_csv(score).drop(columns=['isolated_peak_rss_mb','pre_fit_peak_rss_mb','fit_and_grid_seconds'],errors='ignore');d.merge(new.drop(columns=['solver_converged']),on=['heat_id','model'],validate='many_to_one').to_csv(score,index=False)
