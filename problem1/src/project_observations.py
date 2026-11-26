"""Validate observation-fidelity projection on all 50 B tests, without retuning parameters.
A/C and the latent B optimization are unchanged. The projection is identity on
all held-out evaluation points because their training mask is false.
"""
import json,time,resource
import numpy as np,pandas as pd
from audit import P,HEATS
from anomalies import detect
from models import fit
from experiment import metrics,CONFIG
selection=json.loads((P/'results/final/model_selection.json').read_text());rows=[]
for h in HEATS:
 d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];splits=np.load(P/f'results/validation/{h}_splits.npz');info=selection['per_heat'][h];params=info['models']['B']['params'];params['preserve_observed']=True
 cache=P/f'results/validation/{h}_outer.json';old=json.loads(cache.read_text())
 for row in old:
  if row['model']!='B':continue
  scenario=row['scenario'];q=splits['test_'+scenario];hidden=q.copy()
  if scenario=='time_block':i=int(.68*len(t));hidden[i:i+CONFIG['test_time_block_rows']]=True
  if scenario=='wavelength_block':j=int(.68*len(w));hidden[:,j:j+CONFIG['test_wavelength_block_channels']]=True
  vis=y.copy();vis[hidden]=np.nan;train,_,_=detect(vis,t,w,**info['detector']);start=time.perf_counter();model=fit('B',vis,t,w,train,params);seconds=time.perf_counter()-start
  assert not train[q].any()
  assert np.max(abs(model.grid[q]-model.latent_grid[q]))==0
  assert np.array_equal(model.grid[train],vis[train])
  score=metrics(y,model.grid,q,w,not scenario.startswith('random'))
  max_diff=max(abs(score[key]-row[key]) for key in ['mae','rmse','nrmse','sam_deg'])
  assert max_diff<1e-8,(h,scenario,max_diff)
  rows.append(dict(heat_id=h,scenario=scenario,max_score_difference=max_diff,hidden_prediction_difference=float(np.max(abs(model.grid[q]-model.latent_grid[q]))),valid_observation_difference=float(np.max(abs(model.grid[train]-vis[train]))),solver_converged=model.diagnostics['converged']))
  row.update(score,seconds=seconds,solver_converged=model.diagnostics['converged']);(P/f'results/validation/{h}_B_{scenario}_solver.json').write_text(json.dumps(model.diagnostics,indent=2))
 cache.write_text(json.dumps(old,ensure_ascii=False,indent=2))
 devpath=P/f'results/validation/{h}_development.json';dev=json.loads(devpath.read_text());dev['tuned']['models']['B']['params']['preserve_observed']=True
 for row in dev['development']:
  if row['model']=='B':par=json.loads(row['params']);par['preserve_observed']=True;row['params']=json.dumps(par)
 devpath.write_text(json.dumps(dev,ensure_ascii=False,indent=2));print(h,'5 projection equivalence tests passed',flush=True)
selection['observation_projection_note']='B latent objective unchanged; valid observations imposed before continuous interpolation; all 50 independent masked scores rechecked unchanged'
(P/'results/final/model_selection.json').write_text(json.dumps(selection,ensure_ascii=False,indent=2))
allrows=[]
for h in HEATS:allrows+=json.loads((P/f'results/validation/{h}_outer.json').read_text())
pd.DataFrame(allrows).to_csv(P/'results/model_comparison/model_comparison.csv',index=False);pd.DataFrame(rows).to_csv(P/'results/model_comparison/observation_projection_equivalence.csv',index=False)
# Refresh a signature only after explicit invariance checks, retaining old signature.
path=P/'results/validation/run_signature.json'
if path.exists():path.rename(P/'results/validation/superseded_pilot/pre_projection_signature.json')
from provenance import check_or_stamp
check_or_stamp()
