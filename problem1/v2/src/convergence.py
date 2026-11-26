"""Full native-matrix solver audit, same visibility, rank/eta and stopping thresholds."""
import time,json
import numpy as np,pandas as pd
from sklearn.utils.extmath import randomized_svd
from common import *
from solver import fit as newfit

def main():
 selected=json.loads((V1/'results/final/model_selection.json').read_text());rows=[];summary=[]
 for h in ['6C02057','6C02058','6C02107','6C02108']:
  y,t,w=load(h);rng=np.random.default_rng(811);hidden=rng.random(y.shape)<.01;i=len(t)//2;hidden[i:i+8]=True;visible=y.copy();visible[hidden]=np.nan
  params=selected['per_heat'][h]['models']['B']['params'];cfg=selected['per_heat'][h]['detector'];mask=detect(visible,t,w,**cfg)[0];q=hidden&(y>0)&(y<64000)
  def error(grid):return dict(held_rmse=float(np.sqrt(np.mean((grid[q]-y[q])**2))))
  z=np.where(mask,visible/65535,0);L=fill_training(np.where(mask,visible/65535,np.nan),t,w);start=time.perf_counter();previous=None
  for k in range(300):
   S=np.where(mask,np.sign(z-L)*np.maximum(abs(z-L)-params['eta'],0),0);a=np.where(mask,z-S,L);u,s,v=randomized_svd(a,n_components=params['rank'],n_iter=3,random_state=20261126);new=(u*s)@v;change=float(np.linalg.norm(new-L)/max(np.linalg.norm(L),1e-12));L=new
   residual=np.where(mask,z-L-S,0);obj=.5*np.sum(residual**2)+params['eta']*np.sum(abs(S))
   r=dict(heat_id=h,solver='V1',iteration=k+1,objective=float(obj),relative_change=change,seconds=time.perf_counter()-start,**error(L*65535));rows.append(r)
   if k+1 in [30,120,300]:
    summary.append(r|dict(prediction_change_from_previous_budget=float(np.sqrt(np.mean(((L*65535)[q]-previous)**2))) if previous is not None else np.nan,converged_1e4=change<1e-4,converged_1e5=change<1e-5));previous=(L*65535)[q].copy()
  pd.DataFrame(rows).to_csv(P/'results/convergence_curves.csv',index=False);print('convergence V1',h,flush=True)
  result=newfit(visible,t,w,mask,'B3',rank=params['rank'],eta=params['eta'],max_iter=300,tol=1e-5,callback=error)
  for r in result['history']:rows.append(dict(heat_id=h,solver='V2',**r))
  for budget in [30,120,300]:
   r=result['history'][min(budget,len(result['history']))-1];summary.append(dict(heat_id=h,solver='V2',requested_budget=budget,**r,converged_1e4=r['relative_change']<1e-4 and r['accepted'],converged_1e5=r['relative_change']<1e-5 and r['accepted']))
  pd.DataFrame(rows).to_csv(P/'results/convergence_curves.csv',index=False);pd.DataFrame(summary).to_csv(P/'results/convergence_budgets.csv',index=False);print('convergence V2',h,flush=True)
if __name__=='__main__':main()
