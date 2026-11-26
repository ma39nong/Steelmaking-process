"""Why change is slow: inspect saved V1 full-data latent iterates, not retune models."""
import numpy as np,pandas as pd,json
from sklearn.utils.extmath import randomized_svd
from common import *
from models import load_surface

def main():
 lock=json.loads((V1/'results/final/model_selection.json').read_text());rows=[]
 for h in ['6C02057','6C02058','6C02107','6C02108']:
  y,t,w=load(h);m=np.load(V1/f'data/masks/{h}.npz')['valid'];model=load_surface(V1/f'data/reconstructed/{h}_model.npz');L=model.latent_grid/65535;params=lock['per_heat'][h]['models']['B']['params'];eta=params['eta'];r=np.where(m,L-y/65535,0);g=np.clip(r,-eta,eta);a=L-g;rank=params['rank'];preds=[]
  for seed in [0,1,2]:
   u,s,v=randomized_svd(a,n_components=rank,n_iter=3,random_state=seed);preds.append((u*s)@v)
  jitter=np.mean([np.linalg.norm(x-preds[0])/max(np.linalg.norm(L),1e-12) for x in preds[1:]])
  u,s,v=randomized_svd(a,n_components=min(rank+2,min(L.shape)),n_iter=7,random_state=137)
  rows.append(dict(heat_id=h,rank=rank,eta=eta,huber_tail_fraction=float(np.mean(abs(r[m])>eta)),clipped_to_quadratic_gradient_norm_ratio=float(np.linalg.norm(g)/max(np.linalg.norm(r),1e-12)),one_step_relative_change=float(np.linalg.norm(preds[0]-L)/max(np.linalg.norm(L),1e-12)),svd_seed_relative_jitter=float(jitter),proposal_rank_gap_ratio=float(s[rank-1]/max(s[rank],1e-12)),warning='Saved V1 full-data latent iterate; differs from masked convergence audit visibility'))
 pd.DataFrame(rows).to_csv(P/'results/solver_slow_change_diagnostics.csv',index=False)
if __name__=='__main__':main()
