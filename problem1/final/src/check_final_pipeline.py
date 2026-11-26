"""Fixed selected detector, SAME confirmation cases, exact final V3 stability safeguard.
This checks implementation transfer, not a second independent test or parameter search.
"""
from core import *
from detectors import adaptive,local_mask
from experiment_final import metrics

def main():
 chosen=json.loads((P/'config/SELECTION_LOCKED.json').read_text())['strategy'];assert chosen=='D_adaptive4';rows=[]
 for h,b,y,t,w,i,j in panels('confirmation'):
  for seed in CFG['confirmation_seeds']:
   for kind in ['spike','drop']:
    for strength in [.3,.7]:
     corrupt,truth=inject(y,np.random.default_rng(seed+sum(map(ord,h))+j),kind,strength);e=adaptive(corrupt,t,w);mask=local_mask(e,4);valid=(corrupt!=65535)&~mask
     if mask.any():r=complete_with_stability(corrupt,t,w,valid);D=r['grid'];converged=r['latent_result']['converged'];fallback=int(r['numerical_fallback'][mask].sum())
     else:D=corrupt;converged=True;fallback=0
     rows.append(dict(heat_id=h,band=b,seed=seed,kind=kind,strength=strength,**metrics(y,corrupt,D,truth,mask,w),converged=converged,numerical_fallback_points=fallback))
  pd.DataFrame(rows).to_csv(P/'results/exact_final_pipeline_confirmation.csv',index=False);print('final pipeline',h,b,flush=True)
if __name__=='__main__':main()
