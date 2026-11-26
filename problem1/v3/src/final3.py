import json,time,shutil
import numpy as np,pandas as pd
from ctx import *
from preprocess import apply
from model3 import solve,project_record
from interface import Spectrum

def main():
 assert (P/'results/PREPROCESSING_ACCEPTED.json').exists();assert (P/'config/FINAL_MODEL_LOCKED.json').exists();policy=json.loads((P/'config/final_preprocessing_policy.json').read_text())['policy'];summary=[];origin=pd.read_csv(P/'results/negative_origin_audit.csv').to_dict('records') if (P/'results/negative_origin_audit.csv').exists() else []
 for h in HEATS:
  if (P/f'results/{h}_final.json').exists():summary.append(json.loads((P/f'results/{h}_final.json').read_text()));continue
  y,t,w=load(h);pp=apply(y,t,w,policy,heat=h)
  if h in CFG['stability_heats']:
   if not (P/f'data/{h}_native_fit.npz').exists():print('waiting native fit',h,flush=True);continue
   L=np.load(P/f'data/{h}_native_fit.npz')['latent'];solver=json.loads((P/f'results/{h}_native_solver.json').read_text())
  else:
   r=solve(y,t,w,pp['valid'],weights=pp['weights'],max_iter=300);L=r['latent'];solver={k:v for k,v in r.items() if k not in ['latent','checkpoints']};dump(P/f'results/{h}_native_solver.json',solver)
  D,q,detail=project_record(y,pp['valid'],L,t,w,True,True);q|=pp['types']<<6
  conditional=D.copy();sat=pp['sat'];conv=np.zeros(y.shape,bool)
  if sat.any():
   if not (P/f'data/{h}_conditional_saturation.npz').exists():print('waiting saturation',h,flush=True);continue
   s=np.load(P/f'data/{h}_conditional_saturation.npz');conditional[sat]=s['estimate'][sat];conv=s['solver_converged_at_sat'];q[sat]|=2;q[sat&~conv]|=8192
  model=Spectrum(t,w,D,y,q,conditional);np.savez_compressed(P/f'data/{h}_interface.npz',time=t,wavelength=w,raw=y,cleaned=np.where(pp['invalid'],np.nan,y),reconstructed=D,latent=L,conditional=conditional,valid=pp['valid'],invalid=pp['invalid'],quality=q,quality_grade=pp['grade'],suspected_saturation=sat,saturated_unidentifiable=sat,conditional_solver_converged=conv,offline=True,interpolation='bilinear',rank=4,eta=.003,raw_is_absolute_radiometry=False)
  td=np.linspace(t[0],t[-1],2*len(t)-1);wd=np.linspace(w[0],w[-1],2*len(w)-1);dense=model.evaluate(td,wd).astype(np.float32);np.savez_compressed(P/f'data/{h}_dense.npz',time=td,wavelength=wd,recorded=dense,offline=True)
  old=np.load(V2/f'data/{h}_model.npz');oldmask=np.load(V2/f'data/{h}_masks.npz')['valid'];oldgrid=old['grid'];neg=oldgrid<0
  origin.append(dict(heat_id=h,V2_negative=int(neg.sum()),V2_negative_at_valid=int((neg&oldmask).sum()),V2_negative_at_masked=int((neg&~oldmask).sum()),V3_negative=int((D<0).sum()),V3_max_raw_valid_difference=float(np.max(abs(D[pp['valid']]-y[pp['valid']])))))
  out=dict(heat_id=h,records=len(t),channels=len(w),candidate_frames=int(pp['candidate_rows'].sum()),confirmed_invalid_count=int(pp['invalid'].sum()),recorded_saturation_count=int(sat.sum()),recorded_saturation_max_error=float(np.max(abs(D[sat]-65535))) if sat.any() else 0,conditional_lower_violations=int((conditional[sat]<65535).sum()),conditional_converged_sat_fraction=float(conv[sat].mean()) if sat.any() else np.nan,valid_max_error=float(np.max(abs(D[pp['valid']]-y[pp['valid']]))),negative_count=int((D<0).sum()),dense_negative_count=int((dense<0).sum()),raw_peak_max_relative_change=float(np.max(abs(D.max(1)-y.max(1))/np.maximum(y.max(1),1))),lowrank_converged=solver['converged'],iterations=len(solver['history']),stop_reason=solver['stop_reason'],guarded_missing=int((q&8>0).sum()),dense_shape=str(dense.shape),original_cache_sha256=sha(V1/f'data/processed/{h}.npz'),interface_sha256=sha(P/f'data/{h}_interface.npz'))
  dump(P/f'results/{h}_final.json',out);summary.append(out);pd.DataFrame(summary).to_csv(P/'results/final_summary.csv',index=False);pd.DataFrame(origin).to_csv(P/'results/negative_origin_audit.csv',index=False);print('final',h,solver['converged'],len(solver['history']),flush=True)
 if len(summary)==10:dump(P/'results/FINAL_COMPLETE.json',dict(heats=HEATS,target='recorded clipped observations; conditional estimates separate',policy=policy,model='V3-C',offline=True))
if __name__=='__main__':main()
