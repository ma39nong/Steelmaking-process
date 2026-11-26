"""Fixed post-lock cap check and predeclared numerical fallback audit; no parameter selection."""
from ctx import *
from preprocess import apply
from model3 import solve,project_record
import pandas as pd

def main():
 caps=[];fallback=[]
 for h,stage,y,t,w,i,j in panels('model'):
  if stage!=1:continue
  # One prespecified threshold/heat. Reference is genuinely below the instrument cap.
  eligible=np.isfinite(y)&(y>0)&(y<64000);cap=float(np.quantile(y[eligible],.8));sat=eligible&(y>=cap);visible=y.copy();visible[sat]=cap;valid=eligible&~sat
  r=convex.fit_convex(visible,t,w,valid,sat,cap,tau=.01,rho=1,max_iter=300,tol=1e-4)
  baseline=fill_training(np.where(valid,y,np.nan),t,w)
  caps.append(dict(heat_id=h,cap=cap,count=int(sat.sum()),convex_rmse=float(np.sqrt(np.mean((r['latent'][sat]-y[sat])**2))),mask_linear_rmse=float(np.sqrt(np.mean((baseline[sat]-y[sat])**2))),lower_violation_fraction=float(np.mean(r['latent'][sat]<cap)),converged=r['converged'],iterations=len(r['history'])))
  if h not in CFG['stability_heats']:continue
  hidden=dict(hides(y.shape,9307))['random_0.05'];v=y.copy();v[hidden]=np.nan;pp=apply(v,t,w,'P0');s=solve(v,t,w,pp['valid'],max_iter=1200,checkpoints=[600,1200],wall_seconds=120)
  a=project_record(v,pp['valid'],s['checkpoints'][600],t,w,True,True)[0];b=project_record(v,pp['valid'],s['checkpoints'][1200],t,w,True,True)[0]
  delta=abs(a-b);unstable=hidden&(delta>CFG['stability_output_rms_tolerance_counts']);linear=fill_training(np.where(pp['valid'],v,np.nan),t,w);af=a.copy();bf=b.copy();af[unstable]=linear[unstable];bf[unstable]=linear[unstable]
  # The SAME visible-only fallback is used at both checkpoints, never selected by true error.
  np.savez_compressed(P/f'data/{h}_numerical_fallback_audit.npz',mask=unstable,prediction=bf,hidden=hidden)
  fallback.append(dict(heat_id=h,unstable_hidden_count=int(unstable.sum()),hidden_count=int(hidden.sum()),raw_delta_rms=float(np.sqrt(np.mean((a[hidden]-b[hidden])**2))),protected_delta_rms=float(np.sqrt(np.mean((af[hidden]-bf[hidden])**2))),protected_delta_q95=float(np.quantile(abs(af[hidden]-bf[hidden]),.95)),protected_rmse=score(y,bf,hidden,w)['rmse'],raw_1200_rmse=score(y,b,hidden,w)['rmse'],strict_converged=s['converged']))
 pd.DataFrame(caps).to_csv(P/'results/fixed_artificial_cap_check.csv',index=False);pd.DataFrame(fallback).to_csv(P/'results/numerical_fallback_audit.csv',index=False)
if __name__=='__main__':main()
