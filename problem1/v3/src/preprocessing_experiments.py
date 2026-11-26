"""Finite strategy and synthetic diagnosis validation; all thresholds fixed in protocol."""
import sys,json
import numpy as np,pandas as pd
from scipy.ndimage import label
from ctx import *
from frame_review import frame_metrics
from preprocess import apply
from model3 import solve,project_record

def binary(pred,truth):
 tp=int((pred&truth).sum());fp=int((pred&~truth).sum());fn=int((~pred&truth).sum());tn=int((~pred&~truth).sum());p=tp/max(tp+fp,1);r=tp/max(tp+fn,1)
 return dict(precision=p,recall=r,f1=2*p*r/max(p+r,1e-12),false_positive_rate=fp/max(fp+tn,1),tp=tp,fp=fp,fn=fn,tn=tn)
def event_metrics(pred,truth):
 a,na=label(pred);b,nb=label(truth);matched_p={int(x) for x in np.unique(a[truth]) if x};matched_b={int(x) for x in np.unique(b[pred]) if x};p=len(matched_p)/max(na,1);r=len(matched_b)/max(nb,1)
 return dict(precision=p,recall=r,f1=2*p*r/max(p+r,1e-12),false_positive_rate=(na-len(matched_p))/max(len(pred),1),predicted_events=na,true_events=nb)
def inject(y,kind,rng,amp,dur):
 z=y.copy();truth=np.zeros(y.shape,bool);n,m=y.shape;i=int(rng.integers(5,n-dur-5));truth[i:i+dur]=True;span=max(float(np.quantile(y,.9)),100)
 if kind=='gain':z[truth]*=1+amp
 elif kind=='offset':z[truth]+=amp*span
 elif kind=='attenuation':z[truth]*=1-amp
 elif kind=='flat':z[i:i+dur]=np.median(y[i:i+dur],axis=1)[:,None]
 elif kind=='multiframe':z[truth]*=1+amp
 elif kind=='local_wave':
  truth[:]=False;j=int(rng.integers(10,m-25));truth[i:i+dur,j:j+16]=True;z[truth]+=amp*span
 return np.clip(z,0,65535),truth

def main(split):
 rows=[];injections=[];risk=[];fidelity=[]
 seeds=CFG['preprocessing_seeds'] if split=='development' else CFG['confirmation_seeds']
 if split!='development':assert (P/'config/final_preprocessing_policy.json').exists()
 for h,stage,y,t,w,i,j in panels('preprocessing',split):
  base=frame_metrics(y,t,w);fast=abs(base.integral_relative_change)>np.quantile(abs(base.integral_relative_change),.9)
  risk.append(dict(heat_id=h,split=split,fast_frames=int(fast.sum()),fast_frame_candidates=int((fast&base.candidate).sum()),fast_frame_auto_deleted=int((fast&base.structural_invalid).sum()),all_frame_candidates=int(base.candidate.sum()),confirmed_real_faults=0))
  for seed in seeds:
   for kind in CFG['synthetic_types']:
    for amp,dur in zip(CFG['synthetic_amplitudes'],CFG['synthetic_durations']):
     rng=np.random.default_rng(seed+int(h[-3:]));z,truth=inject(y,kind,rng,amp,dur);fm=frame_metrics(z,t,w);point=~v2_detect(z,t,w,threshold=6,aux_threshold=10)[0];frame=fm.candidate.to_numpy()
     # Candidate localization is diagnostic; confirmed-mask localization remains a separate stricter outcome.
     combined=point.copy();combined[frame]=True;frtruth=truth.any(1)
     for level,sc in [('point',binary(combined,truth)),('frame',binary(frame,frtruth)),('event',event_metrics(frame,frtruth))]:injections.append(dict(heat_id=h,split=split,seed=seed,kind=kind,amplitude=amp,duration=dur,level=level,region_iou=float((combined&truth).sum()/max((combined|truth).sum(),1)),location_error_cells=int((combined^truth).sum()),auto_masked_frame_count=int(fm.structural_invalid.sum()),**sc))
   for scenario,hidden in hides(y.shape,seed):
    visible=y.copy();visible[hidden]=np.nan
    for policy in CFG['policies']:
     pp=apply(visible,t,w,policy);r=solve(visible,t,w,pp['valid'],weights=pp['weights'],max_iter=100);D,_,_=project_record(visible,pp['valid'],r['latent'],t,w,True,True)
     rows.append(dict(heat_id=h,split=split,seed=seed,scenario=scenario,policy=policy,masked_fraction=float(pp['invalid'].mean()),suspected_frames=int(pp['candidate_rows'].sum()),auto_full_rows=int((pp['invalid'].all(1)&np.isfinite(visible).any(1)).sum()),seconds=r['seconds'],**score(y,D,hidden,w)))
   # Real candidate full-row masking is sensitivity, scored separately from random missingness.
   for policy in ['P2','P3','P4']:
    pp=apply(y,t,w,policy);hidden=np.repeat(pp['candidate_rows'][:,None],y.shape[1],axis=1)
    if policy=='P4' and hidden.any():
     r=solve(y,t,w,pp['valid'],weights=pp['weights'],max_iter=100);D,_,_=project_record(y,pp['valid'],r['latent'],t,w,True,True);rows.append(dict(heat_id=h,split=split,seed=seed,scenario='suspect_frame_mask_sensitivity',policy=policy,masked_fraction=float(pp['invalid'].mean()),suspected_frames=int(pp['candidate_rows'].sum()),auto_full_rows=int(pp['invalid'].all(1).sum()),seconds=r['seconds'],**score(y,D,hidden,w)))
  for policy in CFG['policies']:
   pp=apply(y,t,w,policy)
   if not pp['valid'].any():
    fidelity.append(dict(heat_id=h,policy=policy,feasible=False,raw_peak_p95=np.nan));continue
   rr=solve(y,t,w,pp['valid'],weights=pp['weights'],max_iter=100);D,_,_=project_record(y,pp['valid'],rr['latent'],t,w,True,True)
   q=y.max(1)<65535;rel=abs(D.max(1)-y.max(1))/np.maximum(y.max(1),1);area=np.trapezoid(y,w,axis=1);da=np.trapezoid(D,w,axis=1)
   fidelity.append(dict(heat_id=h,policy=policy,feasible=True,raw_peak_p95=float(np.quantile(rel[q],.95)) if q.any() else np.nan,integral_relative_mae=float(np.mean(abs(da-area)/np.maximum(area,1))),temporal_integral_difference_rmse=float(np.sqrt(np.mean((np.diff(da)-np.diff(area))**2))),masked_fraction=float(pp['invalid'].mean())))
  print('preprocess',split,h,flush=True)
 pd.DataFrame(fidelity).to_csv(P/f'results/preprocessing_signal_fidelity_{split}.csv',index=False)
 pd.DataFrame(rows).to_csv(P/f'results/preprocessing_{split}.csv',index=False);pd.DataFrame(injections).to_csv(P/f'results/synthetic_diagnosis_{split}.csv',index=False);pd.DataFrame(risk).to_csv(P/f'results/fast_change_risk_{split}.csv',index=False)
 if split=='development':
  df=pd.DataFrame(rows);q=df[df.scenario!='suspect_frame_mask_sensitivity'];means=q.groupby(['heat_id','scenario','policy']).nrmse.mean().groupby('policy').mean();baseline=means['P0'];eligible=[]
  for p in ['P2','P3']:
   a=q[q.policy==p].groupby('scenario').nrmse.mean();b=q[q.policy=='P0'].groupby('scenario').nrmse.mean()
   if means[p]<=.95*baseline and (a<=1.1*b).all() and q[q.policy==p].masked_fraction.max()<=.01 and pd.DataFrame(fidelity).query('policy==@p').raw_peak_p95.max()<=.01:eligible.append(p)
  selected=min(eligible,key=lambda p:means[p]) if eligible else 'P0'
  dump(P/'config/final_preprocessing_policy.json',dict(policy=selected,macro_nrmse=means.to_dict(),rank=4,eta=.003,whole_frame_faults_confirmed=0,automatic_statistical_frame_deletion=False,selection='Predeclared development prediction/risk rule; no test F1 tuning',protocol_sha256=sha(P/'config/protocol.json'),scope='offline; no online preprocessing claim'))
if __name__=='__main__':main(sys.argv[1])
