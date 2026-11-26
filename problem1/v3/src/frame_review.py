"""Leave-current-frame-out native-coordinate diagnosis; no statistical fault confirmation."""
import numpy as np,pandas as pd
from scipy.ndimage import label
from sklearn.utils.extmath import randomized_svd
from ctx import *

def robust_z(v):
 med=np.median(v);s=max(1e-9,1.4826*np.median(abs(v-med)));return abs(v-med)/s,med,s

def frame_metrics(y,t,w):
 finite=np.isfinite(y)&(y>=0);z=fill_training(np.where(finite,y,np.nan),t,w);n,m=z.shape
 area=np.trapezoid(z,w,axis=1);norm=np.linalg.norm(z,axis=1);ref=np.empty_like(z);coverage=np.zeros(n);recovery=np.zeros(n)
 for i in range(n):
  ix=np.array([k for k in [i-2,i-1,i+1,i+2] if 0<=k<n]);assert i not in ix
  # Actual temporal distance; weighted average excludes the current row by construction.
  weights=1/np.maximum(abs(t[ix]-t[i]),1e-12);weights/=weights.sum();ref[i]=weights@z[ix]
  spread=np.maximum(25,1.4826*np.median(abs(z[ix]-np.median(z[ix],axis=0)),axis=0));coverage[i]=np.mean(abs(z[i]-ref[i])>6*spread)
  if 0<i<n-1:recovery[i]=np.linalg.norm(z[i-1]-z[i+1])/max(np.linalg.norm(ref[i]),1)
 refarea=np.trapezoid(ref,w,axis=1);logchange=np.log(np.maximum(area,1)/np.maximum(refarea,1));brightness=np.median((z-ref)/np.maximum(ref,25),axis=1)
 angle=np.degrees(np.arccos(np.clip(np.sum(z*ref,axis=1)/np.maximum(norm*np.linalg.norm(ref,axis=1),1e-12),-1,1)))
 shape=z/np.maximum(norm[:,None],1);refs=ref/np.maximum(np.linalg.norm(ref,axis=1)[:,None],1);shaperes=np.linalg.norm(shape-refs,axis=1)
 # Cross-fitted spectral basis: the current frame never appears in the basis training set.
 lr=np.zeros(n)
 for parity in [0,1]:
  train=np.arange(n)%2!=parity;test=~train;r=min(4,int(train.sum()),m)
  _,_,v=randomized_svd(z[train],n_components=r,n_iter=5,random_state=137);pred=(z[test]@v.T)@v;lr[test]=np.linalg.norm(z[test]-pred,axis=1)/np.maximum(norm[test],1)
 zint,_,sint=robust_z(logchange);zangle,ma,sa=robust_z(angle);zs,_,_=robust_z(shaperes);zl,_,_=robust_z(lr)
 flat=np.ptp(z,axis=1)==0;invalid=(~finite).all(1)
 candidate=(abs(logchange)>max(np.log1p(.2),6*sint))|(angle>max(5,ma+6*sa))|(zs>6)|(zl>6)|flat|invalid
 # Fault confirmation requires external/structural impossibility evidence; statistical shape alone does not suffice.
 structural=invalid|((z<=0).all(1))
 events,count=label(candidate);dur=np.bincount(events)[events]*np.median(np.diff(t));dur[~candidate]=0
 kind=np.where(flat,'flat_spectrum_candidate',np.where((abs(brightness)>.2)&(angle<2),'brightness_change_candidate',np.where(angle>5,'shape_change_candidate','multimetric_candidate')))
 kind=np.where(~candidate,'valid_or_unremarkable',kind);kind=np.where(structural,'structurally_unavailable',kind)
 seq=np.rint(t/.5).astype(int)+1
 return pd.DataFrame(dict(time_s=t,sequence=seq,integral_intensity=area,reference_integral=refarea,integral_first_difference_per_second=np.r_[np.nan,np.diff(area)/np.diff(t)],integral_relative_first_difference=np.r_[np.nan,np.diff(area)/np.maximum(area[:-1],1)],integral_relative_change=np.expm1(logchange),local_change_per_second=np.expm1(logchange)/np.median(np.diff(t)),brightness_relative=brightness,spectral_angle_deg=angle,normalized_shape_residual=shaperes,time_outlier_wavelength_fraction=coverage,event_id=events,event_duration_s=dur,neighbor_recovery_relative=recovery,saturation_fraction=(y==65535).mean(1),crossfit_lowrank_residual=lr,integral_robust_z=zint,angle_robust_z=zangle,shape_robust_z=zs,lowrank_robust_z=zl,candidate=candidate,structural_invalid=structural,anomaly_type=kind,confidence=np.where(structural,'high_structural_invalid',np.where(candidate,'suspected_unconfirmed','retained')),basis='leave-current-frame-out temporal reference; opposite-parity low-rank basis',recommendation=np.where(structural,'mask_invalid_extent',np.where(candidate,'retain_flag_and_sensitivity','retain')),confirmed_acquisition_fault=False))

def main():
 rows=[];ch=[];sweep=[]
 for h in HEATS:
  y,t,w=load(h);d=frame_metrics(y,t,w);d.insert(0,'heat_id',h);rows.append(d);d.to_csv(P/f'results/{h}_frame_diagnostics.csv',index=False)
  median=np.median(y,axis=0);from scipy.ndimage import median_filter
  contrast=median-median_filter(median,size=11,mode='nearest');scale=max(25,1.4826*np.median(abs(contrast)));persistent=abs(contrast)>8*scale
  ch.append(pd.DataFrame(dict(heat_id=h,wavelength_nm=w,median_intensity=median,contrast=contrast,score=abs(contrast)/scale,candidate=persistent,confidence=np.where(persistent,'suspected_true_line_or_channel','retained'),reason='persistent spectral contrast cannot distinguish real emission line from channel fault',recommendation='retain_and_flag',confirmed_instrument_fault=False)))
  for b in [.1,.2,.3]:
   for a in [1,2,5,10]:sweep.append(dict(heat_id=h,brightness_threshold=b,angle_threshold_deg=a,brightness_candidates=int((abs(d.brightness_relative)>b).sum()),angle_candidates=int((d.spectral_angle_deg>a).sum()),joint_candidates=int(((abs(d.brightness_relative)>b)&(d.spectral_angle_deg>a)).sum()),deletion_validated=False))
  print('frame review',h,int(d.candidate.sum()),flush=True)
 pd.concat(rows).to_csv(P/'results/frame_diagnostics_all.csv',index=False);pd.concat(rows).query('candidate').to_csv(P/'results/frame_anomaly_candidates.csv',index=False);pd.concat(ch).to_csv(P/'results/channel_diagnostics_all.csv',index=False);pd.concat(ch).query('candidate').to_csv(P/'results/channel_anomaly_candidates.csv',index=False);pd.DataFrame(sweep).to_csv(P/'results/frame_threshold_sensitivity.csv',index=False)
if __name__=='__main__':main()
