"""Conservative detector; auxiliary row/channel evidence is exported separately."""
import warnings
import numpy as np
from scipy.ndimage import median_filter,label
from common import detect,fill_training
from solver import project

def diagnostics(y,t,w,threshold=9,aux_threshold=14):
 good=np.isfinite(y);valid,types,base=detect(y,t,w,window=5,threshold=threshold)
 z=fill_training(np.where(valid,y,np.nan),t,w)
 low=project(z,min(4,min(z.shape)-1),z.size<100000)
 residual=z-low;center=np.median(residual,axis=0);scale=np.maximum(50,1.4826*np.median(abs(residual-center),axis=0));rz=abs(residual-center)/scale
 norm=np.linalg.norm(z,axis=1);shape=z/np.maximum(norm[:,None],1)
 angle=np.zeros(len(t));angle[1:]=np.degrees(np.arccos(np.clip(np.sum(shape[1:]*shape[:-1],axis=1),-1,1)))
 brightness=np.median(z,axis=1);reference=median_filter(brightness,size=5,mode='nearest');br=abs(brightness-reference)/np.maximum(reference,50)
 # Long/narrow features remain diagnostic. Automatic additions require two-axis evidence.
 support=good&(rz>aux_threshold)&(base['time_z']>threshold*.7)&(base['wave_z']>threshold*.7)
 components,count=label(support);sizes=np.bincount(components.ravel());connected=support&(sizes[components]>=2)
 added=connected&valid;valid=valid&~added;types[added]|=32
 channel_score=np.median(rz,axis=0);frame_score=np.quantile(rz,.95,axis=1)
 return valid,types,dict(residual_z=rz,component_labels=components,frame_angle_deg=angle,frame_brightness_relative=br,frame_residual=frame_score,channel_residual=channel_score,added=added)

def inject(y,kind,rng,strength,length):
 z=y.copy();truth=np.zeros(y.shape,bool);n,m=y.shape;lo=max(float(np.quantile(y,.9)-np.quantile(y,.1)),1000);amp=lo*strength
 if kind in ['spike','drop']:
  ix=rng.choice(y.size,max(8,int(y.size*.002)),replace=False);truth.ravel()[ix]=True;z.ravel()[ix]+=amp if kind=='spike' else -amp
 elif kind=='short_block':
  i=int(rng.integers(5,n-length-5));j=int(rng.integers(5,m-3*length-5));truth[i:i+length,j:j+3*length]=True;z[truth]+=amp
 elif kind=='row':
  i=int(rng.integers(5,n-length-5));truth[i:i+length]=True;z[truth]+=amp
 elif kind=='column':
  j=int(rng.integers(5,m-length-5));truth[:,j:j+length]=True;z[truth]+=amp
 elif kind=='saturation':
  i=int(rng.integers(5,n-length-5));j=int(rng.integers(5,m-3*length-5));truth[i:i+length,j:j+3*length]=True;z[truth]=65535
 return z,truth

def detection_score(pred,truth,baseline):
 # Incremental detections on pre-existing anomalous observations are not counted as clean false positives.
 eligible=~baseline|truth;tp=int((pred&truth).sum());fp=int((pred&~truth&eligible).sum());fn=int((~pred&truth).sum());tn=int((~pred&~truth&eligible).sum())
 p=tp/max(tp+fp,1);r=tp/max(tp+fn,1)
 return dict(precision=p,recall=r,f1=2*p*r/max(p+r,1e-12),false_positive_rate=fp/max(fp+tn,1),tp=tp,fp=fp,fn=fn,tn=tn,total_deleted_fraction=float(pred.mean()),incremental_reference_exclusion=float(baseline.mean()))
