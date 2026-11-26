"""Conservative bidirectional evidence; statistical row/column candidates are not deleted."""
import numpy as np
from scipy.ndimage import median_filter
from sklearn.decomposition import PCA
from sklearn.covariance import MinCovDet
from scipy.stats import chi2

def fill_training(y,t,w):
 """Two-axis interpolation from observed entries only. Full rows/columns handled separately."""
 y=np.asarray(y,float);a=y.copy();b=y.copy();good=np.isfinite(y)
 if not good.any():raise ValueError('No training observations')
 for j in range(y.shape[1]):
  q=good[:,j]
  if q.any():a[:,j]=np.interp(t,t[q],y[q,j])
 for i in range(y.shape[0]):
  q=good[i]
  if q.any():b[i]=np.interp(w,w[q],y[i,q])
 n=np.isfinite(a).astype(int)+np.isfinite(b).astype(int)
 z=(np.nan_to_num(a)+np.nan_to_num(b))/np.maximum(n,1)
 if np.any(n==0):z[n==0]=np.median(y[good])
 z[good]=y[good]
 return z

def detect(y,t,w,window=5,threshold=8.0):
 good=np.isfinite(y);sat=good&(y==65535);base=fill_training(np.where(good&~sat,y,np.nan),t,w)
 mt=median_filter(base,size=(window,1),mode='nearest');mw=median_filter(base,size=(1,window),mode='nearest')
 rt=base-mt;rw=base-mw
 st=1.4826*median_filter(abs(rt),size=(window,1),mode='nearest')
 sw=1.4826*median_filter(abs(rw),size=(1,window),mode='nearest')
 # Noise floor prevents quantized/flat neighborhoods producing unlimited z-scores.
 floor=max(25.0,float(np.median(abs(rw)))*1.4826*.25)
 zt=abs(rt)/np.maximum(st,floor);zw=abs(rw)/np.maximum(sw,floor)
 local=good&~sat&(zt>threshold)&(zw>threshold)
 types=np.zeros(y.shape,np.uint8);types[~good]|=1;types[sat]|=2;types[local]|=4
 valid=(types&7)==0
 return valid,types,dict(time_z=zt,wave_z=zw,noise_floor=floor)

def spectrum_candidates(y,seed=20261126):
 z=np.log1p(np.maximum(y,0));z=z[:,::4]
 scores=PCA(n_components=6,svd_solver='randomized',random_state=seed).fit_transform(z)
 mcd=MinCovDet(random_state=seed,support_fraction=.75).fit(scores)
 distance=mcd.mahalanobis(scores)
 return distance,distance>chi2.ppf(.997,6)

def channel_candidates(y):
 mean=np.median(y,axis=0);local=median_filter(mean,size=11,mode='nearest');r=mean-local
 s=max(25,1.4826*np.median(abs(r-np.median(r))))
 return abs(r)/s,abs(r)>8*s

def inject(y,kind,rng):
 z=y.copy();truth=np.zeros(y.shape,bool);n,m=y.shape
 if kind in ['spike','drop']:
  count=max(30,int(y.size*.001));idx=rng.choice(y.size,count,replace=False);truth.ravel()[idx]=True
  z.ravel()[idx]+=rng.choice([-1,1],count)*30000 if kind=='spike' else -np.maximum(y.ravel()[idx]*.9,1500)
 elif kind=='short_block':
  for k in range(8):
   i=int(rng.integers(3,n-5));j=int(rng.integers(10,m-15));truth[i:i+3,j:j+7]=True;z[i:i+3,j:j+7]+=25000
 elif kind=='row':
  ix=rng.choice(np.arange(3,n-3),2,replace=False);truth[ix]=True;z[ix]+=25000
 elif kind=='column':
  ix=rng.choice(np.arange(5,m-5),2,replace=False);truth[:,ix]=True;z[:,ix]+=25000
 elif kind=='saturation':
  for k in range(4):
   i=int(rng.integers(3,n-8));j=int(rng.integers(10,m-20));truth[i:i+4,j:j+12]=True;z[i:i+4,j:j+12]=65535
 return z,truth

def detection_metrics(pred,truth,eligible=None):
 q=np.ones(truth.shape,bool) if eligible is None else eligible
 tp=int((pred&truth&q).sum());fp=int((pred&~truth&q).sum());fn=int((~pred&truth&q).sum());p=tp/max(tp+fp,1);r=tp/max(tp+fn,1)
 return dict(precision=p,recall=r,f1=2*p*r/max(p+r,1e-12),tp=tp,fp=fp,fn=fn)
