"""Scale-calibrated multi-evidence point score; original downstream models untouched."""
from core import *
from scipy.ndimage import median_filter,label,find_objects,maximum_filter

def adaptive(y,t,w):
 good=np.isfinite(y)&(y>=0);sat=y==65535;z=fill_training(np.where(good&~sat,y,np.nan),t,w)
 tr=z.copy();wr=z.copy();alpha=(t[1:-1]-t[:-2])/(t[2:]-t[:-2]);tr[1:-1]=(1-alpha[:,None])*z[:-2]+alpha[:,None]*z[2:];beta=(w[1:-1]-w[:-2])/(w[2:]-w[:-2]);wr[:,1:-1]=(1-beta[None,:])*z[:,:-2]+beta[None,:]*z[:,2:]
 rt=z-tr;rw=z-wr
 def scale(r,axis):return np.maximum(25,1.4826*np.median(np.stack([abs(np.roll(r,k,axis=axis)) for k in [-4,-3,-2,2,3,4]]),axis=0))
 st=scale(rt,0);sw=scale(rw,1);zt=abs(rt)/st;zw=abs(rw)/sw
 initial=detect(y,t,w,window=5,threshold=6)[0];base=fill_training(np.where(initial&good&~sat,y,np.nan),t,w);low=legacy.project(base,4,y.size<100000);res=z-low;center=median_filter(res,size=(9,9),mode='nearest');sl=np.maximum(50,1.4826*median_filter(abs(res-center),size=(9,9),mode='nearest'));zl=abs(res-center)/sl
 prominence=z-median_filter(z,size=(1,5),mode='nearest');persistent=np.zeros(y.shape,np.uint8)
 for k in [-2,-1,1,2]:
  q=np.roll(prominence,k,axis=0);persistent+=((q*prominence>0)&(abs(q)>=.5*abs(prominence))).astype(np.uint8)
 consistency=1/(1+abs(tr-wr)/(np.maximum(50,.2*abs(tr))+st+sw));score=np.sqrt(zt*(.5*zw/(1+persistent)+.5*zl))*consistency
 eligible=good&~maximum_filter(sat,size=3,mode='nearest');eligible[:4]=False;eligible[-4:]=False;eligible[:,:4]=False;eligible[:,-4:]=False
 return dict(score=score,time_z=zt,wave_z=zw,lowrank_z=zl,time_scale=st,wave_scale=sw,lowrank_scale=sl,time_reference=tr,wave_reference=wr,persistence=persistent,consistency=consistency,eligible=eligible)

def local_mask(e,threshold):
 raw=(e['score']>=threshold)&e['eligible'];lab,n=label(raw,structure=np.ones((3,3)));allowed=np.zeros(n+1,bool)
 for k,sl in enumerate(find_objects(lab),1):
  if sl:allowed[k]=(sl[0].stop-sl[0].start<=3 and sl[1].stop-sl[1].start<=7)
 return raw&allowed[lab]

def all_masks(y,t,w,h):
 e=adaptive(y,t,w);mad=detect(y,t,w,**V1PAR[h]['detector'])[0];c=conservative_evidence(y,t,w)
 return {'A_P0':np.zeros(y.shape,bool),'B_V1_MAD':~mad&(y!=65535)&np.isfinite(y),'C_V31':c['high'],'D_adaptive4':local_mask(e,4),'D_adaptive6':local_mask(e,6)},e

def reconstruct(y,t,w,mask):
 valid=np.isfinite(y)&(y>=0)&(y!=65535)&~mask
 if not mask.any():return y.copy(),dict(converged=True,iterations=0,seconds=0.)
 r=solve(y,t,w,valid,max_iter=300);D,_,_=project_record(y,valid,r['latent'],t,w,True,True)
 return D,dict(converged=r['converged'],iterations=len(r['history']),seconds=r['seconds'])
