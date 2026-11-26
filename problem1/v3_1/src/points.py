"""Fixed multi-evidence statistical point grading; never confirms hardware faults."""
from common31 import *
from scipy.ndimage import median_filter,label,maximum_filter

def evidence(y,t,w):
 good=np.isfinite(y)&(y>=0);sat=y==65535
 valid,types,aux=v2_detect(y,t,w,threshold=6,aux_threshold=10)
 _,_,mad=detect(y,t,w,window=5,threshold=6)
 candidate=good&~sat&~valid
 # Fit a rank4 spectral background to candidate-masked data, then score ORIGINAL values.
 z=fill_training(np.where(good&valid&~sat,y,np.nan),t,w);low=legacy.project(z,4,y.size<100000)
 res=y-low;center=np.nanmedian(np.where(good&~candidate&~sat,res,np.nan),axis=0);scale=np.maximum(50,1.4826*np.nanmedian(abs(np.where(good&~candidate&~sat,res-center,np.nan)),axis=0));rz=abs(res-center)/scale
 # Leave-center axial neighbors; only interior cells are eligible for automatic repair.
 axial=np.stack([np.roll(y,k,axis=a) for a in [0,1] for k in [-2,-1,1,2]])
 tr=np.median(axial[:4],axis=0);wr=np.median(axial[4:],axis=0);floor=mad['noise_floor'];c=CFG['conservative']
 agree=abs(tr-wr)<=np.maximum(c['neighbor_agreement_noise_multiple']*floor,c['neighbor_agreement_relative']*np.maximum(abs((tr+wr)/2),100))
 signs=(y-tr)*(y-wr)>0
 exterior=(y>np.max(axial,axis=0)+c['outside_neighbor_noise_multiple']*floor)|(y<np.min(axial,axis=0)-c['outside_neighbor_noise_multiple']*floor)
 prominence=y-median_filter(y,size=(1,5),mode='nearest');persistent=np.zeros(y.shape,np.uint8)
 for k in [-2,-1,1,2]:
  p=np.roll(prominence,k,axis=0);persistent+=((p*prominence>0)&(abs(p)>=c['persistent_prominence_ratio']*abs(prominence))).astype(np.uint8)
 components,n=label(candidate,structure=np.ones((3,3)));sizes=np.bincount(components.ravel());small=np.zeros(n+1,bool)
 from scipy.ndimage import find_objects
 for k,sl in enumerate(find_objects(components),1):
  if sl:small[k]=sizes[k]<=c['component_max_cells'] and sl[0].stop-sl[0].start<=c['component_max_time_span'] and sl[1].stop-sl[1].start<=c['component_max_wave_span']
 near_sat=maximum_filter(sat,size=5,mode='nearest');interior=np.ones(y.shape,bool);interior[:2]=False;interior[-2:]=False;interior[:,:2]=False;interior[:,-2:]=False
 high=candidate&(mad['time_z']>=c['two_axis_z'])&(mad['wave_z']>=c['two_axis_z'])&(rz>=c['lowrank_residual_z'])&agree&signs&exterior&(persistent<c['persistent_neighbor_frames'])&small[components]&~near_sat&interior
 high[high.mean(1)>c['max_mask_fraction_per_frame']]=False
 grade=np.zeros(y.shape,np.uint8);grade[candidate]=1;grade[high]=2
 return dict(candidate=candidate,high=high,suspect=candidate&~high,grade=grade,time_z=mad['time_z'],wave_z=mad['wave_z'],residual_z=rz,temporal_reference=tr,wave_reference=wr,neighbor_agreement=agree,exterior=exterior,persistence=persistent,component=components,component_size=sizes[components],near_saturation=near_sat,lowrank_background=low)

def mask_for(e,name):return np.zeros(e['candidate'].shape,bool) if name=='V3-P0' else e['high'] if name=='conservative' else e['candidate']
def repair(y,t,w,e,name,budget=300):
 remove=mask_for(e,name);valid=np.isfinite(y)&(y>=0)&(y!=65535)&~remove
 if not remove.any():return y.copy(),dict(converged=True,iterations=0,seconds=0.,fit_required=False),np.zeros(y.shape,np.uint16)
 r=solve(y,t,w,valid,rank=4,eta=.003,max_iter=budget);D,q,_=project_record(y,valid,r['latent'],t,w,True,True)
 return D,dict(converged=r['converged'],iterations=len(r['history']),seconds=r['seconds'],fit_required=True),q
