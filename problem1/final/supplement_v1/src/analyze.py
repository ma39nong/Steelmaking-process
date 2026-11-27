"""Supplementary diagnostics on immutable recorded and reconstructed spectra."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from scipy.ndimage import maximum_filter,median_filter
S=Path(__file__).resolve().parents[1];F=S.parent;ROOT=F.parents[1]
sys.path.insert(0,str(ROOT/'problem1/v3/src'))
from ctx import fill_training
PROTOCOL=json.loads((S/'protocol.json').read_text())
HEATS=pd.read_csv(F/'results/final_summary.csv').heat_id.tolist()

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for c in iter(lambda:f.read(1024**2),b''):h.update(c)
 return h.hexdigest()

def robust_fit(x,y):
 center=np.median(x);spread=max(float(np.quantile(x,.75)-np.quantile(x,.25)),1.);z=(x-center)/spread;A=np.column_stack([np.ones(len(x)),z]);coef=np.linalg.lstsq(A,y,rcond=None)[0]
 for _ in range(20):
  e=y-A@coef;scale=max(1.,1.4826*np.median(abs(e-np.median(e))));weights=np.minimum(1,1.345*scale/np.maximum(abs(e),1e-12));new=np.linalg.lstsq(A*np.sqrt(weights[:,None]),y*np.sqrt(weights),rcond=None)[0]
  if np.linalg.norm(new-coef)<1e-6*max(np.linalg.norm(coef),1):coef=new;break
  coef=new
 residual=y-A@coef
 return coef[0]-coef[1]*center/spread,coef[1]/spread,float(1.4826*np.median(abs(residual-np.median(residual))))

def boot(x,y,seed):
 a,b,s=robust_fit(x,y);rng=np.random.default_rng(seed);slopes=[];block=16;n=len(x)
 for _ in range(64):
  starts=rng.integers(0,n-block+1,int(np.ceil(n/block)));idx=(starts[:,None]+np.arange(block)).ravel()[:n];slopes.append(robust_fit(x[idx],y[idx])[1])
 lo,hi=np.quantile(slopes,[.025,.975]);return a,b,s,float(lo),float(hi)

def channel_reference(y,w):
 r=y.copy();b=(w[1:-1]-w[:-2])/(w[2:]-w[:-2]);r[:,1:-1]=(1-b)*y[:,:-2]+b*y[:,2:];d=y-r
 scale=np.maximum(25,1.4826*np.median(np.stack([abs(np.roll(d,k,axis=1)) for k in [-4,-3,-2,2,3,4]]),axis=0));return r,d,scale

def neighbors(y,obs,i,j):
 vals=[]
 for v,q,k in [(y[:,j],obs[:,j],i),(y[i],obs[i],j)]:
  ix=np.flatnonzero(q);p=np.searchsorted(ix,k)
  if p:vals.append(v[ix[p-1]])
  if p<len(ix):vals.append(v[ix[p]])
 return max(0,float(min(vals))),float(max(vals))

def main():
 hashes=[]
 for r in pd.read_csv(ROOT/'data_manifest.csv').itertuples():
  p=ROOT/r.source_file;digest=sha(p);assert digest==r.sha256;hashes.append(dict(kind='raw',path=r.source_file,sha256=digest))
 frozen=pd.read_csv(F/'results/artifact_manifest.csv')
 for r in frozen[frozen.artifact.isin(['interface','difference'])].itertuples():
  p=F/'data'/f'{r.heat_id}_{r.artifact}.npz';digest=sha(p);assert digest==r.sha256;hashes.append(dict(kind='frozen',path=p.relative_to(ROOT),sha256=digest))
 pd.DataFrame(hashes).to_csv(S/'results/input_hashes.csv',index=False)
 allchannels=[];follow=[];repair_rows=[];svrows=[];summary=[];counterfactual=[];phys=[];groupcurves=[];targetsnap=[]
 repairtype=pd.read_csv(F/'results/statistical_repairs_evidence.csv').set_index(['heat_id','row','column'])
 for hi,h in enumerate(HEATS):
  with np.load(F/f'data/{h}_interface.npz') as z:
   y=z['raw'];D=z['reconstructed'];w=z['wavelength'];t=z['time'];mask=z['statistical_repair'];fall=z['numerical_fallback'];latent=z['latent'];sat=z['suspected_saturation']
  med=np.median(y,axis=0);q10,q90=np.quantile(y,[.1,.9],axis=0);integral=np.trapezoid(y,w,axis=1);rate=np.r_[0,abs(np.diff(integral)/np.diff(t))];fast=rate>=np.quantile(rate,.85);bright=integral>=np.quantile(integral,.85)
  ref,d,scale=channel_reference(y,w);zscore=abs(d)/scale;md=np.median(d,axis=0);mz=np.median(zscore,axis=0);pos=(d>0).mean(0);neg=(d<0).mean(0);persist=np.maximum(pos,neg)
  support=(np.sign(d[:,1:-1])*np.sign(d[:,:-2])>0).mean(0)+(np.sign(d[:,1:-1])*np.sign(d[:,2:])>0).mean(0);support=np.r_[0,support,0]/2
  candidate=(persist>=.95)&(abs(md)>=250)&(mz>=4);candidate[:4]=False;candidate[-4:]=False
  focus=383.293 if h.startswith('6C') else 835.887;jfocus=int(np.argmin(abs(w-focus)))
  ix=np.flatnonzero(candidate);top=ix[np.argsort(mz[ix])[-8:]];physical=[588.995,589.5924,766.48991,769.89645,760.0]
  selected=np.unique(np.r_[top,jfocus,[np.argmin(abs(w-v)) for v in physical]]).astype(int)
  seg=.5*(y[:,:-1]+y[:,1:])*np.diff(w);prefix=np.column_stack([np.zeros(len(t)),np.cumsum(seg,axis=1)])
  for j in range(len(w)):
   allchannels.append(dict(heat_id=h,group=h[:2],channel=j,wavelength_nm=w[j],median_intensity=med[j],median_residual=md[j],median_abs_z=mz[j],positive_fraction=pos[j],negative_fraction=neg[j],sign_persistence=persist[j],adjacent_sign_support=support[j],candidate=bool(candidate[j]),target_channel=bool(j==jfocus)))
  for j in selected:
   a=max(0,int(np.searchsorted(w,w[j]-2)));b=min(len(w)-1,int(np.searchsorted(w,w[j]+2)))
   B=(integral-(prefix[:,b]-prefix[:,a]))/(w[-1]-w[0]-(w[b]-w[a]));alpha,beta,sigma,lo,up=boot(B,d[:,j],27011+hi*100+j)
   spread=max(float(np.quantile(B,.75)-np.quantile(B,.25)),1.)
   for fraction in [.1,.15,.2]:
    cut=int(np.ceil(len(t)*fraction));order=np.argsort(B);dark=order[:cut];light=order[-cut:]
    follow.append(dict(heat_id=h,channel=j,wavelength_nm=w[j],candidate=bool(candidate[j]),focus=bool(j==jfocus),fraction=fraction,alpha=alpha,beta=beta,beta_ci_low=lo,beta_ci_high=up,residual_mad=sigma,brightness_iqr=spread,slope_iqr_counts=beta*spread,dark_residual=np.median(d[dark,j]),bright_residual=np.median(d[light,j]),dark_positive_fraction=np.mean(d[dark,j]>0),bright_positive_fraction=np.mean(d[light,j]>0),median_residual=md[j],positive_fraction=pos[j],median_abs_z=mz[j],adjacent_sign_support=support[j]))
    if j==jfocus:
     a2=max(0,j-12);b2=min(len(w),j+13)
     for k in range(a2,b2):targetsnap.append(dict(heat_id=h,fraction=fraction,wavelength_nm=w[k],dark=np.median(y[dark,k]),bright=np.median(y[light,k]),reference_dark=np.median(ref[dark,k]),reference_bright=np.median(ref[light,k])))
   if j==jfocus:
    pd.DataFrame(dict(time_s=t,raw=y[:,j],neighbor_reference=ref[:,j],residual=d[:,j],brightness=B)).to_csv(S/f'results/{h}_focus_trace.csv',index=False)
  # Native candidate is explicitly a sensitivity output, not a replacement release.
  C=D.copy();C[:,jfocus]=np.clip(ref[:,jfocus],0,65535);C[sat]=65535
  np.savez_compressed(S/f'data/{h}_channel_sensitivity.npz',time=t,wavelength=w,candidate=C,changed=abs(C-D)>1e-8,formal_release=False)
  changed=abs(C-D)>1e-8;areaD=np.trapezoid(D,w,axis=1);areaC=np.trapezoid(C,w,axis=1)
  counterfactual.append(dict(heat_id=h,channel=jfocus,wavelength_nm=w[jfocus],changed=int(changed.sum()),overlap_with_baseline_repair=int((changed&mask).sum()),integral_relative_mae=float(np.mean(abs(areaC-areaD)/np.maximum(abs(areaD),1))),peak_relative_mean=float(np.mean(abs(C.max(1)-D.max(1))/np.maximum(D.max(1),1))),formal_release=False))
  np.savez_compressed(S/f'data/{h}_summary.npz',time=t,wavelength=w,median=med,q10=q10,q90=q90,integral=integral,brightness=np.median(y,axis=1),change_rate=rate,focus_channel=jfocus)
  # Equal-heat summaries on an explicitly shared grid.
  grid=np.arange(340.535,1023.303+.0001,.5)
  for k,x in enumerate(grid):groupcurves.append(dict(heat_id=h,group=h[:2],wavelength_nm=x,median=np.interp(x,w,med),q10=np.interp(x,w,q10),q90=np.interp(x,w,q90)))
  # Actual-coordinate distances to reference wavelengths; assignment remains a candidate.
  for label,target in [('NaI_D2',588.995),('NaI_D1',589.5924),('KI_D2',766.48991),('KI_D1',769.89645),('O2_A_band_reference',760.0)]:
   j=int(np.argmin(abs(w-target)));a=max(0,j-6);b=min(len(w),j+7);peak=a+int(np.argmax(med[a:b]));phys.append(dict(heat_id=h,feature=label,reference_nm=target,reference_medium='air' if not label.startswith('O2') else 'band_neighborhood',nearest_channel=j,nearest_nm=w[j],grid_offset_nm=w[j]-target,local_median_peak_nm=w[peak],local_step_nm=float(np.median(np.diff(w[a:b]))),resolution_nm='configuration_unknown',identification='candidate_not_composition_confirmation'))
  # Recreate the detection working matrix for references, not the background solver.
  work=fill_training(np.where(~sat,y,np.nan),t,w);rt=work.copy();rw=work.copy();aa=(t[1:-1]-t[:-2])/(t[2:]-t[:-2]);bb=(w[1:-1]-w[:-2])/(w[2:]-w[:-2]);rt[1:-1]=(1-aa[:,None])*work[:-2]+aa[:,None]*work[2:];rw[:,1:-1]=(1-bb)*work[:,:-2]+bb*work[:,2:]
  near=maximum_filter(sat,size=(11,21),mode='nearest')
  with np.load(F/f'data/{h}_masks.npz') as mzdata:
   special=list(zip(*np.where(mask)))
   if h=='6C02108':special.append((608,1211))
   for i,j in special:
    low,high=neighbors(y,~mask,i,j)
    repair_rows.append(dict(heat_id=h,row=i,column=j,time_s=t[i],wavelength_nm=w[j],raw=y[i,j],repaired=D[i,j],difference=D[i,j]-y[i,j],latent=latent[i,j],time_reference=rt[i,j],wavelength_reference=rw[i,j],support_low=low,support_high=high,score=float(mzdata['score'][i,j]),time_z=float(mzdata['time_z'][i,j]),wave_z=float(mzdata['wave_z'][i,j]),lowrank_z=float(mzdata['lowrank_z'][i,j]),persistence=int(mzdata['persistence'][i,j]),fallback=bool(fall[i,j]),statistical_repair=bool(mask[i,j]),high_brightness=bool(bright[i]),rapid_frame=bool(fast[i]),near_top=bool(near[i,j]),high_raw=bool(y[i,j]>=.9*65535),top_code=bool(sat[i,j]),type=repairtype.loc[(h,i,j),'repair_type'] if mask[i,j] else 'retained_line'))
  # Correctly separate original sampled and column-centered decompositions.
  if h in PROTOCOL['svd']['heats']:
   Y=y[::8];ts=t[::8]
   for mode,mat in [('raw',Y),('centered',Y-Y.mean(axis=0))]:
    u,s,v=np.linalg.svd(mat,full_matrices=False);coef=u[:,:4]*s[:4]
    for k in range(4):
     sign=1 if v[k,np.argmax(abs(v[k]))]>=0 else -1;v[k]*=sign;coef[:,k]*=sign
    for r in PROTOCOL['svd']['ranks']:
     energy=float(np.sum(s[:r]**2)/np.sum(s**2));svrows.append(dict(heat_id=h,mode=mode,rank=r,energy=energy,unexplained=1-energy,rows=len(Y),time_stride=8))
    np.savez_compressed(S/f'data/{h}_svd_{mode}.npz',time=ts,wavelength=w,coefficients=coef,basis=v[:4],singular=s,background=coef@v[:4],rank1=np.outer(coef[:,0],v[0]),residual=mat-coef@v[:4],source=mode)
  summary.append(dict(heat_id=h,frames=len(t),observations=y.size,focus_nm=w[jfocus],focus_positive_fraction=pos[jfocus],focus_median_residual=md[jfocus],focus_z=mz[jfocus],persistent_candidates=int(candidate.sum()),baseline_repairs=int(mask.sum()),fallback=int((fall&mask).sum()),nonfallback=int((mask&~fall).sum())))
  print('analyzed',h,summary[-1],flush=True)
 for name,rows in [('channel_diagnostics',allchannels),('channel_followup',follow),('repair_reliability_points',repair_rows),('svd_conventions',svrows),('heat_summary',summary),('channel_sensitivity',counterfactual),('physical_feature_candidates',phys),('dark_bright_focus_spectra',targetsnap)]:pd.DataFrame(rows).to_csv(S/f'results/{name}.csv',index=False)
 diag=pd.DataFrame(allchannels);diag[diag.candidate].to_csv(S/'results/channel_candidates.csv',index=False)
 g=pd.DataFrame(groupcurves);g.to_csv(S/'results/common_grid_heat_summaries.csv',index=False);g.groupby(['group','wavelength_nm'])[['median','q10','q90']].mean().reset_index().to_csv(S/'results/common_grid_group_summaries.csv',index=False)
 rel=pd.DataFrame(repair_rows);a=rel[rel.statistical_repair];assert a.fallback.sum()==465 and (~a.fallback).sum()==578
 a.groupby('fallback')[['high_brightness','rapid_frame','near_top','high_raw']].mean().to_csv(S/'results/repair_reliability_groups.csv')
 a.groupby(['fallback','type']).size().rename('count').to_csv(S/'results/repair_reliability_types.csv')
 follow=pd.DataFrame(follow);follow[(follow.focus)&(follow.fraction==.15)].groupby('heat_id').first().to_csv(S/'results/focus_channels.csv')
 print('COMPLETE',flush=True)
if __name__=='__main__':main()
