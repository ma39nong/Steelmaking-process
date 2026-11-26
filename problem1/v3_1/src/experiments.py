"""End-to-end injected-contamination repair; raw reference never filtered by detector."""
from common31 import *
from points import evidence,mask_for,repair
from scipy.ndimage import maximum_filter1d

def stats(ref,corrupt,pred,truth,mask,w):
 error=pred-ref;before=corrupt-ref;untouched=~truth;tp=int((mask&truth).sum());fp=int((mask&~truth).sum());fn=int((~mask&truth).sum());tn=int((~mask&~truth).sum());precision=tp/max(tp+fp,1);recall=tp/max(tp+fn,1)
 area_ref=np.trapezoid(ref,w,axis=1);area_pred=np.trapezoid(pred,w,axis=1);peaks=ref.max(1);pr=pred.max(1)
 # Explicit significant unsupported NEW local maxima relative to the reference's +/-2 envelope.
 envelope=maximum_filter1d(ref,size=5,axis=1,mode='nearest');local=(pred>=np.roll(pred,1,1))&(pred>=np.roll(pred,-1,1));local[:,:2]=False;local[:,-2:]=False
 unsupported=local&(pred>envelope+np.maximum(500,.1*peaks[:,None]))
 corrupt_local=(corrupt>=np.roll(corrupt,1,1))&(corrupt>=np.roll(corrupt,-1,1));corrupt_unsupported=corrupt_local&(corrupt>envelope+np.maximum(500,.1*peaks[:,None]))
 return dict(introduced_unsupported_peak_count=int((unsupported&~corrupt_unsupported).sum()),tp=tp,fp=fp,fn=fn,tn=tn,precision=precision,recall=recall,f1=2*precision*recall/max(precision+recall,1e-12),false_positive_rate=fp/max(fp+tn,1),pollution_rmse=float(np.sqrt(np.mean(error[truth]**2))) if truth.any() else np.nan,before_pollution_rmse=float(np.sqrt(np.mean(before[truth]**2))) if truth.any() else np.nan,unpolluted_rmse=float(np.sqrt(np.mean(error[untouched]**2))),unpolluted_changed_fraction=float(np.mean(abs(error[untouched])>1e-8)),overall_rmse=float(np.sqrt(np.mean(error**2))),before_overall_rmse=float(np.sqrt(np.mean(before**2))),integral_relative_mae=float(np.mean(abs(area_pred-area_ref)/np.maximum(abs(area_ref),1))),peak_p95=float(np.quantile(abs(pr-peaks)/np.maximum(peaks,1),.95)),peak_position_mae_nm=float(np.mean(abs(w[pred.argmax(1)]-w[ref.argmax(1)]))),unsupported_peak_count=int(unsupported.sum()),negative_count=int((pred<0).sum()),masked_fraction=float(mask.mean()))

def inject(y,rng,kind,strength):
 # No detector used to choose reference positions: includes pre-existing unusual raw data.
 z=y.copy();truth=np.zeros(y.shape,bool);margin=4;eligible=np.zeros(y.shape,bool);eligible[margin:-margin,margin:-margin]=True;eligible&=(y>100)&(y<60000)
 positions=np.argwhere(eligible);positions=positions[rng.permutation(len(positions))];target=max(8,int(y.size*.005));selected=[]
 for i,j in positions:
  if any(abs(i-a)<=2 and abs(j-b)<=2 for a,b in selected):continue
  selected.append((int(i),int(j)))
  if len(selected)>=target:break
 amplitude=max(5000,float(np.quantile(y,.95)-np.quantile(y,.05)))*strength
 for i,j in selected:z[i,j]=min(64000,y[i,j]+amplitude) if kind=='spike' else y[i,j]*(1-strength)
 truth=z!=y
 return z,truth

def main(split):
 if split=='confirmation':assert (P/'config/SELECTION_LOCKED.json').exists()
 rows=[];controls=[]
 seeds=CFG[split+'_seeds']
 for h,y,t,w,i,j in panels(split):
  e=evidence(y,t,w)
  for name in CFG['strategies']:
   D,s,_=repair(y,t,w,e,name);mask=mask_for(e,name);out=stats(y,y,D,np.zeros(y.shape,bool),mask,w)
   # Actual high-prominence lines and real fastest transitions: retention diagnostics, not truth labels.
   from scipy.ndimage import median_filter
   prominent=(y-median_filter(y,size=(1,5),mode='nearest'))>np.maximum(500,.1*y.max(1)[:,None]);area=np.trapezoid(y,w,axis=1);rate=np.r_[0,abs(np.diff(area)/np.diff(t))];fast=rate>=np.quantile(rate,.9)
   controls.append(dict(heat_id=h,column_start=j,strategy=name,**out,real_prominent_cells=int(prominent.sum()),real_prominent_cells_masked=int((prominent&mask).sum()),fast_frame_cells_masked=int(mask[fast].sum()),**s))
  for seed in seeds:
   for kind in ['spike','drop']:
    for strength in [.3,.7]:
     rng=np.random.default_rng(seed+sum(map(ord,h))+j);corrupt,truth=inject(y,rng,kind,strength);e=evidence(corrupt,t,w)
     for name in CFG['strategies']:
      D,s,_=repair(corrupt,t,w,e,name);rows.append(dict(heat_id=h,column_start=j,seed=seed,kind=kind,strength=strength,strategy=name,**stats(y,corrupt,D,truth,mask_for(e,name),w),**s))
  pd.DataFrame(rows).to_csv(P/f'results/contamination_{split}.csv',index=False);pd.DataFrame(controls).to_csv(P/f'results/control_{split}.csv',index=False);print('experiment',split,h,j,flush=True)
if __name__=='__main__':main(sys.argv[1])
