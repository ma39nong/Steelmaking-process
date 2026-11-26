from core import *
from detectors import all_masks,reconstruct

def metrics(ref,corrupt,D,truth,mask,w):
 s=stats(ref,corrupt,D,truth,mask,w);pr=ref.max(1);s['peak_mean_relative']=float(np.mean(abs(D.max(1)-pr)/np.maximum(pr,1)));base=max(s['before_overall_rmse']**2,1.)
 s['selection_loss']=s['overall_rmse']**2/base+.5*s['unpolluted_rmse']**2/base+.1*(1-s['f1'])+.25*s['peak_mean_relative']+.25*s['integral_relative_mae']+s['introduced_unsupported_peak_count']/max(int(truth.sum()),1)
 return s

def main(split):
 if split=='confirmation':assert (P/'config/SELECTION_LOCKED.json').exists()
 rows=[];controls=[]
 for h,band,y,t,w,i,j in panels(split):
  masks,e=all_masks(y,t,w,h)
  for name,mask in masks.items():
   D,sol=reconstruct(y,t,w,mask);out=metrics(y,y,D,np.zeros(y.shape,bool),mask,w);area=np.trapezoid(y,w,axis=1);rate=np.r_[0,abs(np.diff(area)/np.diff(t))];fast=rate>=np.quantile(rate,.9);controls.append(dict(heat_id=h,band=band,strategy=name,**out,fast_frame_cells_modified=int(((D!=y)&fast[:,None]).sum()),**sol))
  for seed in CFG[split+'_seeds']:
   for kind in ['spike','drop']:
    for strength in [.3,.7]:
     rng=np.random.default_rng(seed+sum(map(ord,h))+j);corrupt,truth=inject(y,rng,kind,strength);masks,e=all_masks(corrupt,t,w,h)
     for name,mask in masks.items():
      D,sol=reconstruct(corrupt,t,w,mask);rows.append(dict(heat_id=h,band=band,seed=seed,kind=kind,strength=strength,strategy=name,**metrics(y,corrupt,D,truth,mask,w),**sol))
  pd.DataFrame(rows).to_csv(P/f'results/contamination_{split}.csv',index=False);pd.DataFrame(controls).to_csv(P/f'results/control_{split}.csv',index=False);print(split,h,band,flush=True)
 if split=='development':
  d=pd.DataFrame(rows);by=d.groupby(['heat_id','strategy']).mean(numeric_only=True).reset_index();r=by.groupby('strategy')[['selection_loss','unpolluted_rmse','pollution_rmse','f1','peak_mean_relative','integral_relative_mae']].mean().reset_index();r.to_csv(P/'results/development_selection.csv',index=False);best=r.selection_loss.min();pool=r[r.selection_loss<=best+.01];chosen=pool.sort_values(['unpolluted_rmse','selection_loss']).iloc[0].strategy
  dump(P/'config/SELECTION_LOCKED.json',dict(strategy=chosen,confirmation_unopened=True,protocol_sha256=sha(P/'config/protocol.json'),detector_sha256=sha(P/'src/detectors.py'),experiment_sha256=sha(P/'src/experiment_final.py'),rule='predeclared continuous benefit-cost loss, no absolute gate',selection_table=r.to_dict('records')))
if __name__=='__main__':main(sys.argv[1])
