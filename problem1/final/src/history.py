"""Historical failure attribution, never used as the new confirmation benchmark."""
from core import *
from common31 import panels as old_panels,CFG as OLD
from scipy.ndimage import label,maximum_filter

def main():
 rows=[];counts=[]
 for h in HEATS:
  y,t,w=load(h);q=np.load(V1/f'data/masks/{h}.npz');types=q['types'];counts.append(dict(heat_id=h,V1_local_masked=int(((types&4)>0).sum()),V1_top_codes=int(((types&2)>0).sum()),V31_high_reviewed=int(np.load(V31/f'data/{h}_point_evidence.npz')['high'].sum())))
 pd.DataFrame(counts).to_csv(P/'results/history_counts.csv',index=False)
 for split in ['development','confirmation']:
  for h,y,t,w,i,j in old_panels(split):
   for kind in ['spike','drop']:
    for strength in [.3,.7]:
     rng=np.random.default_rng(OLD[split+'_seeds'][0]+sum(map(ord,h))+j);c,truth=inject(y,rng,kind,strength);e=conservative_evidence(c,t,w);floor=detect(c,t,w,window=5,threshold=6)[2]['noise_floor'];lab,n=label(e['candidate'],structure=np.ones((3,3)));sizes=np.bincount(lab.ravel());candidate=e['candidate'];high=e['high']
     measures=dict(time_z8=e['time_z']>=8,wave_z8=e['wave_z']>=8,lowrank_z6=e['residual_z']>=6,reference_agreement=e['neighbor_agreement'],outside_all_neighbors=e['exterior'],nonpersistent=e['persistence']<2,small_component=sizes[lab]<=3,final_high=high)
     base=dict(split=split,heat_id=h,column_start=j,kind=kind,strength=strength,truth_count=int(truth.sum()),wavelength_mid=float(np.mean(w)),reference_IQR=float(np.quantile(y,.95)-np.quantile(y,.05)),median_reference=float(np.median(y)),noise_floor=floor,median_truth_time_z=float(np.median(e['time_z'][truth])),median_truth_wave_z=float(np.median(e['wave_z'][truth])),median_truth_lowrank_z=float(np.median(e['residual_z'][truth])),median_injected_magnitude=float(np.median(abs(c[truth]-y[truth]))))
     for name,m in measures.items():base[name+'_pass_fraction']=float(m[truth].mean())
     rows.append(base)
 pd.DataFrame(rows).to_csv(P/'results/V31_failure_attribution.csv',index=False)
 pd.DataFrame(rows).groupby('split').mean(numeric_only=True).to_csv(P/'results/V31_failure_summary.csv');print('history done',flush=True)
if __name__=='__main__':main()
