"""Run-length and duplicate-channel audit from cached exact raw matrices."""
import numpy as np,pandas as pd
from audit import HEATS,P
rows=[];intervals=[]
for h in HEATS:
 d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];run=np.ones(y.shape[1],int);maximum=np.ones_like(run);end=np.zeros_like(run)
 for i in range(1,len(t)):
  run=np.where(y[i]==y[i-1],run+1,1);improved=run>maximum;maximum=np.maximum(maximum,run);end[improved]=i
 unique,first,inverse,counts=np.unique(y.T,axis=0,return_index=True,return_inverse=True,return_counts=True)
 duplicate_columns=int(len(w)-len(unique))
 for j in np.flatnonzero(maximum>=10):intervals.append(dict(heat_id=h,wavelength_nm=w[j],constant_samples=int(maximum[j]),start_s=t[end[j]-maximum[j]+1],end_s=t[end[j]],value=y[end[j],j]))
 # Very large one-frame broadband jumps are diagnostics, not automatic exclusion.
 mean=y.mean(1);dr=np.diff(mean);scale=max(1,1.4826*np.median(abs(dr-np.median(dr))));jumps=abs(dr-np.median(dr))>8*scale
 rows.append(dict(heat_id=h,max_constant_run_samples=int(maximum.max()),max_constant_run_s=float((maximum.max()-1)*.5),channels_with_at_least_10_constant_samples=int(np.sum(maximum>=10)),duplicate_channels=duplicate_columns,zero_runs=0 if not np.any(y==0) else np.nan,broadband_jump_candidates=int(jumps.sum())))
 pd.DataFrame(dict(wavelength_nm=w,max_constant_run_samples=maximum,duplicate_group=inverse,duplicate_group_size=counts[inverse])).to_csv(P/f'results/audit/{h}_channel_runs.csv',index=False)
pd.DataFrame(rows).to_csv(P/'results/audit/extended_quality.csv',index=False);pd.DataFrame(intervals).to_csv(P/'results/audit/constant_intervals.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
