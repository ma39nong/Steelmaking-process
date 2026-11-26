"""Descriptive peak/ridge retention checks; raw observations are not asserted to be clean truth."""
import numpy as np,pandas as pd
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from audit import P,HEATS
from models import load_surface
from eda import save
rows=[];ridges=[]
for h in HEATS:
 d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];mask=np.load(P/f'data/masks/{h}.npz')['valid'];model=load_surface(P/f'data/reconstructed/{h}_model.npz');pred=model.grid
 median=np.median(y,axis=0);peaks,_=find_peaks(median,prominence=max(50,float(median.max())*.01),distance=5)
 for j in peaks:
  valid=(y[:,j]<65535)
  ridges.append(dict(heat_id=h,wavelength_nm=w[j],median_intensity=median[j],masked_fraction=float(np.mean(~mask[valid,j])),raw_fit_relative_difference=float(np.mean(abs(pred[valid,j]-y[valid,j])/np.maximum(y[valid,j],1)))))
 rawpeak=np.max(y,1);fitpeak=np.max(pred,1);eligible=rawpeak<65535;relative=abs(fitpeak-rawpeak)/np.maximum(rawpeak,1);candidates=np.flatnonzero(eligible);i=int(candidates[np.argmax(relative[candidates])]);j=np.argmax(y,1)
 rows.append(dict(heat_id=h,raw_peak_masked_fraction=float(np.mean(~mask[np.arange(len(t)),j])),raw_fit_peak_difference_median=float(np.median(relative[eligible])),raw_fit_peak_difference_p95=float(np.quantile(relative[eligible],.95)),worst_unsaturated_raw_peak_difference=float(relative[i]),worst_time_s=t[i],interpretation='descriptive difference from raw; not true-clean error'))
 fig,ax=plt.subplots(figsize=(10,4));ax.plot(w,y[i],lw=.9,label='Raw (may include anomalies)');ax.plot(w,pred[i],lw=.9,label='Final continuous fit');ax.scatter(w[~mask[i]],y[i,~mask[i]],color='red',s=10,label='Masked');ax.set(xlabel='Wavelength (nm)',ylabel='Intensity',title=f'{h}: largest raw peak difference, t={t[i]:g}s');ax.legend(fontsize=8);save(fig,P/f'figures/comparison/{h}_largest_raw_peak_difference')
pd.DataFrame(rows).to_csv(P/'results/validation/peak_retention_descriptive.csv',index=False);pd.DataFrame(ridges).to_csv(P/'results/validation/persistent_ridges.csv',index=False)
