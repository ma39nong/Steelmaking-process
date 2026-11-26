"""Descriptive raw-peak and persistent-ridge audit; raw readings are NOT clean ground truth."""
import numpy as np,pandas as pd
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *
from interpolation import load_surface

def main():
 rows=[];ridges=[]
 for h in HEATS:
  y,t,w=load(h);mask=np.load(P/f'data/{h}_masks.npz')['valid'];pred=load_surface(P/f'data/{h}_model.npz').grid;median=np.median(y,axis=0);peaks,_=find_peaks(median,prominence=max(50,float(median.max())*.01),distance=5)
  for j in peaks:
   good=y[:,j]<65535;ridges.append(dict(heat_id=h,wavelength_nm=w[j],median_intensity=median[j],masked_fraction=float(np.mean(~mask[good,j])),raw_fit_relative_difference=float(np.mean(abs(pred[good,j]-y[good,j])/np.maximum(y[good,j],1)))))
  raw=y.max(1);fitpeak=pred.max(1);q=raw<65535;relative=abs(fitpeak-raw)/np.maximum(raw,1);ii=np.flatnonzero(q);i=int(ii[np.argmax(relative[ii])]);rows.append(dict(heat_id=h,median_raw_peak_relative_change=float(np.median(relative[q])),p95_raw_peak_relative_change=float(np.quantile(relative[q],.95)),worst_unsaturated_raw_peak_relative_change=float(relative[i]),worst_time=t[i],not_clean_ground_truth=True))
  fig,ax=plt.subplots(figsize=(10,4));ax.plot(w,y[i],lw=.8,label='Raw (may contain anomalies)');ax.plot(w,pred[i],lw=.8,label='Locked V2 latent grid');ax.scatter(w[~mask[i]],y[i,~mask[i]],s=8,c='red',label='Masked');ax.set(title=f'{h} largest raw peak change at {t[i]:g}s',xlabel='Wavelength (nm)',ylabel='Intensity');ax.legend();fig.tight_layout();fig.savefig(P/f'figures/{h}_raw_peak_change.png',dpi=150);plt.close(fig)
 pd.DataFrame(rows).to_csv(P/'results/peak_retention_descriptive.csv',index=False);pd.DataFrame(ridges).to_csv(P/'results/persistent_ridge_retention.csv',index=False)
if __name__=='__main__':main()
