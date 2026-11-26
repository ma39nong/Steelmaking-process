from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit import P,HEATS

def save(fig,path):
 fig.tight_layout();fig.savefig(str(path)+'.png',dpi=150);fig.savefig(str(path)+'.svg');plt.close(fig)
def heatmap(y,t,w,path,title,vmax=None):
 fig,ax=plt.subplots(figsize=(10,4));m=ax.pcolormesh(w,t,y,shading='auto',cmap='magma',vmin=0,vmax=vmax or np.nanmax(y),rasterized=True);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)',title=title);fig.colorbar(m,ax=ax,label='Intensity (recorded units)');save(fig,path)
def run():
 if not (P/'results/audit/COMPLETE.json').exists():raise RuntimeError('Audit required')
 out=[]
 for h in HEATS:
  d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];print(h,flush=True)
  heatmap(y,t,w,P/f'figures/raw/{h}_heatmap',h+' raw')
  fig,axs=plt.subplots(2,2,figsize=(12,7))
  for i in [len(t)//10,len(t)//2,len(t)*9//10]:axs[0,0].plot(w,y[i],label=f'{t[i]:g} s',lw=.8)
  axs[0,0].legend();axs[0,0].set(xlabel='Wavelength (nm)',ylabel='Intensity')
  for v in [400,589,770,950]:
   j=np.argmin(abs(w-v));axs[0,1].plot(t,y[:,j],lw=.7,label=f'{w[j]:.2f} nm')
  axs[0,1].legend();axs[0,1].set(xlabel='Time (s)',ylabel='Intensity')
  axs[1,0].hist(y[::5,::5].ravel(),bins=100);axs[1,0].set(xlabel='Intensity',ylabel='Count (sampled 1/25)')
  axs[1,1].plot(t,w[y.argmax(1)],lw=.7);axs[1,1].set(xlabel='Time (s)',ylabel='Maximum wavelength (nm)')
  fig.suptitle(h);save(fig,P/f'figures/raw/{h}_eda')
  # Correlations calculated on actual adjacent pairs, independently by heat.
  def corr(a,b):return float(np.corrcoef(a.ravel(),b.ravel())[0,1])
  out.append(dict(heat_id=h,adjacent_time_correlation=corr(y[:-1,::4],y[1:,::4]),adjacent_wavelength_correlation=corr(y[::4,:-1],y[::4,1:]),median_abs_time_derivative=float(np.median(abs(np.diff(y,axis=0)/np.diff(t)[:,None]))),median_abs_wavelength_derivative=float(np.median(abs(np.diff(y,axis=1)/np.diff(w)[None,:]))),p99=float(np.quantile(y,.99)),median=float(np.median(y)),q95_peak_nm=float(np.quantile(w[y.argmax(1)],.95))))
  sat=y==65535
  if sat.any():
   i,j=np.argwhere(sat)[len(np.argwhere(sat))//2];a,b=max(0,i-12),min(len(t),i+13);c,e=max(0,j-30),min(len(w),j+31)
   heatmap(y[a:b,c:e],t[a:b],w[c:e],P/f'figures/raw/{h}_saturation_detail',h+' possible clipping')
 pd.DataFrame(out).to_csv(P/'results/eda/statistics.csv',index=False)
 (P/'results/eda/COMPLETE').write_text('Ten heats audited and EDA plots generated.\n')
if __name__=='__main__':run()
