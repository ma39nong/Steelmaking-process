"""Stable offline recorded-spectrum interface. Raw causal slicing never returns offline imputations."""
from pathlib import Path
import json
import numpy as np
class Spectrum:
 def __init__(self,time,wavelength,grid,raw=None,quality=None,conditional=None):
  self.time=np.asarray(time);self.wavelength=np.asarray(wavelength);self.grid=np.asarray(grid);self.raw=raw;self.quality=quality;self.conditional=conditional
  if self.grid.shape!=(len(self.time),len(self.wavelength)):raise ValueError('Axis/grid mismatch')
  if np.any(np.diff(self.time)<=0) or np.any(np.diff(self.wavelength)<=0):raise ValueError('Coordinates must increase')
  if not np.isfinite(self.grid).all() or (self.grid<0).any():raise ValueError('Official recorded grid must be finite and nonnegative')
 def evaluate(self,t,w,target='recorded'):
  t=np.atleast_1d(np.asarray(t,float));w=np.atleast_1d(np.asarray(w,float))
  if not (np.isfinite(t).all() and np.isfinite(w).all()):raise ValueError('Nonfinite coordinate')
  if t.min()<self.time[0] or t.max()>self.time[-1] or w.min()<self.wavelength[0] or w.max()>self.wavelength[-1]:raise ValueError('Extrapolation prohibited')
  grid=self.grid if target=='recorded' else self.conditional if target=='conditional' else None
  if grid is None:raise ValueError('Unknown/unavailable target')
  i=np.clip(np.searchsorted(self.time,t,side='right')-1,0,len(self.time)-2);j=np.clip(np.searchsorted(self.wavelength,w,side='right')-1,0,len(self.wavelength)-2);a=(t-self.time[i])/(self.time[i+1]-self.time[i]);b=(w-self.wavelength[j])/(self.wavelength[j+1]-self.wavelength[j]);a=a[:,None];b=b[None,:]
  return (1-a)*(1-b)*grid[np.ix_(i,j)]+a*(1-b)*grid[np.ix_(i+1,j)]+(1-a)*b*grid[np.ix_(i,j+1)]+a*b*grid[np.ix_(i+1,j+1)]
 def raw_until(self,cutoff):
  if self.raw is None:raise ValueError('Raw data unavailable')
  keep=self.time<=cutoff
  return dict(time=self.time[keep].copy(),wavelength=self.wavelength.copy(),raw=self.raw[keep].copy(),causal_preprocessing_performed=False,offline_reconstruction_included=False)
 @classmethod
 def load(cls,path):
  d=np.load(path);return cls(d['time'],d['wavelength'],d['reconstructed'],d['raw'],d['quality'],d['conditional'])

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('heat');p.add_argument('--time',nargs='+',type=float,required=True);p.add_argument('--wave',nargs='+',type=float,required=True);p.add_argument('--target',default='recorded');a=p.parse_args();m=Spectrum.load(Path(__file__).resolve().parents[1]/f'data/{a.heat}_interface.npz');print(m.evaluate(a.time,a.wave,a.target))
