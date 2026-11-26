import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline,make_interp_spline,RectBivariateSpline
class Surface:
 def __init__(self,t,w,grid,method='pchip'):
  self.time=np.asarray(t);self.wave=np.asarray(w);self.grid=grid;self.method=method
 def evaluate(self,t,w):
  t=np.atleast_1d(t);w=np.atleast_1d(w)
  if t.min()<self.time[0] or t.max()>self.time[-1] or w.min()<self.wave[0] or w.max()>self.wave[-1]:raise ValueError('Extrapolation prohibited')
  if self.method=='tensor_bspline':return RectBivariateSpline(self.time,self.wave,self.grid,kx=3,ky=3,s=0)(t,w)
  def interp(x,y,xq,axis):
   if self.method=='pchip':return PchipInterpolator(x,y,axis=axis,extrapolate=False)(xq)
   if self.method=='cubic':return CubicSpline(x,y,axis=axis,extrapolate=False)(xq)
   if self.method=='linear':return make_interp_spline(x,y,k=1,axis=axis)(xq)
   raise ValueError(self.method)
  z=interp(self.time,self.grid,t,0);return interp(self.wave,z,w,1)
 def save(self,path,**metadata):np.savez_compressed(path,time=self.time,wavelength=self.wave,grid=self.grid,interpolation=self.method,**metadata)
def load_surface(path):
 d=np.load(path);return Surface(d['time'],d['wavelength'],d['grid'],str(d['interpolation']))
