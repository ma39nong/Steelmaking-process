"""Three explicit model families; no gigantic flattened tensor design matrix."""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import BSpline,RectBivariateSpline,PchipInterpolator
from scipy.sparse import csr_matrix,diags
from scipy.sparse.linalg import LinearOperator,cg
from scipy.special import huber
from scipy.signal import savgol_filter
from sklearn.utils.extmath import randomized_svd
from anomalies import fill_training
SCALE=65535.

def local_poly_matrix(w,window=7,degree=2):
 n=len(w);rows=[];cols=[];vals=[]
 for j in range(n):
  start=max(0,min(j-window//2,n-window));ix=np.arange(start,start+window);x=w[ix]-w[j];x/=max(abs(x).max(),1e-12);V=np.vander(x,degree+1,increasing=True)
  weights=np.linalg.pinv(V)[0];rows.extend([j]*len(ix));cols.extend(ix);vals.extend(weights)
 return csr_matrix((vals,(rows,cols)),shape=(n,n))

@dataclass
class GridSurface:
 time:np.ndarray
 wave:np.ndarray
 grid:np.ndarray
 family:str
 params:dict
 diagnostics:dict
 latent_grid:object=None
 def evaluate(self,t,w):
  t=np.atleast_1d(t);w=np.atleast_1d(w)
  if min(t)<self.time[0] or max(t)>self.time[-1] or min(w)<self.wave[0] or max(w)>self.wave[-1]:raise ValueError('Extrapolation is prohibited')
  if self.family=='B' and self.params.get('interpolation','pchip')=='pchip':
   temporal=PchipInterpolator(self.time,self.grid,axis=0,extrapolate=False)(t)
   return PchipInterpolator(self.wave,temporal,axis=1,extrapolate=False)(w)
  return RectBivariateSpline(self.time,self.wave,self.grid,kx=3,ky=3,s=0)(t,w)
 def save(self,path):
  np.savez_compressed(path,kind='grid',time=self.time,wavelength=self.wave,grid=self.grid,family=self.family,params_json=__import__('json').dumps(self.params),diagnostics_json=__import__('json').dumps(self.diagnostics),latent_grid=self.latent_grid if self.latent_grid is not None else np.empty((0,0)))

@dataclass
class TensorSurface:
 time:np.ndarray
 wave:np.ndarray
 knots_t:np.ndarray
 knots_w:np.ndarray
 coefficients:np.ndarray
 params:dict
 diagnostics:dict
 family:str='C'
 def evaluate(self,t,w):
  t=np.atleast_1d(t);w=np.atleast_1d(w)
  if min(t)<self.time[0] or max(t)>self.time[-1] or min(w)<self.wave[0] or max(w)>self.wave[-1]:raise ValueError('Extrapolation is prohibited')
  bt=BSpline.design_matrix((t-self.time[0])/np.ptp(self.time),self.knots_t,3,extrapolate=False)
  bw=BSpline.design_matrix((w-self.wave[0])/np.ptp(self.wave),self.knots_w,3,extrapolate=False)
  return (bw@(bt@self.coefficients).T).T*SCALE
 @property
 def grid(self):return self.evaluate(self.time,self.wave)
 def save(self,path):
  np.savez_compressed(path,kind='tensor',time=self.time,wavelength=self.wave,knots_t=self.knots_t,knots_w=self.knots_w,coefficients=self.coefficients,family=self.family,params_json=__import__('json').dumps(self.params),diagnostics_json=__import__('json').dumps(self.diagnostics))

def load_surface(path):
 import json
 d=np.load(path);p=json.loads(str(d['params_json']));diag=json.loads(str(d['diagnostics_json']))
 if str(d['kind'])=='tensor':return TensorSurface(d['time'],d['wavelength'],d['knots_t'],d['knots_w'],d['coefficients'],p,diag)
 return GridSurface(d['time'],d['wavelength'],d['grid'],str(d['family']),p,diag,d['latent_grid'] if 'latent_grid' in d and d['latent_grid'].size else None)

def fit_a(y,t,w,mask,params):
 z=fill_training(np.where(mask,y,np.nan),t,w)
 tw=params['time_window'];ww=params['wave_window']
 # Audited grids are time-uniform. A general missing Sequence must not silently use row SG.
 if not np.allclose(np.diff(t),np.diff(t)[0]):raise ValueError('SG time axis is not uniform')
 z=savgol_filter(z,tw,2,axis=0,mode='interp')
 z=(local_poly_matrix(w,ww)@z.T).T
 return GridSurface(t,w,z,'A',params,{'solver':'actual-coordinate local polynomial + tensor cubic interpolation','converged':True})

def fit_b(y,t,w,mask,params):
 rank=params['rank'];eta=params['eta'];z=np.where(mask,y/SCALE,0);L=fill_training(np.where(mask,y/SCALE,np.nan),t,w);S=np.zeros_like(L);hist=[]
 for k in range(params.get('max_iter',18)):
  S=np.where(mask,np.sign(z-L)*np.maximum(abs(z-L)-eta,0),0)
  a=np.where(mask,z-S,L)
  u,s,v=randomized_svd(a,n_components=rank,n_iter=3,random_state=20261126)
  new=(u*s)@v;change=np.linalg.norm(new-L)/max(np.linalg.norm(L),1e-12);L=new
  residual=np.where(mask,z-L-S,0);obj=.5*np.sum(residual**2)+eta*np.sum(abs(S));hist.append({'objective':float(obj),'relative_change':float(change)})
  if change<params.get('tol',1e-4):break
 latent=L*SCALE
 grid=np.where(mask,y,latent) if params.get('preserve_observed',True) else latent
 return GridSurface(t,w,grid,'B',params,{'solver':'masked projected alternating rank/S shrinkage; randomized SVD approximation','converged':bool(change<params.get('tol',1e-4)),'history':hist,'rank':rank,'observation_projection':params.get('preserve_observed',True)},latent)

def basis(x,k):
 knots=np.r_[np.zeros(4),np.linspace(0,1,k-2)[1:-1],np.ones(4)]
 return BSpline.design_matrix((x-x[0])/np.ptp(x),knots,3,extrapolate=False).tocsr(),knots

def difference(k,order):
 import math
 return diags([np.ones(k-order)*((-1)**(order-j))*math.comb(order,j) for j in range(order+1)],range(order+1),shape=(k-order,k),format='csr')

def fit_c(y,t,w,mask,params):
 kt=params['kt'];kw=params['kw'];bt,nt=basis(t,kt);bw,nw=basis(w,kw);n=kt*kw
 dt=difference(kt,2);dw=difference(kw,2);pt=dt.T@dt;pw=dw.T@dw;g1=difference(kt,1);g2=difference(kw,1);gt=g1.T@g1;gw=g2.T@g2
 z=np.where(mask,y/SCALE,0);bt2=bt.multiply(bt);bw2=bw.multiply(bw)
 cell=mask.size/n;at=params['alpha_t']*cell;aw=params['alpha_w']*cell;ax=params.get('alpha_cross',0)*cell;ridge=1e-8
 def forward(c):return (bw@(bt@c).T).T
 def adj(v):return (bw.T@(bt.T@v).T).T
 def pen(c):return at*(pt@c)+aw*(pw@c.T).T+ax*(gw@(gt@c).T).T+ridge*c
 # Noise scale from only visible wavelength residuals, independent of held-out values.
 visible=np.where(mask,y,np.nan);filled=fill_training(visible,t,w)
 r=filled[:,1:-1]-(filled[:,:-2]+filled[:,2:])*.5
 scale=max(25.,float(1.4826*np.median(abs(r-np.median(r)))))/SCALE
 delta=params['delta'];level=params.get('censor_level',65535.)/SCALE;sat=np.isfinite(y)&(y==params.get('censor_level',65535.));eta=params.get('censor_weight',0)
 c=np.zeros((kt,kw));history=[];W=mask.astype(float)
 max_irls=params.get('irls',4)
 for it in range(max_irls):
  pred=forward(c)
  if it:
   r=(z-pred)/scale;W=mask*np.minimum(1,delta/np.maximum(abs(r),1e-12))
  active=sat&(pred<level) if eta else np.zeros(mask.shape,bool)
  weights=W+eta*active;target=W*z+eta*active*level
  diagonal=adj_diag=(bw2.T@(bt2.T@weights).T).T
  diagonal=diagonal+at*pt.diagonal()[:,None]+aw*pw.diagonal()[None,:]+ax*gt.diagonal()[:,None]*gw.diagonal()[None,:]+ridge
  def op(v):
   a=v.reshape(kt,kw);return (adj(weights*forward(a))+pen(a)).ravel()
  A=LinearOperator((n,n),matvec=op,dtype=float);pre=LinearOperator((n,n),matvec=lambda v:v/diagonal.ravel(),dtype=float)
  count=[0]
  def cb(v):count[0]+=1
  rhs=adj(target).ravel();answer,info=cg(A,rhs,x0=c.ravel(),M=pre,rtol=params.get('cg_tol',2e-4),maxiter=params.get('cg_iter',100),callback=cb)
  new=answer.reshape(kt,kw);pred=forward(new);r=(z-pred)/scale
  objective=float(scale**2*np.sum(mask*huber(delta,r))+.5*np.sum(new*pen(new))+.5*eta*np.sum(np.maximum(level-pred[sat],0)**2))
  change=np.linalg.norm(new-c)/max(np.linalg.norm(new),1e-12);c=new
  history.append(dict(irls=it,cg_info=int(info),cg_iterations=count[0],linear_relative_residual=float(np.linalg.norm(A@answer-rhs)/max(np.linalg.norm(rhs),1e-12)),objective=objective,relative_change=float(change)))
  if change<params.get('irls_tol',1e-3):break
  if len(history)>1 and abs(history[-1]['objective']-history[-2]['objective'])/max(abs(history[-2]['objective']),1e-12)<params.get('objective_tol',1e-4):break
 return TensorSurface(t,w,nt,nw,c,params,{'solver':'matrix-free IRLS/PCG','scale_intensity':scale*SCALE,'history':history,'converged':bool(info==0 and (change<params.get('irls_tol',1e-3) or (len(history)>1 and abs(history[-1]['objective']-history[-2]['objective'])/max(abs(history[-2]['objective']),1e-12)<params.get('objective_tol',1e-4)))),'cg_converged':bool(info==0),'penalty_scale':cell})

def fit(family,y,t,w,mask,params):return {'A':fit_a,'B':fit_b,'C':fit_c}[family](y,t,w,mask,params)
