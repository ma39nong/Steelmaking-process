import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from solver import fit,loss_grad,project,time_weights
from interpolation import Surface
from common import panels,CFG,P,digest
from convex_censor import fit_convex

def test_gradient():
 r=np.random.default_rng(2);L=r.normal(size=(9,11));z=r.normal(size=L.shape);mask=r.random(L.shape)>.2;sat=~mask;weights=r.random(7)
 val,g=loss_grad(L,z,mask,.03,sat,1.3,.2,weights);v=r.normal(size=L.shape);e=1e-6
 numerical=(loss_grad(L+e*v,z,mask,.03,sat,1.3,.2,weights)[0]-loss_grad(L-e*v,z,mask,.03,sat,1.3,.2,weights)[0])/(2*e)
 assert np.isclose(numerical,(g*v).sum(),rtol=1e-6,atol=1e-6)
def test_sparse_elimination():
 x=np.linspace(-2,2,101);eta=.2;s=np.sign(x)*np.maximum(abs(x)-eta,0);direct=.5*(x-s)**2+eta*abs(s);hub=np.where(abs(x)<=eta,.5*x*x,eta*abs(x)-.5*eta*eta);assert np.allclose(direct,hub)
def test_b2_b3_hidden_prediction():
 t=np.arange(12.);w=np.arange(15.);y=np.outer(1+t,2+w);mask=np.ones(y.shape,bool);mask[3:5,5:9]=False;v=np.where(mask,y,np.nan)
 a=fit(v,t,w,mask,'B2',rank=2);b=fit(v,t,w,mask,'B3',rank=2)
 assert np.array_equal(a['grid'][~mask],b['grid'][~mask]);assert np.array_equal(b['grid'][mask],y[mask]);assert all(x['objective']>=z['objective']-1e-10 for x,z in zip(a['history'],a['history'][1:]))
def test_hidden_values_cannot_affect_fit():
 t=np.arange(12.);w=np.arange(15.);y=np.outer(1+t,2+w);mask=np.ones(y.shape,bool);mask[4:6]=False
 a=y.copy();a[~mask]=np.nan;b=y.copy();b[~mask]=np.nan
 assert np.array_equal(fit(a,t,w,mask,'B5',rank=2)['grid'],fit(b,t,w,mask,'B5',rank=2)['grid'])
def test_convex_true_bound():
 t=np.arange(12.);w=np.arange(15.);y=np.outer(1+t,2+w);cap=100.;sat=y>cap;visible=np.minimum(y,cap);mask=~sat
 r=fit_convex(visible,t,w,mask,sat,cap,max_iter=100);assert (r['grid'][sat]>=cap).all()
@pytest.mark.parametrize('method',['linear','pchip','cubic','tensor_bspline'])
def test_interpolation_affine(method):
 t=np.array([0,.7,1.2,2,3.]);w=np.array([300,302,307,308,315.]);z=3*t[:,None]+w;model=Surface(t,w,z,method);tt=np.linspace(0,3,14);ww=np.linspace(300,315,17);assert np.allclose(model.evaluate(tt,ww),3*tt[:,None]+ww)
 with pytest.raises(ValueError):model.evaluate([-1],w)
def test_cubic_representation_equivalence():
 rng=np.random.default_rng(4);t=np.arange(7.);w=np.arange(8.)**1.1;z=rng.random((7,8));a=Surface(t,w,z,'cubic');b=Surface(t,w,z,'tensor_bspline');assert np.allclose(a.evaluate(t+.0,w),b.evaluate(t,w),atol=1e-10)
def test_panel_independence():
 dev={(h,s):(t,w) for h,s,y,t,w in panels('development')};test=list(panels('test'))
 for h,s,y,t,w in test:
  for ds in range(3):assert not np.intersect1d(t,dev[h,ds][0]).size

def test_v1_immutable():
 import json
 for f,sha in json.loads((P/'results/v1_immutable_manifest.json').read_text()).items():assert digest(P.parents[1]/f)==sha,f
