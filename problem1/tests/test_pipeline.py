from pathlib import Path
import sys,json,hashlib
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'problem1/src'))
from audit import HEATS,NEW_HASH,P,source_files,read_heat
from anomalies import fill_training,detect
from models import basis,difference,fit,load_surface,local_poly_matrix

@pytest.mark.parametrize('h',HEATS)
def test_ten_heat_read(h):
 file=next(p for p in source_files() if h in p.name)
 sheet,header,sequence,t,w,y,extra=read_heat(file)
 assert sheet==h and y.shape==(len(t),2047) and len(header)==2051
 assert extra.shape==(len(t),3)
 assert np.allclose(np.diff(t),.5) and np.all(np.diff(w)>0)
 assert not np.allclose(np.diff(w),np.diff(w)[0])
 if h==HEATS[0]:
  assert hashlib.sha256(file.read_bytes()).hexdigest()==NEW_HASH
  assert y.shape[0]==1479 and w[-1]==1027.963
  assert np.count_nonzero(extra[:,1])>0

def test_version_replacement():
 files=source_files();assert len(files)==10
 assert hashlib.sha256(files[0].read_bytes()).hexdigest()==NEW_HASH
 assert len([p for p in files if HEATS[0] in p.name])==1

def test_nonuniform_local_polynomial():
 w=np.array([0,.2,.7,1.2,1.8,2.1,2.8,3.2,4.0]);y=2+3*w-2*w*w
 np.testing.assert_allclose(local_poly_matrix(w,5)@y,y,atol=1e-10)

def test_mask_and_real_peak_persistence():
 t=np.arange(31)*.5;w=np.linspace(400,700,41);y=np.ones((31,41))*1000;y[:,20]=3000;y[12,9]=50000;y[20:23,30:33]=65535
 valid,types,_=detect(y,t,w,window=5,threshold=6)
 assert not valid[12,9] and np.all((types[20:23,30:33]&2)>0)
 assert np.all(valid[:,20]) # Narrow persistent line is kept.

def test_visibility_and_full_missing_rows_columns():
 t=np.arange(12)*.5;w=np.array([1,2,4,7,10,14.]);a=t[:,None]+w[None,:];x=a.copy();x[3:6,:]=np.nan;x[:,2]=np.nan
 out=fill_training(x,t,w);assert np.isfinite(out).all();np.testing.assert_allclose(out[np.isfinite(x)],x[np.isfinite(x)])

def test_tensor_normal_equation_against_dense():
 t=np.linspace(0,1,20);w=np.linspace(2,4,25)**1.1;y=65535*(.1+.03*np.sin(3*t[:,None])+0*w);mask=np.ones(y.shape,bool);mask[4:6,7:9]=False
 p=dict(kt=7,kw=9,alpha_t=.03,alpha_w=.02,delta=1e8,irls=1,cg_tol=1e-11,cg_iter=500)
 m=fit('C',np.where(mask,y,np.nan),t,w,mask,p);bt,_=basis(t,7);bw,_=basis(w,9);A=np.kron(bt.toarray(),bw.toarray());dt=difference(7,2).toarray();dw=difference(9,2).toarray();cell=mask.size/(7*9)
 penalty=.03*cell*np.kron(dt.T@dt,np.eye(9))+.02*cell*np.kron(np.eye(7),dw.T@dw)+1e-8*np.eye(63)
 good=mask.ravel();expected=np.linalg.solve(A[good].T@A[good]+penalty,A[good].T@(y.ravel()[good]/65535))
 np.testing.assert_allclose(m.coefficients.ravel(),expected,rtol=1e-7,atol=1e-8)
 assert np.isfinite(m.evaluate([.2,.21],[w[2],w[2]+.01])).all()
 with pytest.raises(ValueError):m.evaluate([-1],w)

def test_model_parameter_roundtrip(tmp_path):
 t=np.arange(12)*.5;w=np.linspace(400,600,16);y=1000+5*t[:,None]+np.sin(w[None,:]);p=dict(time_window=3,wave_window=5)
 m=fit('A',y,t,w,np.ones(y.shape,bool),p);path=tmp_path/'model.npz';m.save(path);other=load_surface(path)
 np.testing.assert_allclose(m.evaluate(t,w),other.evaluate(t,w));assert other.params==p

def test_results_and_ten_plots():
 if not (P/'results/final/COMPLETE').exists():pytest.skip('Final run not yet completed; rerun after final.py')
 for h in HEATS:
  assert (P/f'figures/reconstructed/{h}_surface.png').exists()
  d=np.load(P/f'data/reconstructed/{h}_dense.npz');assert d['intensity'].shape==(len(d['time']),len(d['wavelength'])) and np.isfinite(d['intensity']).all()
  m=load_surface(P/f'data/reconstructed/{h}_model.npz');assert np.isfinite(m.evaluate(m.time[::50],m.wave[::80])).all()
  if m.family=='B':
   original=np.load(P/f'data/processed/{h}.npz');valid=np.load(P/f'data/masks/{h}.npz')['valid'];np.testing.assert_array_equal(m.grid[valid],original['Y'][valid]);assert m.params['interpolation']=='pchip'
 import pandas as pd
 r=pd.read_csv(P/'results/model_comparison/model_comparison.csv');assert len(r)==150 and r.heat_id.nunique()==10 and set(r.model)=={'A','B','C'}

def test_low_rank_observation_projection(tmp_path):
 t=np.arange(14)*.5;w=np.linspace(400,700,22);y=1000+30*np.sin(t[:,None])+10*np.cos(w[None,:]);y[5,7]=3500
 mask=np.ones(y.shape,bool);mask[3:5,10:12]=False;visible=np.where(mask,y,np.nan)
 model=fit('B',visible,t,w,mask,dict(rank=2,eta=.003,max_iter=20,preserve_observed=True))
 np.testing.assert_array_equal(model.grid[mask],y[mask])
 np.testing.assert_array_equal(model.grid[~mask],model.latent_grid[~mask])
 path=tmp_path/'projected.npz';model.save(path);loaded=load_surface(path)
 np.testing.assert_array_equal(loaded.latent_grid,model.latent_grid)
 np.testing.assert_allclose(loaded.evaluate(t,w)[mask],y[mask],atol=1e-8)

def test_shape_preserving_continuous_evaluation():
 from models import GridSurface
 t=np.arange(7,dtype=float);w=np.array([400,401,404,410,420.]);grid=np.array([1,1,100,1,90,1,1.])[:,None]*np.array([1,2,1,4,2.])[None,:]
 model=GridSurface(t,w,grid,'B',{'interpolation':'pchip'}, {'converged':True})
 dense=model.evaluate(np.linspace(t[0],t[-1],121),np.linspace(w[0],w[-1],101))
 assert np.min(dense)>=np.min(grid)-1e-10 and np.max(dense)<=np.max(grid)+1e-10
 np.testing.assert_allclose(model.evaluate(t,w),grid,atol=1e-12)
