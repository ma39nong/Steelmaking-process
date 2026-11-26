import sys,json,hashlib
from pathlib import Path
import numpy as np,pytest,pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from interface import Spectrum
from model3 import project_record

def test_bilinear_convex_hull_and_nodes():
 t=np.array([0.,1.,3.]);w=np.array([4.,5.,9.]);y=np.array([[0.,2.,3.],[1.,5.,8.],[7.,9.,10.]])
 s=Spectrum(t,w,y,y);assert np.allclose(s.evaluate(t,w),y)
 p=s.evaluate([.5],[4.5]);assert 0<=p[0,0]<=5
 with pytest.raises(ValueError):s.evaluate([-1],[4])
 with pytest.raises(ValueError):s.evaluate([1],[np.nan])
 assert not s.raw_until(.5)['offline_reconstruction_included']

def test_projection_retains_observed_and_bounds_hidden():
 y=np.array([[1.,2.,65535.],[3.,np.nan,6.],[4.,5.,7.]]);m=np.isfinite(y)&(y!=65535);L=np.full(y.shape,-100.)
 D,q,b=project_record(y,m,L,np.arange(3),np.arange(3),True,True)
 assert np.array_equal(D[m],y[m]);assert D[0,2]==65535;assert D[1,1]>=2;assert (D>=0).all();assert q[1,1]&16

def test_historical_and_numerical_immutability():
 for name,h in json.loads((P/'results/history_immutable.json').read_text()).items():assert hashlib.sha256((P.parents[1]/name).read_bytes()).hexdigest()==h,name
 for name,h in json.loads((P/'config/NUMERICAL_LOCKED.json').read_text()).items():assert hashlib.sha256((P/'src'/name).read_bytes()).hexdigest()==h,name
 lock=json.loads((P/'config/PROTOCOL_LOCKED.json').read_text());assert hashlib.sha256((P/'config/protocol.json').read_bytes()).hexdigest()==lock['protocol_sha256']

def test_all_final_artifacts():
 summary=pd.read_csv(P/'results/final_summary.csv');assert len(summary)==10
 for h in summary.heat_id:
  d=np.load(P/f'data/{h}_interface.npz');raw=np.load(P.parent/f'data/processed/{h}.npz')['Y'];assert np.array_equal(d['raw'],raw)
  assert np.isfinite(d['reconstructed']).all();assert (d['reconstructed']>=0).all();assert np.array_equal(d['reconstructed'][d['valid']],raw[d['valid']])
  sat=d['suspected_saturation'];assert (d['reconstructed'][sat]==65535).all();assert (d['conditional'][sat]>=65535).all();assert d['conditional_solver_converged'][sat].all()
  s=Spectrum.load(P/f'data/{h}_interface.npz');assert np.array_equal(s.evaluate(s.time[::90],s.wavelength[::90]),d['reconstructed'][::90,::90])
  dense=np.load(P/f'data/{h}_dense.npz')['recorded'];assert np.isfinite(dense).all() and (dense>=0).all()

def test_new_panel_disjoint_and_equal_latent():
 from ctx import panels,load
 for h,stage,y,t,w,i,j in panels('model'):
  yy,tt,ww=load(h);old=np.linspace(0,len(tt)-96,6).astype(int)
  assert all(i+32<=k or i>=k+96 for k in old)
 q=pd.read_csv(P/'results/projection_equivalence_confirmation.csv');assert q.select_dtypes('number').iloc[:,-1].max()==0
