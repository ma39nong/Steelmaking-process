import sys,json
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ctx import P,CFG
from frame_review import frame_metrics
import preprocess
from preprocess import apply
from model3 import solve

def toy():
 t=np.arange(32)*.5;w=np.linspace(350,1000,64);y=np.outer(1000+np.arange(32)*2,1+np.sin(np.linspace(0,2,64))*.2);return y,t,w

def test_reference_excludes_current():
 y,t,w=toy();a=frame_metrics(y,t,w);z=y.copy();z[15]*=1.8;b=frame_metrics(z,t,w);assert a.reference_integral.iloc[15]==b.reference_integral.iloc[15]
def test_suspect_is_not_confirmed_fault():
 y,t,w=toy();y[15]*=1.7;p=apply(y,t,w,'P0');assert p['frames'].candidate.iloc[15];assert p['valid'][15].all();assert not p['invalid'].any();assert not p['frames'].confirmed_acquisition_fault.any()
def test_extent_confirmation_ledger(tmp_path,monkeypatch):
 (tmp_path/'config').mkdir();(tmp_path/'config/confirmed_invalid.json').write_text(json.dumps(dict(frames={'H':[4]},cells={'H':[[6,8]]})));monkeypatch.setattr(preprocess,'P',tmp_path);y,t,w=toy();p=apply(y,t,w,'P0',heat='H');assert p['invalid'][4].all();assert p['invalid'][6].sum()==1;assert p['invalid'][6,8]
def test_weighted_loss_changes_actual_fit():
 from scipy.linalg import svd
 rng=np.random.default_rng(9);y=rng.normal(size=(12,9))*5000+15000;t=np.arange(12.);w=np.arange(9.);mask=np.ones(y.shape,bool);weight=np.ones(y.shape);weight[3]=.25
 r=solve(y,t,w,mask,eta=0,rank=1,weights=weight,max_iter=1000,wall_seconds=30);L=r['latent']/65535;A=np.sqrt(weight)*(y/65535);u,s,v=svd(A,full_matrices=False);optimal=(u[:,:1]*s[:1])@v[:1]/np.sqrt(weight)
 assert np.linalg.norm(np.sqrt(weight)*(L-optimal))<.002

def test_preprocessing_records_complete():
 d=pd.read_csv(P/'results/frame_diagnostics_all.csv');assert len(d)==13589;assert d.heat_id.nunique()==10;assert not d.confirmed_acquisition_fault.any();a=json.loads((P/'config/final_preprocessing_policy.json').read_text());assert a['policy']=='P0'
 for split in ['development','confirmation']:
  d=pd.read_csv(P/f'results/synthetic_diagnosis_{split}.csv');assert d.kind.nunique()==6;assert set(d.level)=={'point','frame','event'}
