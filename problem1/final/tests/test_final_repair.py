import sys,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from core import ROOT,V3,V31,HEATS,panels,load
from detectors import adaptive,local_mask,reconstruct
from interface import Spectrum

def test_adaptive_repairs_isolated_but_keeps_persistent_line():
 t=np.arange(32)*.5;w=np.linspace(400,500,64);y=500+2*t[:,None]+np.cos(w[None,:]);y[15,20]+=5000;y[:,40]+=3000
 e=adaptive(y,t,w);mask=local_mask(e,4);assert mask[15,20];assert not mask[:,40].any();D,_=reconstruct(y,t,w,mask);assert D[15,20]<y[15,20];assert np.array_equal(D[~mask],y[~mask])

def test_coherent_frame_brightness_not_deleted():
 t=np.arange(32)*.5;w=np.linspace(400,500,64);y=np.ones((32,64))*1000;y[15]*=2
 mask=local_mask(adaptive(y,t,w),4);assert not mask.any()

def test_top_code_separate_and_domain_edges_diagnostic():
 t=np.arange(32)*.5;w=np.arange(64.);y=np.ones((32,64))*500;y[15,20]=65535;y[15,21]=5000;y[1,40]=5000
 mask=local_mask(adaptive(y,t,w),4);assert not mask[15,20] and not mask[15,21] and not mask[1,40]

def test_history_and_numeric_locks_unchanged():
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();lock=json.loads((P/'config/SELECTION_LOCKED.json').read_text());assert sha(P/'src/detectors.py')==lock['detector_sha256'];assert sha(P/'src/experiment_final.py')==lock['experiment_sha256'];assert sha(P/'config/protocol.json')==lock['protocol_sha256']
 # The frozen manifest describes the original Git snapshot, including editable README/AGENTS.
 # Verify that snapshot itself, then require live historical scientific artifacts to remain unchanged.
 import subprocess,tarfile,io
 frozen=json.loads((P/'results/history_immutable.json').read_text());base=json.loads((P/'config/protocol.json').read_text())['base_commit']
 blob=subprocess.check_output(['git','archive',base],cwd=ROOT)
 with tarfile.open(fileobj=io.BytesIO(blob)) as archive:
  members={m.name:m for m in archive.getmembers() if m.isfile()}
  assert set(frozen)<=set(members)
  for path,h in frozen.items():
   raw=archive.extractfile(members[path]).read();variants=[raw]
   # CSV checkout bytes may be CRLF although Git stores LF (root .gitattributes).
   if Path(path).suffix=='.csv':variants.append(raw.replace(b'\r\n',b'\n').replace(b'\n',b'\r\n'))
   assert h in {hashlib.sha256(v).hexdigest() for v in variants},path
   if path.startswith('problem1/') and not path.startswith('problem1/final/') and Path(path).suffix in {'.py','.json','.csv','.npz'}:
    assert sha(ROOT/path)==h,path

def test_new_regions_disjoint_in_raw_coordinates():
 prior=pd.concat([pd.read_csv(V3/'results/new_region_index.csv'),pd.read_csv(V31/'results/region_index.csv')])
 for split in ['development','confirmation']:
  for h,b,y,t,w,i,j in panels(split):
   raw,tt,ww=load(h);starts=np.linspace(0,len(tt)-96,6).astype(int);assert all(i+32<=s or i>=s+96 for s in starts)
   for r in prior[prior.heat_id==h].itertuples():assert i+31<r.row_start or i>r.row_end or j+127<r.column_start or j>r.column_end

def test_four_reviewed_points_reidentified_and_persistent_veto_retained():
 review=pd.read_csv(V31/'results/high_confidence_statistical_points.csv')
 for r in review.itertuples():
  m=np.load(P/f'data/{r.heat_id}_masks.npz');assert m['repair'][r.row,r.column];assert m['automatic'][r.row,r.column]
 assert not np.load(P/'data/6C02108_masks.npz')['repair'][608,1211]

def test_all_native_outputs_and_dense_nonnegative():
 summary=pd.read_csv(P/'results/final_summary.csv');assert len(summary)==10;assert summary.repair_points.sum()==1043;assert summary.actually_changed.sum()>0
 for h in HEATS:
  d=np.load(P/f'data/{h}_interface.npz');raw,t,w=load(h);mask=d['statistical_repair'];D=d['reconstructed'];assert np.array_equal(raw,d['raw']);assert np.array_equal(D[~mask],raw[~mask]);assert (D>=0).all() and np.isfinite(D).all();assert not mask.all(1).any();sat=d['suspected_saturation'];assert (D[sat]==65535).all();assert (d['conditional'][sat]>=65535).all()
  assert np.array_equal(np.load(P/f'data/{h}_difference.npz')['difference'],D-raw);s=Spectrum.load(P/f'data/{h}_interface.npz');assert np.array_equal(s.evaluate(t[::100],w[::100]),D[::100,::100]);assert not s.raw_until(2)['offline_reconstruction_included'];dense=np.load(P/f'data/{h}_dense.npz')['recorded'];assert np.isfinite(dense).all() and (dense>=0).all()

def test_comparison_case_pairing():
 d=pd.read_csv(P/'results/contamination_confirmation.csv');assert len(d)==1200;assert d.groupby(['heat_id','band','seed','kind','strength']).strategy.nunique().eq(5).all();assert set(d.strategy)=={'A_P0','B_V1_MAD','C_V31','D_adaptive4','D_adaptive6'}
