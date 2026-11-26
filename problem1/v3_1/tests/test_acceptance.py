import sys,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from common31 import ROOT,HEATS,panels,load
from points import evidence,mask_for,repair
from experiments import stats,inject
from interface import Spectrum

def test_isolated_spike_and_persistent_line():
 t=np.arange(32)*.5;w=np.linspace(400,500,64);y=500+2*t[:,None]+np.cos(w[None,:]);y[15,20]+=5000
 e=evidence(y,t,w);assert e['high'][15,20]
 y[:,40]+=3000;e=evidence(y,t,w);assert not e['high'][:,40].any()
 D,_,_=repair(y,t,w,e,'conservative');assert np.array_equal(D[~e['high']],y[~e['high']]);assert D[15,20]<y[15,20];assert (D>=0).all()

def test_top_codes_and_neighborhood_never_automatically_repaired():
 t=np.arange(32)*.5;w=np.arange(64);y=np.full((32,64),500.);y[15,20]=65535;y[15,21]=5000
 e=evidence(y,t,w);assert not e['high'][15,20];assert not e['high'][15,21]

def test_absolute_false_positives_not_reference_excluded():
 ref=np.ones((8,8))*100;corrupt=ref.copy();corrupt[2,2]=200;pred=ref.copy();pred[3,3]=110;truth=corrupt!=ref;mask=truth.copy();mask[3,3]=True
 s=stats(ref,corrupt,pred,truth,mask,np.arange(8));assert s['tp']==1 and s['fp']==1 and s['precision']==.5;assert s['unpolluted_rmse']>0

def test_frozen_rules_and_history():
 protocol=json.loads((P/'config/PROTOCOL_LOCKED.json').read_text());lock=json.loads((P/'config/SELECTION_LOCKED.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 assert sha(P/'config/protocol.json')==protocol['protocol_sha256']==lock['protocol_sha256'];assert sha(P/'src/points.py')==protocol['detector_sha256']==lock['points_source_sha256'];assert sha(P/'src/experiments.py')==lock['experiments_source_sha256'];assert sha(P/'src/select.py')==lock['selection_source_sha256']
 for path,h in json.loads((P/'results/history_immutable.json').read_text()).items():assert sha(ROOT/path)==h,path


def test_new_regions_disjoint_from_old_local_panels():
 old3=pd.read_csv(P.parent/'v3/results/new_region_index.csv')
 for split in ['development','confirmation']:
  for h,y,t,w,i,j in panels(split):
   raw,tt,ww=load(h);starts=np.linspace(0,len(tt)-96,6).astype(int);assert all(i+32<=s or i>=s+96 for s in starts)
   for r in old3[old3.heat_id==h].itertuples():assert i+31<r.row_start or i>r.row_end or j+255<r.column_start or j>r.column_end


def test_point_lists_partition_and_review_veto():
 a=pd.read_csv(P/'results/raw_point_audit.csv');rle=pd.read_csv(P/'results/retained_points_rle.csv');high=pd.read_csv(P/'results/high_confidence_statistical_points.csv');sus=pd.read_csv(P/'results/general_suspect_points.csv');assert len(high)==4
 for row in a.itertuples():
  assert rle[rle.heat_id==row.heat_id]['count'].sum()+len(high[high.heat_id==row.heat_id])+len(sus[sus.heat_id==row.heat_id])==row.observations
 veto=sus[(sus.heat_id=='6C02108')&(sus.row==608)&(sus.column==1211)];assert len(veto)==1 and veto.native_review.iloc[0]=='retain_persistent_line'
 assert not high.hardware_fault_confirmed.any()


def test_confirmation_failure_honestly_falls_back():
 frozen=json.loads((P/'config/FINAL_FROZEN.json').read_text());g=pd.read_csv(P/'results/selection_gates_confirmation.csv');assert not frozen['confirmation_passed'];assert frozen['strategy']=='V3-P0';assert not g[g.strategy=='conservative'].passed_added_processing.iloc[0]


def test_ten_final_interfaces_and_nonnegative_surfaces():
 df=pd.read_csv(P/'results/final_summary.csv');assert len(df)==10;assert df.repaired_points.sum()==0
 for h in HEATS:
  d=np.load(P/f'data/{h}_interface.npz');raw,t,w=load(h);assert np.array_equal(raw,d['raw']);assert np.array_equal(raw,d['reconstructed']);assert not d['statistical_repair'].any();assert np.isfinite(d['reconstructed']).all();assert (d['reconstructed']>=0).all()
  sat=d['suspected_saturation'];assert (d['reconstructed'][sat]==65535).all();assert (d['conditional'][sat]>=65535).all()
  model=Spectrum.load(P/f'data/{h}_interface.npz');assert np.array_equal(model.evaluate(t[::100],w[::100]),raw[::100,::100]);assert not model.raw_until(2)['offline_reconstruction_included']
  dense=np.load(P/f'data/{h}_dense.npz')['recorded'];assert np.isfinite(dense).all() and (dense>=0).all()
