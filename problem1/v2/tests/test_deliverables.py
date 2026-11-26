import sys,json
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import P,HEATS,load,digest
from interpolation import load_surface

def test_locked_protocol():
 lock=json.loads((P/'config/LOCKED.json').read_text());assert digest(P/'config/protocol.json')==lock['protocol_sha256']
 for f in ['solver.py','detector.py','common.py']:assert digest(P/'src'/f)==lock['source_sha256'][f]
def test_complete_independent_design():
 d=pd.read_csv(P/'results/test_completion.csv');assert len(d)==3960;assert d.heat_id.nunique()==10;assert d.model.nunique()==6;assert d.scenario.nunique()==11;assert (d.n>0).all();assert np.isfinite(d.rmse).all()
 for (h,s,seed,scene),g in d.groupby(['heat_id','stage','seed','scenario']):assert g.n.nunique()==1
 d=pd.read_csv(P/'results/interpolation_test.csv');assert len(d)==600;assert set(d.method)=={'linear','cubic','pchip','tensor_bspline'}
def test_final_model_exports():
 total=0
 for h in HEATS:
  y,t,w=load(h);m=load_surface(P/f'data/{h}_model.npz');assert np.array_equal(m.time,t);assert np.array_equal(m.wave,w);assert np.isfinite(m.grid).all();assert np.isfinite(m.evaluate(t[1:3]+.1,w[10:12]+.01)).all();d=np.load(P/f'data/{h}_masks.npz');assert np.array_equal(d['saturated_unidentifiable'],y==65535);total+=d['saturated_unidentifiable'].sum();assert (P/f'figures/{h}_surface.png').exists()
 assert total==14805

def test_literature_correction():
 d=pd.read_csv(P/'references/literature_index_v2.csv').set_index('id');assert d.loc['V17','doi']=='10.1093/imaiai/iau006';assert (P/'references/V17_supplied_wrong_DOI_metadata.json').exists()

def test_censor_actual_constraints():
 d=pd.read_csv(P/'results/artificial_censoring.csv');assert len(d)==1800;assert (d[(d.method=='convex_nuclear_bound')&(d.n>0)].lower_bound_violation==0).all();assert (d.n>0).any()

def test_hidden_value_pipeline_invariance():
 from detector import diagnostics
 from solver import fit
 rng=np.random.default_rng(43);t=np.arange(20.);w=np.linspace(350,400,30);y=rng.uniform(100,1000,(20,30));hidden=rng.random(y.shape)<.1;other=y.copy();other[hidden]=65535
 outs=[]
 for z in [y,other]:
  visible=z.copy();visible[hidden]=np.nan;mask=diagnostics(visible,t,w,threshold=6,aux_threshold=10)[0];outs.append(fit(visible,t,w,mask,'B5',rank=4,max_iter=30)['grid'])
 assert np.array_equal(outs[0],outs[1])
