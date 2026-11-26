from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np
P=Path(__file__).resolve().parents[1];V1=P.parent;V2=V1/'v2';ROOT=V1.parent
sys.path.insert(0,str(V1/'src'));sys.path.insert(0,str(V2/'src'))
from audit import HEATS
from anomalies import fill_training,detect
from experiment import metrics
CFG=json.loads((P/'config/protocol.json').read_text())
def external(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
legacy=external('legacy_v2_solver',V2/'src/solver.py')
convex=external('legacy_v2_convex',V2/'src/convex_censor.py')
from models import fit_b
from detector import diagnostics as v2_detect
from interpolation import Surface

def dump(path,x):Path(path).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=True))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(h):
 d=np.load(V1/f'data/processed/{h}.npz');return d['Y'],d['time'],d['wavelength']
def panels(role,split=None):
 for h in HEATS:
  y,t,w=load(h);starts=np.linspace(0,len(t)-96,6).astype(int)
  gaps=[1,2] if role=='preprocessing' else [0,4]
  for stage,k in enumerate(gaps):
   if split is not None and role=='preprocessing' and stage!=(0 if split=='development' else 1):continue
   lo=starts[k]+96;hi=starts[k+1];i=int((lo+hi-32)//2);j=int((.66 if role=='preprocessing' else [.28,.65][stage])*(len(w)-256))
   assert i>=lo and i+32<=hi
   yield h,stage,y[i:i+32,j:j+256].copy(),t[i:i+32],w[j:j+256],i,j

def hides(shape,seed):
 rng=np.random.default_rng(seed);n,m=shape;yield 'random_0.05',rng.random(shape)<.05
 q=np.zeros(shape,bool);i=int(rng.integers(5,n-8));q[i:i+4]=True;yield 'time_4',q
 q=np.zeros(shape,bool);j=int(rng.integers(10,m-18));q[:,j:j+8]=True;yield 'wave_8',q

def score(y,pred,hidden,w):
 q=hidden&np.isfinite(y)&(y>0)&(y<64000);r=metrics(y,pred,q,w,True);q2=hidden&np.isfinite(y)&(y>0)&(y<65535);r['all_unsaturated_rmse']=float(np.sqrt(np.mean((y[q2]-pred[q2])**2))) if q2.any() else np.nan;r['reference_fraction']=float(q.sum()/max(hidden.sum(),1));return r
