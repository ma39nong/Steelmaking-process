from pathlib import Path
import sys,json,hashlib
import numpy as np
P=Path(__file__).resolve().parents[1]; V1=P.parent; ROOT=V1.parent
sys.path.insert(0,str(V1/'src'))
from audit import HEATS
from anomalies import fill_training,detect
from experiment import metrics
CFG=json.loads((P/'config/protocol.json').read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=True))
def load(h):
 d=np.load(V1/f'data/processed/{h}.npz');return d['Y'],d['time'],d['wavelength']
def panels(split):
 for h in HEATS:
  y,t,w=load(h);n,m=CFG['panel_shape'];starts=np.linspace(0,len(t)-n,6).astype(int)
  for stage,k in enumerate(range(0 if split=='development' else 1,6,2)):
   i=int(starts[k]);j=int([.08,.40,.72][stage]*(len(w)-m));yield h,stage,y[i:i+n,j:j+m].copy(),t[i:i+n],w[j:j+m]
def hides(shape,seed):
 rng=np.random.default_rng(seed);n,m=shape
 for f in CFG['random_fractions']:yield f'random_{f:g}',rng.random(shape)<f
 for width in CFG['time_widths']:
  q=np.zeros(shape,bool);i=int(rng.integers(4,n-width-3));q[i:i+width]=True;yield f'time_{width}',q
 for width in CFG['wave_widths']:
  q=np.zeros(shape,bool);j=int(rng.integers(4,m-width-3));q[:,j:j+width]=True;yield f'wave_{width}',q

def score(y,pred,hidden,w):
 q=hidden&np.isfinite(y)&(y>0)&(y<64000);r=metrics(y,pred,q,w,True)
 q2=hidden&np.isfinite(y)&(y>0)&(y<65535);r['all_unsaturated_rmse']=float(np.sqrt(np.mean((pred[q2]-y[q2])**2))) if q2.any() else np.nan
 r['reference_retained_fraction']=float(q.sum()/max(hidden.sum(),1))
 r['valid_observation_rmse_not_prediction']=np.nan
 return r
