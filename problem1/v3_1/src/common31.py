from pathlib import Path
import sys,json,hashlib,subprocess
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];V3=P.parent/'v3';ROOT=P.parents[1]
sys.path.insert(0,str(V3/'src'))
from ctx import HEATS,load,legacy,fill_training,detect,v2_detect
from model3 import solve,project_record
from interface import Spectrum
CFG=json.loads((P/'config/protocol.json').read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=True))
def panels(split):
 for h in HEATS:
  y,t,w=load(h);old=np.linspace(0,len(t)-96,6).astype(int);i=(old[3]+96+old[4]-32)//2
  assert old[3]+96<=i and i+32<=old[4]
  for j in ([64] if split=='development' else [900,1720]):yield h,y[i:i+32,j:j+256].copy(),t[i:i+32],w[j:j+256],int(i),j
