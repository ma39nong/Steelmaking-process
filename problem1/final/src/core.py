from pathlib import Path
import sys,json,hashlib,subprocess
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];V1=P.parent;ROOT=V1.parent;V3=V1/'v3';V31=V1/'v3_1'
sys.path.insert(0,str(V3/'src'))
from ctx import HEATS,load,legacy,detect,fill_training
from model3 import solve,project_record
from stable_completion import complete_with_stability
from interface import Spectrum
sys.path.insert(0,str(V31/'src'))
from points import evidence as conservative_evidence
from experiments import inject,stats
CFG=json.loads((P/'config/protocol.json').read_text())
V1PAR=json.loads((V1/'results/final/model_selection.json').read_text())['per_heat']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=True))
def panels(split):
 for h in HEATS:
  y,t,w=load(h);old=np.linspace(0,len(t)-96,6).astype(int);k=1 if split=='development' else 2;i=int((old[k]+96+old[k+1]-32)//2)
  for band,j in zip(['blue','middle','red'],[64,1050,1740]):yield h,band,y[i:i+32,j:j+128].copy(),t[i:i+32],w[j:j+128],i,j
