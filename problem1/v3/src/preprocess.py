import json
import numpy as np
from ctx import *
from frame_review import frame_metrics

def apply(y,t,w,policy='P0',heat=None,row_offset=0,col_offset=0):
 finite=np.isfinite(y)&(y>=0);sat=finite&(y==65535);confirmed=~finite;ledger=json.loads((P/'config/confirmed_invalid.json').read_text())
 if heat:
  for row in ledger['frames'].get(heat,[]):
   i=int(row)-row_offset
   if 0<=i<len(y):confirmed[i]=True
  for row,col in ledger['cells'].get(heat,[]):
   i=int(row)-row_offset;j=int(col)-col_offset
   if 0<=i<len(y) and 0<=j<y.shape[1]:confirmed[i,j]=True
 frame=frame_metrics(np.where(finite,y,np.nan),t,w)
 # Whole-frame structural unavailability alone is enough; statistical candidates remain unconfirmed.
 confirmed[frame.structural_invalid.to_numpy()]=True
 candidate_local=~v2_detect(np.where(finite,y,np.nan),t,w,threshold=6,aux_threshold=10)[0]&finite&~sat
 local=np.zeros(y.shape,bool) if policy=='P0' else candidate_local
 invalid=confirmed|local
 if policy=='P4':invalid[frame.candidate.to_numpy()]=True
 valid=finite&~invalid&~sat
 weights=np.ones(y.shape)
 if policy=='P3':weights[frame.candidate.to_numpy()]=.25
 weights[~valid]=0
 types=np.zeros(y.shape,np.uint16);types[~finite]|=1;types[sat]|=2;types[candidate_local]|=4;types[confirmed&finite]|=8;types[frame.candidate.to_numpy()]|=16
 grade=np.zeros(y.shape,np.uint8);grade[candidate_local]=1;grade[frame.candidate.to_numpy()]=1;grade[invalid]=2;grade[sat]=1
 return dict(valid=valid,invalid=invalid,sat=sat,weights=weights,types=types,grade=grade,frames=frame,candidate_rows=frame.candidate.to_numpy())
