import subprocess,json
import pandas as pd
from ctx import *
assert subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()==CFG['base_commit']
files=subprocess.check_output(['git','ls-tree','-r','--name-only','-z',CFG['base_commit']]).decode().rstrip('\0').split('\0');dump(P/'results/history_immutable.json',{f:sha(ROOT/f) for f in files})
rows=[]
for role in ['preprocessing','model']:
 for h,s,y,t,w,i,j in panels(role):rows.append(dict(role=role,heat_id=h,stage=s,row_start=i,row_end=i+31,column_start=j,column_end=j+255,time_start=t[0],time_end=t[-1],wave_start=w[0],wave_end=w[-1],panel_sha256=hashlib.sha256(y.astype('<f8').tobytes()).hexdigest()))
pd.DataFrame(rows).to_csv(P/'results/new_region_index.csv',index=False)
dump(P/'config/PROTOCOL_LOCKED.json',dict(protocol_sha256=sha(P/'config/protocol.json'),base_commit=CFG['base_commit'],utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),stage='before any V3 scores'))
