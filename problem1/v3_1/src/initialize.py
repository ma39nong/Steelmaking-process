from common31 import *
from datetime import datetime,timezone
assert subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()==CFG['base_commit']
assert not (P/'config/PROTOCOL_LOCKED.json').exists()
files=subprocess.check_output(['git','ls-tree','-r','--name-only','-z','HEAD']).decode().rstrip('\0').split('\0')
dump(P/'results/history_immutable.json',{f:sha(ROOT/f) for f in files})
rows=[]
for split in ['development','confirmation']:
 for h,y,t,w,i,j in panels(split):rows.append(dict(split=split,heat_id=h,row_start=i,row_end=i+31,column_start=j,column_end=j+255,time_start=t[0],time_end=t[-1],wave_start=w[0],wave_end=w[-1],panel_sha256=hashlib.sha256(y.astype('<f8').tobytes()).hexdigest()))
pd.DataFrame(rows).to_csv(P/'results/region_index.csv',index=False)
dump(P/'config/PROTOCOL_LOCKED.json',dict(utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=sha(P/'config/protocol.json'),detector_sha256=sha(P/'src/points.py'),base_commit=CFG['base_commit'],before_all_V31_scores=True))
