"""Assemble complete independent-heat caches from serial and isolated worker runs."""
import pandas as pd,json
from common import *
rows=[]
for h in HEATS:
 path=P/f'results/{h}_full_confirmation.json'
 if not path.exists():raise RuntimeError('Not complete: '+h)
 local=json.loads(path.read_text());assert len(local)==55 and {r['heat_id'] for r in local}=={h};rows.extend(local)
df=pd.DataFrame(rows);assert len(df)==550
assert not df.duplicated(['heat_id','model','seed','scenario']).any()
df.to_csv(P/'results/full_confirmation.csv',index=False)
paired=df[df.seed==811].pivot(index=['heat_id','scenario'],columns='model',values='rmse');paired['B2_minus_V1']=paired.B2-paired.V1_B3_matched;paired.to_csv(P/'results/full_confirmation_paired_seed811.csv')
dump(P/'results/FULL_CONFIRMATION_COMPLETE.json',dict(rows=len(rows),heat_count=10,config_sha256=digest(P/'config/full_confirmation.json'),assembly='Five serial caches and five isolated worker caches; same code/parameters/seeds. Only complete caches assembled.'))
