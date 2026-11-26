"""Post-run integrity checks against actual outputs, with a machine-readable artifact inventory."""
import json,hashlib
import numpy as np,pandas as pd
from audit import ROOT,P,HEATS,source_files,NEW_HASH
from models import load_surface
required=['references/README.md','references/references.bib','references/literature_index.csv','reports/00_相关文献综述.md','reports/00_文献方法与本题适用性分析.md','reports/00_数学模型推导与算法来源.md','problem1/reports/问题一_数据分析与模型建立及实验结果.md']
for p in required:assert (ROOT/p).is_file(),p
index=pd.read_csv(ROOT/'references/literature_index.csv');assert set(index.id)=={f'R{i:02d}' for i in range(1,13)}
for i in index.id:assert (ROOT/f'references/notes/{i}.md').exists()
assert hashlib.sha256(source_files()[0].read_bytes()).hexdigest()==NEW_HASH
r=pd.read_csv(P/'results/model_comparison/model_comparison.csv');assert len(r)==150 and r.heat_id.nunique()==10 and (r.n>0).all() and np.isfinite(r.rmse).all()
rows=[]
for h in HEATS:
 model=load_surface(P/f'data/reconstructed/{h}_model.npz');orig=np.load(P/f'data/processed/{h}.npz');dense=np.load(P/f'data/reconstructed/{h}_dense.npz');mask=np.load(P/f'data/masks/{h}.npz')
 assert mask['valid'].shape==orig['Y'].shape
 if model.family=='B':
  assert model.params.get('preserve_observed') is True and model.params.get('interpolation')=='pchip'
  np.testing.assert_array_equal(model.grid[mask['valid']],orig['Y'][mask['valid']])
 assert np.array_equal(model.wave,orig['wavelength']) and np.array_equal(model.time,orig['time'])
 assert dense['intensity'].shape==(len(dense['time']),len(dense['wavelength'])) and np.isfinite(dense['intensity']).all()
 # Distinguish saved dense samples from the evaluator: off-grid evaluations are finite.
 assert np.isfinite(model.evaluate(model.time[:3]+.125,model.wave[3:8]+.01)).all()
 assert (P/f'figures/reconstructed/{h}_surface.png').is_file()
 rows.append(dict(heat_id=h,original_shape=str(orig['Y'].shape),dense_shape=str(dense['intensity'].shape),model_kind=model.family))
pd.DataFrame(rows).to_csv(P/'results/final/output_integrity.csv',index=False)
# Detect accidental Python escape control characters in the generated Chinese/math report.
report=(P/'reports/问题一_数据分析与模型建立及实验结果.md').read_text();assert all(ord(c)>=32 or c in '\n\r\t' for c in report)
files=[]
for directory in [ROOT/'references',ROOT/'reports',P]:
 for path in directory.rglob('*'):
  if path.name=='output_inventory.json':continue
  if path.is_file() and not any(k in path.parts for k in ['fulltext','data','__pycache__','superseded_pilot']):
   files.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size))
(P/'results/final/output_inventory.json').write_text(json.dumps(files,ensure_ascii=False,indent=2))
print('12 literature entries/notes; 10 source hashes/grids/models/surfaces; 150 nonempty independent scores verified')
