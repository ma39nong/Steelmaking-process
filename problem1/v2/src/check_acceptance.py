"""Engineering completeness checks are independent from scientific performance acceptance."""
import json
import numpy as np,pandas as pd
from common import *
required=['reports/00_V1模型审计.md','reports/01_数学推导.md','reports/问题一_V2模型深化与严格验证.md','references/literature_index_v2.csv','results/TEST_COMPLETE.json','results/FINAL_COMPLETE.json','results/FULL_CONFIRMATION_COMPLETE.json','results/FULL_INTERPOLATION_COMPLETE.json','results/tests_final.log']
for f in required:assert (P/f).exists(),f
for f,s in json.loads((P/'results/v1_immutable_manifest.json').read_text()).items():assert digest(ROOT/f)==s,f
for h,s in json.loads((P/'results/source_cache_sha256.json').read_text()).items():assert digest(V1/f'data/processed/{h}.npz')==s,h
lock=json.loads((P/'config/LOCKED.json').read_text());assert digest(P/'config/protocol.json')==lock['protocol_sha256']
for f in ['solver.py','detector.py','common.py']:assert digest(P/'src'/f)==lock['source_sha256'][f],f
assert digest(P/'src/interpolation.py')==json.loads((P/'config/INTERPOLATION_LOCKED.json').read_text())['source_sha256']
d=pd.read_csv(P/'results/test_completion.csv');assert len(d)==3960
cols=['heat_id','stage','seed','scenario'];a=d[d.model=='B2'].set_index(cols);b=d[d.model=='B3'].set_index(cols);assert np.array_equal(a.rmse,b.rmse);assert (b.valid_observation_rmse_not_prediction==0).all()
full=pd.read_csv(P/'results/full_confirmation.csv');assert len(full)==550 and full.heat_id.nunique()==10
cont=pd.read_csv(P/'results/full_interpolation.csv');assert len(cont)==250 and cont.heat_id.nunique()==10
for h in HEATS:
 for f in [f'data/{h}_model.npz',f'data/{h}_masks.npz',f'figures/{h}_surface.png',f'figures/{h}_masks.png']:assert (P/f).exists(),f
assert '37 passed' in (P/'results/tests_final.log').read_text()
report=(P/'reports/问题一_V2模型深化与严格验证.md').read_text();assert all(f'## {i} ' in report for i in range(1,12));assert '228.9%' in report
assert not any(ord(c)<32 and c not in '\n\t\r' for c in report)
status=json.loads((P/'results/acceptance_status.json').read_text());assert status['scientific_all_requirements_pass'] is False
summary=dict(engineering_checks_pass=True,scientific_performance_all_pass=False,v1_tracked_files_preserved=len(json.loads((P/'results/v1_immutable_manifest.json').read_text())),native_completion_rows=3960,c_rows=660,v1_matched_local_rows=660,full_native_confirmation_rows=550,local_interpolation_rows=600,full_interpolation_rows=250,artificial_censoring_rows=1800,tests_passed=37,B2_B3_hidden_rmse_identical=True,report_sha256=digest(P/'reports/问题一_V2模型深化与严格验证.md'))
dump(P/'results/final_engineering_check.json',summary);dump(P/'results/final_source_manifest.json',{str(p.relative_to(P)):digest(p) for d in ['src','config','tests'] for p in (P/d).glob('*') if p.is_file()});print(json.dumps(summary,ensure_ascii=False,indent=2))
