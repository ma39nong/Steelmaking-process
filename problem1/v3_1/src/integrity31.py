"""Separate engineering success from failed actual anomaly-removal confirmation."""
from common31 import *
import re

def main():
 historical=json.loads((P/'results/history_immutable.json').read_text());assert all(sha(ROOT/f)==h for f,h in historical.items())
 lock=json.loads((P/'config/SELECTION_LOCKED.json').read_text());assert sha(P/'config/protocol.json')==lock['protocol_sha256'];assert sha(P/'src/points.py')==lock['points_source_sha256'];assert sha(P/'src/experiments.py')==lock['experiments_source_sha256'];assert sha(P/'src/select.py')==lock['selection_source_sha256'];assert sha(P/'config/native_review_vetoes.json')==lock['review_ledger_sha256']
 old_artifacts=pd.read_csv(V3/'results/artifact_manifest.csv')
 for r in old_artifacts.itertuples():assert sha(V3/f'data/{r.heat_id}_{r.artifact}.npz')==r.sha256
 expected=pd.read_csv(V3/'results/data_version_verification.csv');sources=[]
 for r in expected.itertuples():assert sha(ROOT/r.source_file)==r.sha256;sources.append(dict(heat_id=r.heat_id,source_sha256_verified=True))
 pd.DataFrame(sources).to_csv(P/'results/data_version_verification.csv',index=False)
 manifests=[]
 for h in HEATS:
  d=np.load(P/f'data/{h}_interface.npz');v=np.load(V3/f'data/{h}_interface.npz');assert np.array_equal(d['raw'],v['raw']);assert np.array_equal(d['reconstructed'],d['raw']);assert np.array_equal(d['conditional'],v['conditional']);assert not d['statistical_repair'].any();assert np.isfinite(d['reconstructed']).all();assert (d['reconstructed']>=0).all();assert (d['difference_from_raw']==0).all()
  assert (P/f'figures/appendix/{h}_最终连续光谱.png').exists();assert (P/f'figures/appendix/{h}_处理前后与局部热图.png').exists();assert (P/f'figures/appendix/{h}_全矩阵差值.png').exists()
  for artifact in ['interface','dense','difference']:path=P/f'data/{h}_{artifact}.npz';manifests.append(dict(heat_id=h,artifact=artifact,sha256=sha(path),bytes=path.stat().st_size))
 pd.DataFrame(manifests).to_csv(P/'results/artifact_manifest.csv',index=False)
 a=pd.read_csv(P/'results/raw_point_audit.csv');f=pd.read_csv(P/'results/final_summary.csv');c=pd.read_csv(P/'results/counterfactual_repairs.csv');assert len(c)==4;assert c.within_neighbor_support.all();assert not c.formal_repair_applied.any();assert f.repaired_points.sum()==0 and f.whole_frames_deleted.sum()==0
 # Resolve every local Markdown destination in the two final reports.
 for path in (P/'reports').glob('*.md'):
  for dest in re.findall(r'\]\(([^)]+)\)',path.read_text()):
   if not dest.startswith(('http','#')):assert (path.parent/dest).exists(),dest
 for name in ['tests.log','v3_regression_tests.log']:assert 'passed' in (P/'results'/name).read_text()
 frozen=json.loads((P/'config/FINAL_FROZEN.json').read_text());g=pd.read_csv(P/'results/selection_gates_confirmation.csv');cg=g[g.strategy=='conservative'].iloc[0]
 dump(P/'results/ACCEPTANCE.json',dict(engineering_passed=True,historical_files_unchanged=len(historical),raw_files_verified=10,final_heats=10,statistical_candidates=int(a.statistical_candidates.sum()),automatic_strong_proposals=5,review_supported_points=4,persistent_line_vetoes=1,formal_points_repaired=int(f.repaired_points.sum()),whole_frames_deleted=0,formal_negative_count=0,scientific_actual_anomaly_removal_passed=False,formal_scope='anomaly diagnosis and faithful continuous reconstruction; NOT completed actual anomaly elimination',development_strategy='conservative',confirmation_recall=float(cg.recall_micro),confirmation_precision=float(cg.precision_micro),confirmation_pollution_rmse_reduction=float(cg.pollution_reduction),fallback_frozen=frozen['strategy'],thresholds_retuned_after_confirmation=False,engineering_tests=(P/'results/tests.log').read_text().strip().splitlines()[-1],V3_regression_tests=(P/'results/v3_regression_tests.log').read_text().strip().splitlines()[-1],limitations=['same ten historical heats, new regions not external data','reference raw spectra are noisy','very low confirmation recall','persistent-line ambiguity','four supported native statistical points retained despite counterfactual feasibility; deployment gates not met']))
if __name__=='__main__':main()
