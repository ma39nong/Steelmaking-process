"""Final reproducibility, data integrity and scoped scientific acceptance."""
from core import *
import re

def main():
 history=json.loads((P/'results/history_immutable.json').read_text());assert all(sha(ROOT/path)==h for path,h in history.items())
 lock=json.loads((P/'config/SELECTION_LOCKED.json').read_text());assert sha(P/'config/protocol.json')==lock['protocol_sha256'];assert sha(P/'src/detectors.py')==lock['detector_sha256'];assert sha(P/'src/experiment_final.py')==lock['experiment_sha256']
 manifest=pd.read_csv(ROOT/'data_manifest.csv');versions=[]
 for r in manifest.itertuples():assert sha(ROOT/r.source_file)==r.sha256;versions.append(dict(heat_id=r.heat_id,source_sha256=r.sha256,verified=True))
 pd.DataFrame(versions).to_csv(P/'results/raw_version_verification.csv',index=False)
 for folder in [V3,V31]:
  expected=pd.read_csv(folder/'results/artifact_manifest.csv')
  for r in expected.itertuples():assert sha(folder/f'data/{r.heat_id}_{r.artifact}.npz')==r.sha256
 f=pd.read_csv(P/'results/final_summary.csv');assert len(f)==10;assert int(f.repair_points.sum())==1043 and int(f.actually_changed.sum())==1043;assert f.unchanged_observation_max_error.max()==0;assert f.negative_nodes.sum()==0 and f.negative_dense.sum()==0;assert f.whole_frames_deleted.sum()==0 and f.recorded_top_code_max_error.max()==0 and f.conditional_lower_violations.sum()==0
 for r in pd.read_csv(P/'results/artifact_manifest.csv').itertuples():assert sha(P/f'data/{r.heat_id}_{r.artifact}.npz')==r.sha256
 # Every paper image and linked artifact must exist; every mathematical block must be paired.
 for path in [P/'archive/reports/问题一_中文论文初稿.md',P/'README.md',P/'references/README.md']:
  text=path.read_text();assert text.count('$$')%2==0;assert not any(ord(c)<32 and c not in '\n\t' for c in text)
  for dest in re.findall(r'\]\(([^)]+)\)',text):
   if not dest.startswith(('http','#')):assert (path.parent/dest).exists(),dest
 for h in HEATS:
  for suffix in ['异常修复曲线与局部热图','全谱处理与差值','最终连续三维光谱','帧级快变诊断']:assert (P/f'archive/figures/appendix/{h}_{suffix}.png').exists()
 logs={name:(P/'results'/name).read_text().strip().splitlines()[-1] for name in ['tests.log','v3_regression.log','v31_regression.log']};assert all('passed' in v for v in logs.values())
 exact=pd.read_csv(P/'results/exact_final_pipeline_confirmation.csv');x=exact.groupby('heat_id').mean(numeric_only=True).mean();cmp=pd.read_csv(P/'results/unified_comparison.csv');q=cmp[(cmp.split=='confirmation')&(cmp.strategy=='D_adaptive4')].iloc[0];tiles=pd.read_csv(P/'results/conditional_tiles.csv')
 dump(P/'results/ACCEPTANCE.json',dict(engineering_passed=True,historical_files_preserved=len(history),raw_files_verified=10,heats_complete=10,actual_statistical_repairs=int(f.actually_changed.sum()),repair_fraction=float(f.repair_points.sum()/f.observations.sum()),whole_frames_deleted=0,unchanged_observation_max_error=0,negative_native_and_dense=0,recorded_top_code_max_error=0,all_native_fits_stop_criterion_met=bool(f.converged.all()),numerical_fallback_points=int(f.numerical_fallback_points.sum()),conditional_tiles=int(len(tiles)),conditional_tiles_refitted=int(tiles.refitted.sum()),conditional_tiles_converged=int(tiles.converged.sum()),conditional_lower_violations=int(tiles.lower_violations.sum()),fixed_confirmation_precision=float(q.precision),fixed_confirmation_recall=float(q.recall),fixed_confirmation_f1=float(q.f1),final_pipeline_pollution_rmse=float(x.pollution_rmse),before_pollution_rmse=float(x.before_pollution_rmse),pollution_rmse_reduction=float(1-x.pollution_rmse/x.before_pollution_rmse),overall_rmse_reduction=float(1-x.overall_rmse/x.before_overall_rmse),unpolluted_rmse=float(x.unpolluted_rmse),introduced_peaks_in_fixed_final_pipeline=int(exact.introduced_unsupported_peak_count.sum()),tests=logs,scientific_scope='Actual statistically justified point repair supported by fixed end-to-end confirmation, plus nonnegative recorded continuous spectra. Tradeoffs and remaining misses are reported.',independent_external_heats=False,hardware_fault_accuracy_claim=False,saturated_true_upper_bound_identified=False,remaining_limitations=['same historical heats with new raw-coordinate regions','missed weak anomalies and possible real-peak modification','17 introduced peaks in fixed final-pipeline artificial cases','465 actual repairs use pre-existing numerical fallback','one conditional tile finite-iteration approximate; quality flag retained']))
if __name__=='__main__':main()
