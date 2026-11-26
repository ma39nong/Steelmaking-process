"""Apply the predeclared gates once. Confirmation cannot retune thresholds or model."""
from common31 import *
from datetime import datetime,timezone

def assess(split):
 d=pd.read_csv(P/f'results/contamination_{split}.csv');control=pd.concat([pd.read_csv(P/f'results/control_{split}.csv'),pd.read_csv(P/f'results/persistent_line_control_{split}.csv')]);gates=CFG['selection_rule'];out=[]
 for name in CFG['strategies']:
  q=d[d.strategy==name];c=control[control.strategy==name];s=q[['tp','fp','fn','tn']].sum();precision=float(s.tp/max(s.tp+s.fp,1));recall=float(s.tp/max(s.tp+s.fn,1));macro=q.groupby('heat_id').mean(numeric_only=True).mean();by=q.groupby(['kind','strength']).mean(numeric_only=True)
  details=dict(strategy=name,precision_micro=precision,recall_micro=recall,f1_micro=2*precision*recall/max(precision+recall,1e-12),pollution_rmse=float(macro.pollution_rmse),before_pollution_rmse=float(macro.before_pollution_rmse),pollution_reduction=float(1-macro.pollution_rmse/macro.before_pollution_rmse),overall_reduction=float(1-macro.overall_rmse/macro.before_overall_rmse),unpolluted_changed_fraction=float(macro.unpolluted_changed_fraction),unpolluted_to_before_ratio=float(macro.unpolluted_rmse/macro.before_overall_rmse),clean_control_peak_p95_max=float(c.peak_p95.max()),clean_control_mask_fraction_max=float(c.masked_fraction.max()),introduced_unsupported_peak_count=int(q.introduced_unsupported_peak_count.sum()),all_kind_strength_not_worse=bool((by.overall_rmse<=by.before_overall_rmse+1e-8).all()))
  checks=[precision>=gates['minimum_precision_micro'],recall>=gates['minimum_recall_micro'],details['pollution_reduction']>=gates['minimum_pollution_rmse_reduction'],details['overall_reduction']>=gates['minimum_overall_rmse_reduction'],details['unpolluted_changed_fraction']<=gates['maximum_unpolluted_changed_fraction'],details['unpolluted_to_before_ratio']<=gates['maximum_unpolluted_to_before_rmse_ratio'],details['clean_control_peak_p95_max']<=gates['maximum_clean_control_peak_p95'],details['clean_control_mask_fraction_max']<=gates['maximum_clean_control_mask_fraction'],details['introduced_unsupported_peak_count']==0,details['all_kind_strength_not_worse']]
  details['passed_added_processing']=bool(all(checks));details['checks']=checks;out.append(details)
 dump(P/f'results/selection_gates_{split}.json',out);pd.DataFrame([{k:v for k,v in x.items() if k!='checks'} for x in out]).to_csv(P/f'results/selection_gates_{split}.csv',index=False);return out

def main(split):
 out=assess(split)
 if split=='development':
  assert not (P/'config/SELECTION_LOCKED.json').exists();eligible=[x for x in out if x['strategy']!='V3-P0' and x['passed_added_processing']];name=min(eligible,key=lambda x:(x['unpolluted_changed_fraction'],x['pollution_rmse']))['strategy'] if eligible else 'V3-P0'
  dump(P/'config/SELECTION_LOCKED.json',dict(strategy=name,utc=datetime.now(timezone.utc).isoformat(),confirmation_unopened=True,protocol_sha256=sha(P/'config/protocol.json'),points_source_sha256=sha(P/'src/points.py'),experiments_source_sha256=sha(P/'src/experiments.py'),selection_source_sha256=sha(P/'src/select.py'),review_ledger_sha256=sha(P/'config/native_review_vetoes.json'),rule='all gates then least unpolluted change; conservative original proposals additionally veto-reviewed, never added'))
 else:
  lock=json.loads((P/'config/SELECTION_LOCKED.json').read_text());selected=next(x for x in out if x['strategy']==lock['strategy']);accepted=lock['strategy']=='V3-P0' or selected['passed_added_processing'];name=lock['strategy'] if accepted else 'V3-P0';dump(P/'config/FINAL_FROZEN.json',dict(strategy=name,development_strategy=lock['strategy'],confirmation_passed=bool(accepted),utc=datetime.now(timezone.utc).isoformat(),rank=4,eta=.003,interpolation='bilinear',sat='unchanged V3',point_rule='points.py and native veto review',hardware_faults_confirmed=0,further_parameter_search=False))
 print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main(sys.argv[1])
