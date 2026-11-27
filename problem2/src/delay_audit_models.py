"""Descriptive reuse of existing fits; window samples are not independent events."""
import csv
import json
import numpy as np
from .delay_v2 import fit_delay_models
from .v2_validation import leave_one_heat_out


def assess(out,cfg):
    from .delay_audit import write_json,write_csv,compare_tree,ROOT
    paths=sorted((out/'baseline/data').glob('*_q2_processed.npz'))
    if len(paths)!=10: raise ValueError('all ten baseline interfaces required')
    samples={};pooled=[];folds=[]
    for gas in ['co','co2']:
        samples[gas]={}
        for p in paths:
            with np.load(p,allow_pickle=False) as z:
                centers=z[f'delay_{gas}_center_index'];conf=z[f'delay_{gas}_confidence'];delay=z[f'delay_{gas}_observed'];q=z['oxygen_flow'][centers]
                valid=(conf>0)&np.isfinite(delay)&(q>0)
                samples[gas][p.stem.split('_q2')[0]]=(q[valid],delay[valid],conf[valid])
        cat=lambda values,index: np.concatenate([v[index] for v in values])
        vs=list(samples[gas].values());q,y,w=(cat(vs,i) for i in range(3))
        model=fit_delay_models(q,y,w,cfg['delay']['max_lag_s'])
        pooled.append(dict(gas=gas,heats_with_windows=sum(len(v[0])>0 for v in vs),
                           flow_min_m3_s=float(q.min()),flow_max_m3_s=float(q.max()),
                           flow_quantiles_m3_s=np.quantile(q,[0,.25,.5,.75,1]),models=model,
                           evidence_unit='overlapping accepted windows; descriptive training fit only'))
        for heat,(testq,testy,testw) in samples[gas].items():
            train=[v for h,v in samples[gas].items() if h!=heat]
            tq,ty,tw=(cat(train,i) for i in range(3));m=fit_delay_models(tq,ty,tw,cfg['delay']['max_lag_s'])
            fixed=dynamic=None
            if len(testy) and 'fixed' in m:
                fixed=float(np.mean(abs(testy-m['fixed']['c_s'])))
                p=m['constrained_inverse_sqrt'];pred=np.clip(p['a_s']+p['b_s_sqrt_m3']/np.sqrt(testq),0,cfg['delay']['max_lag_s'])
                dynamic=float(np.mean(abs(pred-testy)))
            folds.append(dict(gas=gas,test_heat=heat,train_windows=len(ty),test_windows=len(testy),
                              fixed_test_mae_s=fixed,constrained_inverse_sqrt_test_mae_s=dynamic,
                              selected=m.get('selected'),
                              training_fixed_rmse_s=m.get('fixed',{}).get('rmse_s'),
                              training_dynamic_rmse_s=m.get('constrained_inverse_sqrt',{}).get('rmse_s')))
    reproduced=leave_one_heat_out(out/'baseline/data')
    historical=json.loads((ROOT/'problem2/final/results/cross_validation.json').read_text())['leave_one_heat_out']
    errors=compare_tree(reproduced,historical)
    if errors: raise AssertionError('historical LOHO mismatch: '+repr(errors))
    write_json(out/'flow_model_evidence.json',dict(pooled_training=pooled,
        historical_selected_model_loho_matches=True,
        historical_selected_model_loho=reproduced,
        independent_event_count=None,
        conclusion='insufficient evidence for a general decreasing flow-delay law; no independently confirmed events',
        selection_limitation='existing selection compares nested models using training RMSE; it is not evidence of a general flow effect'))
    write_csv(out/'flow_model_cross_validation.csv',folds)
    accepted=list(csv.DictReader((out/'accepted_windows.csv').open()))
    coverage=[]
    for gas in ['co','co2']:
        for heat in sorted(samples[gas]):
            rs=[r for r in accepted if r['gas']==gas and r['heat_id']==heat]
            coverage.append(dict(heat_id=heat,gas=gas,accepted_windows=len(rs),
                                  stages=';'.join(sorted({r['stage'] for r in rs})),
                                  flow_min_m3_s=min([float(r['flow_m3_s']) for r in rs],default=None),
                                  flow_max_m3_s=max([float(r['flow_m3_s']) for r in rs],default=None),
                                  source_overlap_groups=len({r['source_overlap_group'] for r in rs}),
                                  search_overlap_groups=len({r['search_overlap_group'] for r in rs}),
                                  cross_gas_same_center_windows=sum(int(r['other_gas_same_center_count'])>0 for r in rs),
                                  confirmed_independent_events=None))
    write_csv(out/'flow_stage_coverage.csv',coverage)
