"""Conditional real-data exploration. Called only after the synthetic gate passes."""
import copy
import numpy as np
from ...src.delay_audit import write_csv,write_json
from ..signal_fidelity.run import load_heats
from .model import Design,bootstrap,predict


def execute(cfg,out):
    heats,manifest=load_heats();write_json(out/'input_manifest.json',manifest)
    rows=[];profiles=[];objects={};draws=[]
    for hi,z in enumerate(heats):
        t=z['time'];end=t[-1];local=copy.deepcopy(cfg);p=local['integral_probe']
        train_end=np.floor(.6*end/5)*5
        val_start=np.ceil(max(.7*end,train_end+65)*2)/2;val_end=np.floor((end-1)*2)/2
        p['training_bounds_s']=[60,train_end];p['validation_bounds_s']=[val_start,val_end]
        assert val_end-val_start>=50
        for gas in ['co','co2']:
            o=z['oxygen_cumulative_raw'];y=z[gas+'_raw']
            assert np.isfinite(o).all() and np.isfinite(y).all()
            for method in ['ordinary','weak']:
                design=Design(t,o,local,method);r,co,loss,cond=design.fit(y)
                ci,samples=bootstrap(design,y,r,local,cfg['bootstrap_seed']+hi*10+(gas=='co2'),block=True)
                r.update(ci);r['accepted_with_interval']=bool(r['accepted'] and ci['ci_width_s']<=5)
                label=dict(heat_id=z['heat_id'],gas=gas,method=method,train_start_s=60,train_end_s=train_end,
                           validation_start_s=val_start,validation_end_s=val_end,
                           mean_source_flow_m3_s=float((np.interp(train_end,t,o)-np.interp(60,t,o))/(train_end-60)),
                           input='raw cumulative oxygen and raw recorded gas; no SG derivative or Hampel modification')
                rows.append(dict(**label,**r));objects[(z['heat_id'],gas,method)]=(z,local,r)
                for lag,cc,ll,cn in zip(design.lags,co,loss,cond):profiles.append(dict(heat_id=z['heat_id'],gas=gas,method=method,lag_s=lag,loss=ll,T_s=cc[0],gain=cc[1],offset=cc[2],scaled_design_condition=cn))
                for j,s in enumerate(samples):draws.append(dict(heat_id=z['heat_id'],gas=gas,method=method,replicate=j,lag_s=s[0],T_s=s[1],gain=s[2],offset=s[3]))
            print('real',z['heat_id'],gas,flush=True)
    loho=[]
    for gas in ['co','co2']:
        for held in heats:
            donors=[objects[(z['heat_id'],gas,'weak')] for z in heats if z['heat_id']!=held['heat_id']]
            T,K,c,nt,nc=np.median([[r['T_s'],r['gain'],r['offset'],r['null_T_s'],r['null_offset']] for _,_,r in donors],axis=0)
            scores=[]
            for lag in np.arange(0,60.25,.5):
                terms=[]
                for z,local,r in donors:
                    t=z['time'];a,b=local['integral_probe']['training_bounds_s'];ix=(t>=a)&(t<=b);tt=t[ix];y=z[gas+'_raw'][ix]
                    pred=predict(tt,t,z['oxygen_cumulative_raw'],y[0],T,K,c,lag)
                    terms.append(float(np.mean((pred[1:]-y[1:])**2)/max(np.var(y[1:]),1e-12)))
                scores.append(np.mean(terms))
            lag=float(np.argmin(scores)*.5)
            z,local,_=objects[(held['heat_id'],gas,'weak')];t=z['time'];a,b=local['integral_probe']['validation_bounds_s'];ix=(t>=a)&(t<=b);tt=t[ix];y=z[gas+'_raw'][ix]
            pred=predict(tt,t,z['oxygen_cumulative_raw'],y[0],T,K,c,lag);null=predict(tt,t,z['oxygen_cumulative_raw'],y[0],nt,0,nc,0)
            rmse=float(np.sqrt(np.mean((pred[1:]-y[1:])**2)));nr=float(np.sqrt(np.mean((null[1:]-y[1:])**2)))
            loho.append(dict(heat_id=z['heat_id'],gas=gas,donor_heats=len(donors),delay_s=lag,T_s=float(T),gain=float(K),offset=float(c),
                             validation_rmse=rmse,null_validation_rmse=nr,validation_improvement=1-rmse/nr if nr>1e-8 else 0,
                             held_out_heat_parameter_fitting=False))
    write_csv(out/'real_trials.csv',rows);write_csv(out/'real_profiles.csv',profiles);write_csv(out/'real_bootstrap_draws.csv',draws);write_csv(out/'real_loho.csv',loho)
    import pandas as pd
    frame=pd.DataFrame(rows);loo=pd.DataFrame(loho);gates=[]
    for gas in ['co','co2']:
        q=frame[(frame.gas==gas)&(frame.method=='weak')];l=loo[loo.gas==gas]
        accepted=int(q.accepted_with_interval.sum());positive=int((l.validation_improvement>0).sum());med=float(q.validation_improvement.median())
        passed=accepted>=6 and positive>=8 and med>=.1
        gates.append(dict(gas=gas,accepted_heats=accepted,positive_loho_heats=positive,median_within_heat_heldout_improvement=med,
                          passed=passed,delay_range_s=[float(q.estimate_s.min()),float(q.estimate_s.max())],
                          T_range_s=[float(q.T_s.min()),float(q.T_s.max())]))
    write_json(out/'real_gate.json',gates)
    return dict(real_data_status='completed',real_validation=gates,
                flow_delay_analysis='not_supported' if not all(g['passed'] for g in gates) else 'eligible_for_separate_source_interval_design',
                adoption='do_not_adopt_as_physical_delay_model' if not all(g['passed'] for g in gates) else 'research_candidate_requires_independent_event_or_clock_validation')
