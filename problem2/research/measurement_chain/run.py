"""Replay the raw chain and measure timing, without changing formal outputs."""
import csv, hashlib, json, math, platform, subprocess
from collections import Counter
from pathlib import Path
import numpy as np
from pyxlsb import open_workbook
from scipy.signal import savgol_filter
from ...src.delay_audit import write_csv, write_json
from ...src.q2_pipeline import hampel_mask, pchip_clean, oxygen_flow
from .core import integrate, reconstruct, pulse, held, errors, runs, lag_row,delay_features

ROOT=Path(__file__).resolve().parents[3]; HERE=Path(__file__).resolve().parent
BASE=ROOT/'problem2/final/delay_audit'


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def raw_audit(out,cfg):
    inventory=json.loads((BASE/'input_audit.json').read_text())
    manifest=[]; headers=[]; stats=[]; intervals=[]; plateaus=[]; increments=[]
    sync=[]; replay=[]; flows=[]; repairs=[]; heats=[]; periodic=[]; grid=[]
    for info in sorted(inventory,key=lambda r:r['heat_id']):
        h=info['heat_id']; p=ROOT/'Data'/Path(info['path']).name
        assert sha(p)==info['sha256']
        bp=BASE/'baseline/data'/f'{h}_q2_processed.npz'
        with np.load(bp,allow_pickle=False) as zz: z={k:zz[k].copy() for k in zz.files}
        z['heat_id']=h; heats.append(z)
        manifest.extend([dict(path=str(p.relative_to(ROOT)),sha256=sha(p)),dict(path=str(bp.relative_to(ROOT)),sha256=sha(bp))])
        # Actual XLSB rows and all headers are examined, not just cached arrays.
        with open_workbook(p) as wb:
            assert len(wb.sheets)==1
            sheet_name=wb.sheets[0]
            with wb.get_sheet(1) as sheet:
                rr=sheet.rows(); header=[c.v for c in next(rr)]
                cols={name:header.index(info['columns'][name]) for name in ['oxygen_cumulative','co','co2']}
                raw=np.array([[r[0].v]+[r[cols[k]].v for k in cols] for r in rr],float)
        assert len(header)==2051 and np.array_equal(raw[:,0],z['sequence'])
        headers.append(dict(heat_id=h,sheet=sheet_name,column_count=len(header),sequence_header=str(header[0]),
                            process_headers={k:str(header[v]) for k,v in cols.items()},
                            other_columns_are_numeric_wavelengths=all(isinstance(v,(float,int)) for v in header[1:2048]),
                            independent_time_or_operation_columns=[]))
        t=z['time']; dt=cfg['dt_s']; active=np.flatnonzero(np.diff(z['oxygen_cumulative_raw'])>0)+1
        active_slice=slice(int(active[0]),int(active[-1])+1)
        changes={}
        for j,k in enumerate(cols,1):
            x=z[k+'_raw']; assert np.array_equal(raw[:,j],x,equal_nan=True)
            mask=hampel_mask(x); clean=pchip_clean(t,x,mask)
            assert np.array_equal(mask,z['outlier_oxygen' if k=='oxygen_cumulative' else 'outlier_'+k])
            assert np.allclose(clean,z[k+'_clean'],rtol=1e-12,atol=1e-12)
            replay.append(dict(heat_id=h,signal=k,raw_exact=True,mask_exact=True,clean_match=True))
            for a,b in runs(mask.astype(int)):
                if mask[a]:
                    left=max(0,a-1);right=min(len(x)-1,b)
                    repairs.append(dict(heat_id=h,signal=k,start_s=float(t[a]),end_s=float(t[b-1]),
                                        removed_samples=b-a,visible_support_span_s=float(t[right]-t[left]),edge=a==0 or b==len(x)))
            d=np.diff(x); ix=np.flatnonzero(np.isfinite(d)&(d!=0))+1; changes[k]=np.isfinite(d)&(d!=0)
            pos=d[d>0]
            for scope,xx in [('whole',x),('oxygen_active',x[active_slice])]:
                dd=np.diff(xx); ii=np.flatnonzero(np.isfinite(dd)&(dd!=0))+1
                spacing=np.diff(ii)*dt
                mode=Counter(spacing).most_common(1)
                stats.append(dict(heat_id=h,signal=k,scope=scope,samples=len(xx),adjacent_pairs=len(dd),
                                  repeated_pairs=int(sum(dd==0)),repeat_fraction=float(np.mean(dd==0)),
                                  positive_updates=int(sum(dd>0)),negative_updates=int(sum(dd<0)),
                                  update_interval_mode_s=float(mode[0][0]) if mode else None,
                                  mode_share=mode[0][1]/len(spacing) if mode else None,
                                  update_interval_p50_s=float(np.median(spacing)) if len(spacing) else None,
                                  update_interval_p95_s=float(np.quantile(spacing,.95)) if len(spacing) else None,
                                  positive_increment_median=float(np.median(dd[dd>0])) if np.any(dd>0) else None))
                for gap,count in sorted(Counter(spacing).items()):
                    intervals.append(dict(heat_id=h,signal=k,scope=scope,interval_s=float(gap),count=count))
            for a,b in runs(x):
                plateaus.append(dict(heat_id=h,signal=k,start_index=a,end_exclusive=b,samples=b-a,
                                     equal_transition_duration_s=(b-a-1)*dt,grid_occupancy_s=(b-a)*dt,
                                     edge_censored=a==0 or b==len(x),value=float(x[a])))
            for i in ix:
                increments.append(dict(heat_id=h,signal=k,index=int(i),time_s=float(t[i]),increment=float(d[i-1])))
            indicator=changes[k].astype(float)
            for lag in range(1,41):
                periodic.append(dict(heat_id=h,signal=k,lag_s=lag*dt,
                                     update_indicator_correlation=float(np.corrcoef(indicator[:-lag],indicator[lag:])[0,1])))
            for quantum in [1.,.1,.01,.001,.0001,.000001]:
                residual=abs(x/quantum-np.rint(x/quantum))
                grid.append(dict(heat_id=h,signal=k,tested_grid=quantum,
                                 fraction_within_1e_6_grid_unit=float(np.mean(residual<1e-6)),
                                 meaning='representation test, not device quantization identification'))
        for a,b in [('oxygen_cumulative','co'),('oxygen_cumulative','co2'),('co','co2')]:
            aa,bb=changes[a],changes[b];both=int(sum(aa&bb))
            sync.append(dict(heat_id=h,a=a,b=b,changes_a=int(sum(aa)),changes_b=int(sum(bb)),
                             simultaneous_updates=both,jaccard=both/int(sum(aa|bb)),
                             observed_joint_fraction=float(np.mean(aa&bb)),
                             product_of_marginal_fractions=float(np.mean(aa)*np.mean(bb))))
        q=oxygen_flow(t,z['oxygen_cumulative_clean']); assert np.allclose(q,z['oxygen_flow'],rtol=1e-12,atol=1e-12)
        oc=z['oxygen_cumulative_clean']; area=integrate(t,q);span=oc[-1]-oc[0]; res=oc[0]+area-oc
        raw_der=savgol_filter(oc,31,3,deriv=1,delta=dt,mode='interp')
        flows.append(dict(heat_id=h,oxygen_span_m3=float(span),integral_m3=float(area[-1]),
                          final_relative_error=float((area[-1]-span)/span),max_path_residual_m3=float(max(abs(res))),
                          path_rmse_m3=float(np.sqrt(np.mean(res**2))),negative_unclipped_derivative_fraction=float(np.mean(raw_der<0)),
                          zero_clipping_added_m3=float(integrate(t,q-raw_der)[-1]),
                          first_positive_increment_s=float(t[active[0]]),last_positive_increment_s=float(t[active[-1]]),
                          first_positive_is_record_start=bool(active[0]==1),last_positive_is_record_end=bool(active[-1]==len(t)-1)))
        print('raw and reconstruction verified',h,flush=True)
    for name,rr in [('raw_signal_summary',stats),('update_intervals',intervals),('constant_runs',plateaus),
                    ('raw_increments',increments),('update_coincidence',sync),('reconstruction_replay',replay),
                    ('flow_integral',flows),('repair_spans',repairs),('update_autocorrelation',periodic),('representation_grids',grid)]:
        write_csv(out/(name+'.csv'),rr)
    write_json(out/'raw_headers.json',headers);write_json(out/'input_manifest.json',manifest)
    return heats,stats,flows


def model_chain(out,cfg):
    sc=cfg['simulation'];dt=cfg['dt_s'];t=np.arange(sc['samples'])*dt;center=sc['source_center_s'];rows=[];lags=[];examples=[]
    for shape in sc['shapes']:
        for width in sc['durations_s']:
            truth=sc['amplitude_m3_s']*pulse(t,center,width,shape); oxygen=integrate(t,truth)
            # Gas processing is cached by acquisition period/phase/delay and shared.
            gases=[]
            for gp in sc['gas_update_periods_s']:
                for phase in range(round(gp/dt)):
                    for delay in sc['delays_s']:
                        y=sc['amplitude_m3_s']*pulse(t,center+delay,width,shape)
                        observed=held(y,round(gp/dt),phase);gm=hampel_mask(observed)
                        yc=pchip_clean(t,observed,gm)
                        gases.append((gp,phase,delay,yc,errors(t,y,yc,'gas_')))
            for period in sc['counter_update_periods_s']:
                for phase in range(round(period/dt)):
                    for quantum in sc['counter_quanta_m3']:
                        measured=held(oxygen,round(period/dt),phase)
                        if quantum:measured=np.rint(measured/quantum)*quantum
                        q,oc,mask=reconstruct(t,measured)
                        label=dict(shape=shape,duration_s=width,counter_period_s=period,counter_phase_s=phase*dt,quantum_m3=quantum)
                        rows.append(dict(**label,**errors(t,truth,q),flagged_oxygen_points=int(mask.sum()),
                                         SG_only_vs_full_max_difference_m3_s=float(np.max(abs(oxygen_flow(t,measured)-q))),
                                         final_integral_error_relative=float((integrate(t,q)[-1]-oxygen[-1])/oxygen[-1])))
                        for gp,gphase,delay,yc,ge in gases:
                            lags.append(lag_row(q,yc,round(center/dt),delay,cfg,**label,
                                                gas_period_s=gp,gas_phase_s=gphase*dt,**ge,
                                                **delay_features(t,q,yc,delay,shape=='triangle')))
                        if shape=='triangle' and width in [2,15] and period in [.5,2] and phase==0 and quantum==0:
                            examples.append(dict(t=t,truth=truth,measured=measured,flow=q,duration=width,period=period))
            print('controlled acquisition chain',shape,width,flush=True)
    write_csv(out/'controlled_reconstruction.csv',rows);write_csv(out/'controlled_delay.csv',lags)
    return rows,lags,examples


def background_chain(out,cfg,heats):
    rows=[];timing=[];dt=cfg['dt_s']
    for hi,z in enumerate(heats):
        t=z['time'];o=z['oxygen_cumulative_raw'];base=z['oxygen_flow'];amp=float(np.median(base))
        for gi,gas in enumerate(['co','co2']):
            y=z[gas+'_raw'];ag=float(np.std(y));baseline_gas=z[gas+'_clean']
            for seed in cfg['seeds']:
                rng=np.random.default_rng(seed+10000*hi+1000*gi)
                center=int(rng.choice(z['delay_co_center_index']));ct=t[center]
                for width in cfg['background']['durations_s']:
                    planted=amp*pulse(t,ct,width)
                    # Full chain starts at raw cumulative O, rather than estimated q.
                    q,oc,mask=reconstruct(t,o+integrate(t,planted));direct=base+planted
                    labels=dict(heat_id=z['heat_id'],gas=gas,seed=seed,center_index=center,duration_s=width,
                                source_amplitude_m3_s=amp,target_amplitude_m3_s=ag)
                    timing.append(dict(**labels,**errors(t,planted,q-base),
                                       changed_oxygen_mask_points=int(np.sum(mask!=z['outlier_oxygen']))))
                    for method,source in [('direct_flow',direct),('raw_counter_chain',q)]:
                        rows.append(lag_row(source,baseline_gas,center,None,cfg,**labels,method=method,scenario='source_only'))
                    for delay in cfg['background']['delays_s']:
                        yy=y+ag*pulse(t,ct+delay,width);gm=hampel_mask(yy);yc=pchip_clean(t,yy,gm)
                        for method,source in [('direct_flow',direct),('raw_counter_chain',q)]:
                            rows.append(lag_row(source,yc,center,delay,cfg,**labels,method=method,scenario='known_addition'))
                            rows[-1].update(delay_features(t,source-base,yc-baseline_gas,delay))
            print('raw-background chain',z['heat_id'],gas,flush=True)
    write_csv(out/'background_trials.csv',rows);write_csv(out/'background_source_timing.csv',timing)
    return rows,timing


def prior_cluster(out,cfg):
    prior=ROOT/'problem2/research/target_excitation/results/injection_by_heat.csv'
    with prior.open() as f:rr=list(csv.DictReader(f))
    heats=sorted({r['heat_id'] for r in rr});ar=[]
    for heat in heats:
        sel=[r for r in rr if r['heat_id']==heat];n=sum(int(r['paired_conditions']) for r in sel)
        a=sum(int(r['A_correct']) for r in sel);b=sum(int(r['B_correct']) for r in sel)
        ar.append(dict(heat_id=heat,n=n,A_correct=a,B_correct=b,net_correct=b-a,A_rate=a/n,B_rate=b/n,gain_rate=(b-a)/n))
    write_csv(out/'prior_AB_by_heat.csv',ar)
    totals=np.array([[r['n'],r['A_correct'],r['B_correct']] for r in ar]);total=totals.sum(axis=0)
    loo=[]
    for i,h in enumerate(heats):
        n,a,b=total-totals[i];loo.append(dict(omitted_heat=h,n=int(n),A_rate=a/n,B_rate=b/n,gain_rate=(b-a)/n))
    write_csv(out/'prior_AB_leave_one_heat_out.csv',loo)
    rng=np.random.default_rng(cfg['cluster_bootstrap']['seed']);s=totals[rng.integers(0,len(heats),(cfg['cluster_bootstrap']['replicates'],len(heats)))].sum(axis=1)
    delta=(s[:,2]-s[:,1])/s[:,0]; ci=np.quantile(delta,[.025,.975])
    result=dict(paired_unit='whole heat; all gas/seed/12 conditions retained together',clusters=len(heats),
                conditions=int(total[0]),distinct_heat_gas_seed_draws=600,
                caveat='600 background draws need not be 600 distinct time indices; centers are sampled with replacement',
                gain_rate=float((total[2]-total[1])/total[0]),percentile_95_lower=float(ci[0]),percentile_95_upper=float(ci[1]),
                interpretation='descriptive heat-cluster bootstrap of this injection design; not a physical-event confidence interval',
                input=dict(path=str(prior.relative_to(ROOT)),sha256=sha(prior)))
    write_json(out/'prior_AB_cluster_bootstrap.json',result)
    return ar,loo,result


def cached_numeric(out,cfg):
    """Rerender verified numeric CSVs without repeating raw acquisition experiments."""
    old=json.loads((out/'provenance.json').read_text())
    assert json.loads((out/'effective_config.json').read_text())==cfg
    for r in old['outputs']:assert sha(out/r['path'])==r['sha256'],'cached output modified: '+r['path']
    for r in json.loads((out/'input_manifest.json').read_text()):assert sha(ROOT/r['path'])==r['sha256']
    def read(name):
        def value(s):
            if s=='':return None
            if s in ['True','False']:return s=='True'
            try:return float(s)
            except ValueError:return s
        with (out/(name+'.csv')).open() as f:return [{k:value(v) for k,v in r.items()} for r in csv.DictReader(f)]
    heats=[]
    for p in sorted((BASE/'baseline/data').glob('*_q2_processed.npz')):
        with np.load(p,allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
        d['heat_id']=p.name.split('_')[0];heats.append(d)
    sc=cfg['simulation'];t=np.arange(sc['samples'])*cfg['dt_s'];examples=[]
    for width in [2.,15.]:
        for period in [.5,2.]:
            truth=sc['amplitude_m3_s']*pulse(t,sc['source_center_s'],width)
            measured=held(integrate(t,truth),round(period/cfg['dt_s']))
            q,_,_=reconstruct(t,measured)
            examples.append(dict(t=t,truth=truth,measured=measured,flow=q,duration=width,period=period))
    return (heats,read('raw_signal_summary'),read('flow_integral'),read('controlled_reconstruction'),
            read('controlled_delay'),examples,read('background_trials'),read('background_source_timing')),old


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--render-only',action='store_true');args=parser.parse_args()
    cfg=json.loads((HERE/'config.json').read_text());out=HERE/'results';out.mkdir(exist_ok=True)
    protected={p:sha(ROOT/p) for p in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
               if (ROOT/p).is_file() and not p.startswith('problem2/research/measurement_chain/')}
    lock=ROOT/'problem2/final/relationships/paper_freeze_manifest.json';frozen=json.loads(lock.read_text())
    for r in frozen['files']:assert sha(ROOT/r['path'])==r['sha256']
    old=None
    if args.render_only:
        (heats,stats,flows,model,lags,examples,background,timing),old=cached_numeric(out,cfg)
    else:
        write_json(out/'effective_config.json',cfg)
        heats,stats,flows=raw_audit(out,cfg)
        model,lags,examples=model_chain(out,cfg)
        background,timing=background_chain(out,cfg,heats)
    ar,loo,boot=prior_cluster(out,cfg)
    from .present import deliver
    deliver(out,cfg,heats,stats,flows,model,lags,examples,background,timing,ar,loo,boot)
    for p,digest in protected.items():assert sha(ROOT/p)==digest,'protected tracked file changed: '+p
    for r in frozen['files']:assert sha(ROOT/r['path'])==r['sha256']
    write_json(out/'protection_check.json',dict(protected_tracked_files=len(protected),unchanged=True,frozen_files=len(frozen['files']),frozen_unchanged=True))
    write_json(out/'provenance.json',dict(execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        python=platform.python_version(),formal_baseline=cfg['formal_baseline'],research_baseline=cfg['research_baseline'],
        numeric_execution_source_hashes=(old.get('numeric_execution_source_hashes') or old['source_hashes']) if old else None,
        configuration_sha256=sha(HERE/'config.json'),source_hashes={str(p.relative_to(ROOT)):sha(p) for p in HERE.glob('*.py')},
        outputs=[dict(path=str(p.relative_to(out)),sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='provenance.json']))

if __name__=='__main__':main()
