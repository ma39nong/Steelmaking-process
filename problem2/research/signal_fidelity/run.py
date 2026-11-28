"""Fixed benchmark protocol. Every output stays in this research directory."""
import argparse,csv,hashlib,importlib.metadata,json,platform,subprocess
from pathlib import Path
import numpy as np
from ...src.delay_audit import write_csv,write_json
from ...src.q2_pipeline import hampel_mask,oxygen_flow,_standardized_change
from .operators import acquire,source_methods,gas_methods,pulse,changes,fast_hampel
from .metrics import metrics,cumulative_metrics

ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent;OUT=HERE/'results'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_heats():
    manifest=json.loads((ROOT/'problem2/research/measurement_chain/results/input_manifest.json').read_text())
    for r in manifest:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    heats=[]
    for r in manifest:
        p=ROOT/r['path']
        if p.suffix!='.npz':continue
        with np.load(p,allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
        d['heat_id']=p.name.split('_')[0];heats.append(d)
    assert len(heats)==10
    return sorted(heats,key=lambda d:d['heat_id']),manifest


def verify_baseline(heats,cfg):
    rows=[]
    for z in heats:
        t=z['time'];src=source_methods(t,z['oxygen_cumulative_raw'],cfg['source_methods'])['formal']
        assert np.array_equal(src['mask'],hampel_mask(z['oxygen_cumulative_raw']))
        assert np.allclose(src['flow'],z['oxygen_flow'],rtol=1e-11,atol=1e-11)
        assert np.allclose(src['flow'],oxygen_flow(t,src['clean']),rtol=1e-11,atol=1e-11)
        for gas in ['co','co2']:
            raw=z[gas+'_raw'];formal=gas_methods(t,raw,cfg['gas_methods'])['formal']
            assert np.array_equal(formal['mask'],hampel_mask(raw))
            assert np.array_equal(formal['clean'],z[gas+'_clean'])
            assert np.array_equal(changes(src['flow'],'formal_median21'),_standardized_change(src['flow']))
            rows.append(dict(heat_id=z['heat_id'],gas=gas,source_match=True,gas_match=True,change_match=True))
    return rows


def source_case(t,base,center,width,amplitude,shape,mechanism,schedule,cfg,baseline=None,**labels):
    truth,added=pulse(t,center,width,amplitude,shape);volume=amplitude*width/(2 if shape=='triangle' else 1)
    ob=acquire(base,t,mechanism,schedule);oa=acquire(base+added,t,mechanism,schedule)
    baseline=source_methods(t,ob,cfg['source_methods']) if baseline is None else baseline
    aug=source_methods(t,oa,cfg['source_methods']);rows=[]
    common=dict(**labels,shape=shape,duration_s=width,amplitude=amplitude,acquisition=mechanism,event_center_s=center)
    for name in cfg['source_methods']:
        delta=aug[name]['flow']-baseline[name]['flow']
        rows.append(dict(**common,method=name,**metrics(t,delta,center,width,amplitude,shape,cfg),
                         changed_mask_samples=int(np.sum(aug[name]['mask']!=baseline[name]['mask'])),
                         sampled_total_baseline=float(np.trapezoid(baseline[name]['flow'],t)),
                         sampled_total_augmented=float(np.trapezoid(aug[name]['flow'],t))))
    stage=[]
    for name,delta in [('acquisition',oa-ob),('cleaning',aug['formal']['clean']-baseline['formal']['clean'])]:
        stage.append(dict(**common,stage=name,**cumulative_metrics(t,delta,added,center,width,volume)))
    change=[];truediff=np.diff(truth,prepend=truth[0]);den=float(np.sum(truediff**2))
    for name in cfg['change_methods']:
        delta=changes(aug['formal']['flow'],name)-changes(baseline['formal']['flow'],name)
        change.append(dict(**common,method=name,retained_energy_ratio=float(np.sum(delta**2)/den),
                           change_waveform_relative_l2_error=float(np.sqrt(np.sum((delta-truediff)**2)/den))))
    return rows,stage,change,dict(t=t,truth=truth,augmented=aug,baseline=baseline,observed=oa-ob)


def controlled(cfg,heat):
    t=np.arange(801)*cfg['dt_s'];source=[];stage=[];change=[];gasrows=[];examples=[]
    for shape in cfg['controlled_shapes']:
        for width in cfg['durations_s']:
            for base_flow in cfg['controlled_source_base_flows_m3_s']:
                for mechanism in cfg['source_acquisition']:
                    for phase in cfg['controlled_event_phases_s']:
                        rr,ss,cc,ex=source_case(t,100+base_flow*t,200+phase,width,5.,shape,mechanism,
                            heat['oxygen_cumulative_raw'][:len(t)],cfg,domain='controlled',heat_id='synthetic',seed=0,
                            control_base_flow=base_flow,control_event_phase_s=phase)
                        source.extend(rr);stage.extend(ss);change.extend(cc)
                        if phase==.15 and base_flow==5 and shape=='triangle' and width in [4,16] and mechanism in ['ideal','hold1_phase0']:
                            examples.append((f'source_{width}_{mechanism}',ex))
            for mechanism in cfg['gas_acquisition']:
                for gas in ['co','co2']:
                    for phase in cfg['controlled_event_phases_s']:
                        center=200+phase;baseline=np.full_like(t,10.);expected=pulse(t,center,width,1.,shape)[0]
                        ob=acquire(baseline,t,mechanism,heat[gas+'_raw'][:len(t)])
                        oa=acquire(baseline+expected,t,mechanism,heat[gas+'_raw'][:len(t)])
                        b=gas_methods(t,ob,cfg['gas_methods']);a=gas_methods(t,oa,cfg['gas_methods'])
                        gasrows.append(dict(domain='controlled',heat_id='synthetic',gas=gas,seed=0,shape=shape,duration_s=width,
                                            amplitude=1.,amplitude_multiplier=1.,acquisition=mechanism,method='acquisition_only',
                                            control_event_phase_s=phase,**metrics(t,oa-ob,center,width,1.,shape,cfg),changed_mask_samples=0))
                        for method in cfg['gas_methods']:
                            gasrows.append(dict(domain='controlled',heat_id='synthetic',gas=gas,seed=0,shape=shape,duration_s=width,
                                                amplitude=1.,amplitude_multiplier=1.,acquisition=mechanism,method=method,control_event_phase_s=phase,
                                                **metrics(t,a[method]['clean']-b[method]['clean'],center,width,1.,shape,cfg),
                                                changed_mask_samples=int(np.sum(a[method]['mask']!=b[method]['mask']))))
    return source,stage,change,gasrows,examples


def backgrounds(cfg,heats):
    source=[];stage=[];change=[];gasrows=[];remeasure=[]
    for hi,z in enumerate(heats):
        t=z['time'];o=z['oxygen_cumulative_raw'];n=len(t);eligible=np.flatnonzero((t>=90)&(t<=t[-1]-90))
        amp=cfg['source_amplitude_flow_median_multiplier']*float(np.median(z['oxygen_flow']))
        bases={}
        for mechanism in cfg['source_acquisition']:
            ob=acquire(o,t,mechanism,o);bases[mechanism]=source_methods(t,ob,cfg['source_methods'])
            remeasure.append(dict(heat_id=z['heat_id'],signal='oxygen',acquisition=mechanism,
                                 max_remeasurement_difference=float(max(abs(ob-o))),rmse_remeasurement=float(np.sqrt(np.mean((ob-o)**2)))))
        for seed in cfg['seeds']:
            rng=np.random.default_rng(seed+10000*hi);index=int(rng.choice(eligible));center=t[index]+rng.uniform(*cfg['event_subgrid_phase_range_s'])
            for width in cfg['durations_s']:
                for mechanism in cfg['source_acquisition']:
                    rr,ss,cc,_=source_case(t,o,center,width,amp,cfg['background_shape'],mechanism,o,cfg,bases[mechanism],
                                          domain='raw_background',heat_id=z['heat_id'],seed=seed)
                    source.extend(rr);stage.extend(ss);change.extend(cc)
        print('source fidelity',z['heat_id'],flush=True)
        for gi,gas in enumerate(['co','co2']):
            raw=z[gas+'_raw'];scale=float(np.std(raw));gb={}
            for mechanism in cfg['gas_acquisition']:
                ob=acquire(raw,t,mechanism,raw);gb[mechanism]=gas_methods(t,ob,cfg['gas_methods'])
                remeasure.append(dict(heat_id=z['heat_id'],signal=gas,acquisition=mechanism,
                                     max_remeasurement_difference=float(max(abs(ob-raw))),rmse_remeasurement=float(np.sqrt(np.mean((ob-raw)**2)))))
            for seed in cfg['seeds']:
                rng=np.random.default_rng(seed+10000*hi+1000*(gi+1));index=int(rng.choice(eligible));center=t[index]+rng.uniform(*cfg['event_subgrid_phase_range_s'])
                for width in cfg['durations_s']:
                    for multiplier in cfg['gas_amplitude_std_multipliers']:
                        amp=multiplier*scale;planted=pulse(t,center,width,amp,cfg['background_shape'])[0]
                        for mechanism in cfg['gas_acquisition']:
                            oa=acquire(raw+planted,t,mechanism,raw);a=gas_methods(t,oa,cfg['gas_methods']);b=gb[mechanism]
                            ob=acquire(raw,t,mechanism,raw)
                            gasrows.append(dict(domain='raw_background',heat_id=z['heat_id'],gas=gas,seed=seed,shape=cfg['background_shape'],
                                                duration_s=width,amplitude=amp,amplitude_multiplier=multiplier,acquisition=mechanism,
                                                method='acquisition_only',event_center_s=center,**metrics(t,oa-ob,center,width,amp,cfg['background_shape'],cfg),
                                                changed_mask_samples=0))
                            for method in cfg['gas_methods']:
                                gasrows.append(dict(domain='raw_background',heat_id=z['heat_id'],gas=gas,seed=seed,shape=cfg['background_shape'],
                                                    duration_s=width,amplitude=amp,amplitude_multiplier=multiplier,acquisition=mechanism,
                                                    method=method,event_center_s=center,**metrics(t,a[method]['clean']-b[method]['clean'],center,width,amp,cfg['background_shape'],cfg),
                                                    changed_mask_samples=int(np.sum(a[method]['mask']!=b[method]['mask']))))
            print('gas fidelity',z['heat_id'],gas,flush=True)
    return source,stage,change,gasrows,remeasure


def controls(cfg,heats):
    rows=[]
    for z in heats:
        t=z['time'];o=100+5*t
        for mechanism in cfg['source_acquisition']:
            # Transfer each heat's recorded update schedule to a known constant-flow process.
            ob=acquire(o,t,mechanism,z['oxygen_cumulative_raw']);b=source_methods(t,ob,cfg['source_methods'])
            corrupt=ob.copy();index=len(t)//2
            step=np.diff(z['oxygen_cumulative_raw']);spike=cfg['source_logging_spike_increment_multiplier']*float(np.median(step[step>0]))
            corrupt[index]+=spike;a=source_methods(t,corrupt,cfg['source_methods'])
            keep=(t>=10)&(t<=t[-1]-10)
            for method in cfg['source_methods']:
                q=b[method]['flow'];delta=a[method]['flow']-q
                rows.append(dict(heat_id=z['heat_id'],signal='oxygen',method=method,acquisition=mechanism,
                                 constant_flow_relative_rmse=float(np.sqrt(np.mean((q[keep]-5)**2))/5),
                                 constant_flow_ripple_cv=float(np.std(q[keep])/5),constant_flow_mean_ratio=float(np.mean(q[keep])/5),
                                 constant_flow_false_excursion_fraction=float(np.mean(abs(q[keep]-5)>=cfg['null_ripple_fraction']*5)),
                                 constant_flow_total_relative_error=float(np.trapezoid(q,t)/(5*(t[-1]-t[0]))-1),
                                 logging_spike_amplitude=spike,artifact_peak_response=float(max(abs(delta))),
                                 artifact_absolute_area=float(np.trapezoid(abs(delta),t)),
                                 corrupted_point_flagged=bool(a[method]['mask'][index])))
        for gas in ['co','co2']:
            base=z[gas+'_raw'];spike=cfg['gas_logging_spike_std_multiplier']*float(np.std(base))
            for mechanism in cfg['gas_acquisition']:
                ob=acquire(base,t,mechanism,base);corrupt=ob.copy();index=len(t)//2;corrupt[index]+=spike
                b=gas_methods(t,ob,cfg['gas_methods']);a=gas_methods(t,corrupt,cfg['gas_methods'])
                for method in cfg['gas_methods']:
                    delta=a[method]['clean']-b[method]['clean']
                    rows.append(dict(heat_id=z['heat_id'],signal=gas,method=method,acquisition=mechanism,
                                     logging_spike_amplitude=spike,artifact_peak_response=float(max(abs(delta))),
                                     artifact_absolute_area=float(np.trapezoid(abs(delta),t)),
                                     corrupted_point_flagged=bool(a[method]['mask'][index])))
    return rows


def integral(cfg):
    from .integral import execute
    rows,profiles=execute(cfg)
    write_csv(OUT/'integral_trials.csv',rows);write_csv(OUT/'integral_profiles.csv',profiles)
    write_json(OUT/'integral_config.json',cfg['integral_probe'])
    return rows,profiles


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--integral-only',action='store_true');args=parser.parse_args()
    cfg=json.loads((HERE/'config.json').read_text());OUT.mkdir(exist_ok=True)
    protected={p:sha(ROOT/p) for p in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
               if (ROOT/p).is_file() and not p.startswith('problem2/research/signal_fidelity/')}
    lock=json.loads((ROOT/'problem2/final/relationships/paper_freeze_manifest.json').read_text())
    for r in lock['files']:assert sha(ROOT/r['path'])==r['sha256']
    write_json(OUT/'effective_config.json',cfg)
    ir,ip=integral(cfg)
    if not args.integral_only:
        heats,manifest=load_heats();write_json(OUT/'input_manifest.json',manifest)
        write_csv(OUT/'baseline_replay.csv',verify_baseline(heats,cfg))
        cs,st,cc,cg,examples=controlled(cfg,heats[0])
        bs,bst,bc,bg,rm=backgrounds(cfg,heats)
        write_csv(OUT/'source_fidelity.csv',cs+bs);write_csv(OUT/'source_stages.csv',st+bst)
        write_csv(OUT/'change_fidelity.csv',cc+bc);write_csv(OUT/'gas_fidelity.csv',cg+bg)
        write_csv(OUT/'baseline_remeasurement.csv',rm);ctrl=controls(cfg,heats);write_csv(OUT/'negative_controls.csv',ctrl)
        oracle=[];tt=np.arange(801)*cfg['dt_s']
        for shape in cfg['controlled_shapes']:
            for width in cfg['durations_s']:
                q=pulse(tt,200.15,width,1.,shape)[0];diff=np.diff(q,prepend=q[0]);den=float(np.sum(diff**2))
                for method in cfg['change_methods']:
                    oracle.append(dict(shape=shape,duration_s=width,method=method,
                                       retained_energy_ratio=float(np.sum(changes(q,method)**2)/den)))
        write_csv(OUT/'oracle_change_fidelity.csv',oracle)
        from .present import deliver
        deliver(OUT,cfg,cs+bs,cc+bc,cg+bg,ctrl,ir,ip,examples,heats)
    for p,h in protected.items():assert sha(ROOT/p)==h,'protected file changed: '+p
    for r in lock['files']:assert sha(ROOT/r['path'])==r['sha256']
    write_json(OUT/'protection_check.json',dict(protected_tracked_files=len(protected),unchanged=True,frozen_files=len(lock['files']),frozen_unchanged=True))
    write_json(OUT/'provenance.json',dict(execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        python=platform.python_version(),packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','matplotlib']},
        config_sha256=sha(HERE/'config.json'),source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.glob('*.py'))},
        automatic_selection=False,phase='integral_only' if args.integral_only else 'complete_benchmark',
        outputs=[dict(path=str(p.relative_to(OUT)),sha256=sha(p)) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='provenance.json']))

if __name__=='__main__':main()
