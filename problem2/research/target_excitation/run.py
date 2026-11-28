"""Fixed-protocol A/B audit, separate from the frozen final deliverable."""
from __future__ import annotations
import copy,hashlib,json,platform,subprocess
from pathlib import Path
import numpy as np
from ...src.q2_pipeline import _standardized_change
from ...src.delay_v2 import local_delay_reliable
from ...src.delay_audit import write_json,write_csv
from .core import compare,packet,shift,trial_row

ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
PRIOR=ROOT/'problem2/final/delay_audit'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def baseline(out,cfg):
    rows=[];replay=[];inputs=[];heats=[]
    known={r['path']:r['sha256'] for r in json.loads((PRIOR/'provenance.json').read_text())['outputs']}
    for info in sorted(json.loads((PRIOR/'input_audit.json').read_text()),key=lambda d:d['heat_id']):
        heat=info['heat_id'];raw=ROOT/'Data'/Path(info['path']).name;p=PRIOR/'baseline/data'/f'{heat}_q2_processed.npz'
        assert sha(raw)==info['sha256'] and sha(p)==known['baseline/data/'+p.name]
        inputs.extend(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in [raw,p])
        with np.load(p,allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
        d['heat_id']=heat;heats.append(d)
        for gas in ['co','co2']:
            pair=compare(d['oxygen_flow'],d[gas+'_clean'],cfg);a=pair['A']
            formal=local_delay_reliable(d['oxygen_flow'],d[gas+'_clean'],cfg['dt_s'])
            expected=dict(center_index='center_index',center_time='center_time',observed='delay_observed',accepted='reliable',score='max_corr',
                separation='peak_separation',source_std='source_std',unshifted_target_std='target_std',boundary='boundary_hit',confidence='confidence',prominence='peak_prominence',width='peak_width_samples')
            for key,field in expected.items():
                exact=np.array_equal(a[key],formal[field],equal_nan=True)
                close=np.allclose(a[key],formal[field],rtol=1e-9,atol=1e-10,equal_nan=True)
                assert close,(heat,gas,key)
                replay.append(dict(heat_id=heat,gas=gas,field=field,exact_match=exact,numerical_match=close))
            assert np.array_equal(a['accepted'],d[f'delay_{gas}_valid_mask'])
            assert np.allclose(a['observed'],d[f'delay_{gas}_observed'],equal_nan=True)
            for method in ['A','B']:
                for i in range(len(a['center_index'])):rows.append(trial_row(pair,method,i,None,cfg,heat_id=heat,gas=gas,scenario='original_real'))
    print('formal A replay verified for all ten heats',flush=True)
    write_csv(out/'baseline_replay.csv',replay);write_csv(out/'real_windows.csv',rows);write_json(out/'input_manifest.json',inputs)
    return heats,rows


def synthetic(out,cfg):
    sc=cfg['synthetic'];n=sc['samples'];dt=cfg['dt_s'];center=int(sc['evaluation_center_s']/dt);rng0=np.random.default_rng(0);rows=[];examples=[]
    for seed in cfg['seeds']:
        rng=np.random.default_rng(seed)
        for offset in sc['source_offsets_s']:
            template=packet(n,center+int(offset/dt),rng,dt)
            for noise in sc['noise_ratios']:
                for lag in sc['delays_s']:
                    source=1000+template+rng.normal(size=n)*noise*np.std(template)
                    target=1000+shift(template,lag,dt)+rng.normal(size=n)*noise*np.std(template)
                    pair=compare(source,target,cfg,[center]);scenario='centered_response' if offset==0 else 'late_response'
                    for method in ['A','B']:rows.append(trial_row(pair,method,0,lag,cfg,scenario=scenario,seed=seed,noise_ratio=noise,delay_s=lag))
                    if seed==cfg['seeds'][0] and noise==0 and lag==30 and offset==20:
                        examples.append(('late_response',source,target,pair,center,lag))
        template=packet(n,center-int(20/dt),rng,dt);other=packet(n,center-int(25/dt),rng,dt)
        source=1000+template;target=1000+40*other+sc['weak_copy_gain']*shift(template,sc['weak_copy_delay_s'],dt)
        pair=compare(source,target,cfg,[center])
        for method in ['A','B']:rows.append(trial_row(pair,method,0,sc['weak_copy_delay_s'],cfg,scenario='weak_real_copy',seed=seed,noise_ratio=0.,delay_s=sc['weak_copy_delay_s']))
        if seed==cfg['seeds'][0]:examples.append(('weak_copy',source,target,pair,center,sc['weak_copy_delay_s']))
        for scenario in ['stationary_unrelated','quiet_future_unrelated']:
            template=packet(n,center+int(20/dt),rng,dt);source=1000+template
            noise=np.convolve(rng.normal(size=n),np.ones(5)/5,mode='same')
            target=1000+noise if scenario=='stationary_unrelated' else 1000+40*other+.001*noise
            pair=compare(source,target,cfg,[center])
            for method in ['A','B']:rows.append(trial_row(pair,method,0,None,cfg,scenario=scenario,seed=seed,noise_ratio=1.,delay_s=None))
        for lag in sc['boundary_delays_s']:
            template=packet(n,center,rng,dt);pair=compare(1000+template,1000+shift(template,lag,dt),cfg,[center])
            for method in ['A','B']:rows.append(trial_row(pair,method,0,lag,cfg,scenario='boundary_reference',seed=seed,noise_ratio=0.,delay_s=lag))
    write_csv(out/'synthetic_trials.csv',rows)
    return rows,examples


def injection(out,cfg,heats):
    rows=[];ic=cfg['injection'];dt=cfg['dt_s'];half=cfg['window']//2
    for heat_index,d in enumerate(heats):
        flow=d['oxygen_flow'];n=len(flow);centers=d['delay_co_center_index'];base_source_scale=np.std(_standardized_change(flow))
        for gas_index,gas in enumerate(['co','co2']):
            original=d[gas+'_clean'];base_target_scale=np.std(_standardized_change(original))
            for seed in cfg['seeds']:
                rng=np.random.default_rng(seed+10000*heat_index+1000*gas_index);center=int(rng.choice(centers))
                template=packet(n,center+int(ic['source_offset_from_window_center_s']/dt),rng,dt)
                template_scale=np.std(_standardized_change(template)[center-half:center+half+1])
                aq=ic['source_change_snr']*base_source_scale/template_scale
                source=flow+aq*template
                null=compare(source,original,cfg,[center])
                for method in ['A','B']:rows.append(trial_row(null,method,0,None,cfg,scenario='real_source_only',heat_id=d['heat_id'],gas=gas,seed=seed,target_snr=None,delay_s=None,source_amplitude=aq,target_amplitude=0.))
                for gain in ic['target_change_snr']:
                    ag=gain*base_target_scale/template_scale
                    for lag in ic['delays_s']:
                        target=original+ag*shift(template,lag,dt);pair=compare(source,target,cfg,[center])
                        for method in ['A','B']:rows.append(trial_row(pair,method,0,lag,cfg,scenario='real_paired_injection',heat_id=d['heat_id'],gas=gas,seed=seed,target_snr=gain,delay_s=lag,source_amplitude=aq,target_amplitude=ag))
            print('injection complete',d['heat_id'],gas,flush=True)
    write_csv(out/'injection_trials.csv',rows)
    return rows


def summarize(rows,groupfields):
    groups={}
    for r in rows:
        key=tuple(r.get(k) for k in groupfields)+ (r['method'],)
        groups.setdefault(key,[]).append(r)
    result=[]
    for key,rr in groups.items():
        count=len(rr);accepted=[r for r in rr if r['accepted']];errors=[r['absolute_error_s'] for r in accepted if r['absolute_error_s'] is not None]
        true=[r for r in rr if r['truth_s'] is not None];correct=sum(r['correct'] for r in true)
        result.append(dict(zip(groupfields+['method'],key),trials=count,accepted=len(accepted),acceptance_rate=len(accepted)/count,rejection_rate=1-len(accepted)/count,
            known_delay_trials=len(true),correct_recoveries=correct,correct_recovery_rate=correct/len(true) if true else None,
            accepted_wrong_delay=sum(r['false_delay'] for r in true),wrong_delay_rate_all=sum(r['false_delay'] for r in true)/len(true) if true else None,
            conditional_mae_s=float(np.mean(errors)) if errors else None,conditional_p90_error_s=float(np.quantile(errors,.9)) if errors else None,conditional_p99_error_s=float(np.quantile(errors,.99)) if errors else None,conditional_max_error_s=float(max(errors)) if errors else None,
            B_full_competitor_acceptances=sum(r['full_competitor_accept'] for r in rr),guarded_correct_recoveries=sum(r['full_competitor_accept'] and r['correct'] for r in rr),guarded_wrong_delays=sum(r['full_competitor_accept'] and r['false_delay'] for r in rr),unshifted_target_gate_failures=sum(not r['unshifted_target_gate'] and r['source_gate'] for r in rr),
            accepted_selected_target_below_threshold=sum(r['accepted'] and r['selected_target_ratio']<.05 for r in rr),
            accepted_due_to_masked_competitors=sum(r['accepted'] and not r['full_competitor_accept'] for r in rr)))
    return result


def main():
    cfg=json.loads((HERE/'config.json').read_text());out=HERE/'results';out.mkdir(exist_ok=True)
    lock=ROOT/'problem2/final/relationships/paper_freeze_manifest.json';frozen=json.loads(lock.read_text())
    for r in frozen['files']:assert sha(ROOT/r['path'])==r['sha256'],'frozen deliverable changed: '+r['path']
    write_json(out/'effective_config.json',cfg);write_json(out/'protocol.json',dict(research_hypothesis='unshifted target gate can reject a later response or admit a weak shifted fluctuation',
        only_primary_change='target excitation location; all correlations, centers, thresholds and search range shared',
        locked_config_sha256=sha(HERE/'config.json'),formal_freeze_manifest_sha256=sha(lock),
        real_injection_unit='one designated center per heat-gas-seed; repeats across doses are paired, not independent industrial events',
        safety_diagnostic='B separation also checked against all raw competitors; not proposed as a replacement model',
        limits='No physical event truth on original data; real-source-only acceptance is no-planted-response acceptance, not an industrial false-positive rate'))
    heats,real=baseline(out,cfg);syn,examples=synthetic(out,cfg);inj=injection(out,cfg,heats)
    write_csv(out/'real_summary.csv',summarize(real,['gas']))
    write_csv(out/'synthetic_summary.csv',summarize(syn,['scenario','noise_ratio','delay_s']))
    write_csv(out/'injection_summary.csv',summarize(inj,['scenario','gas','target_snr','delay_s']))
    write_csv(out/'scenario_summary.csv',summarize(syn+inj,['scenario']))
    from .present import deliver
    deliver(out,cfg,examples,real,syn,inj)
    for r in frozen['files']:assert sha(ROOT/r['path'])==r['sha256'],'frozen deliverable changed during experiment'
    write_json(out/'provenance.json',dict(formal_baseline=cfg['formal_baseline_commit'],execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        python=platform.python_version(),source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.glob('*.py'))},
        config_sha256=sha(HERE/'config.json'),outputs=[dict(path=str(p.relative_to(out)),sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='provenance.json']))

if __name__=='__main__':main()
