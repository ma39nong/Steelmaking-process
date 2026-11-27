"""Reproduce and audit Q2 delays without replacing the formal model.

Run: python -m problem2.src.delay_audit --data-dir Data
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from .delay_audit_core import scan, group_windows, error_metrics
from .q2_pipeline import hampel_mask, pchip_clean, oxygen_flow, _standardized_change
from .q2_pipeline_v2 import process_one_v2
from .delay_v2 import fit_delay_models
from .data_validation import ALIASES, _find_column
from pyxlsb import open_workbook

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serial(value):
    if isinstance(value, dict): return {k:serial(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [serial(v) for v in value]
    if isinstance(value, np.ndarray): return serial(value.tolist())
    if isinstance(value, np.generic): return serial(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(serial(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def write_csv(path, rows):
    if not rows: return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=fields, lineterminator="\n"); writer.writeheader()
        writer.writerows(serial(rows))


def compare_tree(actual, expected, path=''):
    failures=[]
    if isinstance(expected, dict):
        if set(actual) != set(expected): failures.append(path+': keys')
        for k in expected.keys() & actual.keys(): failures += compare_tree(actual[k],expected[k],path+'/'+k)
    elif isinstance(expected, list):
        if len(actual)!=len(expected): failures.append(path+': length')
        else:
            for i,(a,b) in enumerate(zip(actual,expected)): failures += compare_tree(a,b,path+'/'+str(i))
    elif isinstance(expected, (int,float)) and not isinstance(expected,bool):
        if not np.isclose(actual,expected,rtol=1e-7,atol=1e-8,equal_nan=True): failures.append(path)
    elif actual != expected: failures.append(path)
    return failures


def settings(cfg, window=None):
    p=cfg['delay']; dt=cfg['time_step_s']
    return dict(dt=dt,window=window or p['window'],max_lag=int(round(p['max_lag_s']/dt)),
                stride=p['stride'], min_corr=p['min_corr'],min_separation=p['min_separation'],
                min_strength_ratio=p['min_strength_ratio'])


def global_scan(q,y,dt,maxlag):
    # Pairwise normalized correlations on detrended changes; lag>0 means target follows source.
    a,b=_standardized_change(q),_standardized_change(y)
    lags=np.arange(-maxlag,maxlag+1); scores=[]
    for lag in lags:
        aa,bb=(a[:len(a)-lag],b[lag:]) if lag>=0 else (a[-lag:],b[:len(b)+lag])
        scores.append(np.corrcoef(aa,bb)[0,1] if min(aa.std(),bb.std())>0 else np.nan)
    best=int(np.nanargmax(scores)) if np.isfinite(scores).any() else None
    return dict(best_delay_s=float(lags[best]*dt) if best is not None else None,
                max_corr=float(scores[best]) if best is not None else None,
                boundary_hit=best in (0,len(lags)-1),
                lags_s=lags*dt, correlations=scores,
                interpretation='descriptive global peak; no event identification or significance claim')


def real_audit(data_dir,out,cfg,audit_cfg):
    with (ROOT/'problem2/audit/data_manifest.csv').open() as f:
        manifest={r['heat_id']:r for r in csv.DictReader(f)}
    saved={r['heat_id']:r for r in json.loads((ROOT/'problem2/final/results/q2_summary.json').read_text())}
    saved_validation={r['heat_id']:r for r in json.loads((ROOT/'problem2/final/results/process_validation.json').read_text())}
    files=sorted(data_dir.glob('*.xlsb*'))
    if len(files)!=10: raise ValueError('audit requires exactly ten XLSB files')
    baselines=out/'baseline'; rows=[]; windows=[]; accepted=[]; sensitivity=[]; comparisons=[]; sources=[]; globals_=[]; contexts={}
    for path in files:
        heat=next((h for h in manifest if h in path.name),None)
        if heat is None: raise ValueError('unknown heat: '+path.name)
        sha=digest(path)
        if sha!=manifest[heat]['sha256']: raise ValueError('input hash mismatch: '+heat)
        with open_workbook(str(path)) as book:
            if book.sheets != [saved_validation[heat]['sheet_name']]: raise AssertionError('sheet mismatch')
            with book.get_sheet(1) as sheet:
                header=[cell.v for cell in next(sheet.rows())]
        columns={k:str(_find_column(header,v)) for k,v in ALIASES.items()}
        if len(header)!=2051 or columns!=saved_validation[heat]['process_columns']:
            raise AssertionError('actual process header differs from saved validation')
        summary=process_one_v2(path, baselines/'data', baselines/'results')
        failures=compare_tree(summary,saved[heat])
        if failures: raise AssertionError('saved baseline summary mismatch: '+repr(failures))
        with np.load(baselines/'data'/f'{heat}_q2_processed.npz',allow_pickle=False) as z:
            data={k:z[k].copy() for k in z.files}
        t=data['time']; dt=cfg['time_step_s']; seq=data['sequence']
        if not np.array_equal(seq,np.arange(1,len(t)+1)) or not np.array_equal(t,(seq-1)*dt): raise AssertionError('time mismatch')
        if len(t)!=saved_validation[heat]['records']: raise AssertionError('record count mismatch')
        clean={}; array_checks={}
        for key,maskkey in [('oxygen_cumulative','oxygen'),('co','co'),('co2','co2')]:
            raw=data[key+'_raw']; mask=hampel_mask(raw,cfg['hampel_window'],cfg['hampel_threshold'])
            clean[key]=pchip_clean(t,raw,mask)
            array_checks['outlier_'+maskkey]=bool(np.array_equal(mask,data['outlier_'+maskkey]))
            array_checks[key+'_clean']=bool(np.array_equal(clean[key],data[key+'_clean']))
        q=oxygen_flow(t,clean['oxygen_cumulative'],cfg['flow']['window'],cfg['flow']['polyorder'])
        array_checks['oxygen_flow']=bool(np.array_equal(q,data['oxygen_flow']))
        if not all(array_checks.values()): raise AssertionError('configured preprocessing differs from frozen baseline')
        total=clean['oxygen_cumulative'][-1]-clean['oxygen_cumulative'][0]
        sources.append(dict(heat_id=heat,path=str(path.resolve()),sha256=sha,manifest_match=True,records=len(t),
                            missing={k:int((~np.isfinite(data[k+'_raw'])).sum()) for k in clean},
                            negative={k:int((data[k+'_raw']<0).sum()) for k in clean},
                            columns=columns,column_count=len(header),
                            column_evidence='actual XLSB header and sheet checked against saved validation',
                            time_step_s=dt,oxygen_total_m3=total,flow_integral_m3=trapezoid(q,t),
                            integral_relative_error=(trapezoid(q,t)-total)/total))
        contexts[heat]=data
        for gas in ['co','co2']:
            y=clean[gas]; p=settings(cfg); s=scan(q,y,**p)
            prefix='delay_'+gas+'_'
            for key in ['center_index','delay_observed','max_corr','source_std','target_std','peak_separation','peak_prominence','peak_width_samples','boundary_hit','confidence','support_left_index','support_right_index']:
                dest={'delay_observed':'observed'}.get(key,key)
                ok=np.allclose(s[key],data[prefix+dest],rtol=1e-10,atol=1e-12,equal_nan=True)
                array_checks[gas+'_'+key]=bool(ok)
            array_checks[gas+'_reliable']=bool(np.array_equal(s['reliable'],data[prefix+'valid_mask']))
            m=fit_delay_models(q[s['center_index']],s['delay_observed'],s['confidence'],cfg['delay']['max_lag_s'])
            if compare_tree(m,summary['models'][gas]): raise AssertionError('configured model differs from baseline')
            if not all(array_checks.values()): raise AssertionError('experimental scan differs from production')
            groups=group_windows(s['center_index'],s['support_left_index'],s['support_right_index'],s['reliable'],p['stride'])
            checkcounts={k:dict(passed=int(v.sum()),rejected=int((~v).sum())) for k,v in s['checks'].items()}
            cumulative=np.ones(len(s['reliable']),bool); funnel={}
            for key,v in s['checks'].items(): cumulative &= v; funnel[key]=int(cumulative.sum())
            sr=dict(heat_id=heat,gas=gas,candidate_windows=len(s['reliable']),accepted_windows=int(s['reliable'].sum()),
                    check_counts=checkcounts,sequential_funnel=funnel,
                    consecutive_groups=len(set(groups[0][groups[0]>=0])),
                    source_overlap_groups=len(set(groups[1][groups[1]>=0])),
                    search_overlap_groups=len(set(groups[2][groups[2]>=0])),
                    zero_peak_count=int((s['best_delay_s']==0).sum()),
                    upper_boundary_peak_count=int((s['best_delay_s']==cfg['delay']['max_lag_s']).sum()),
                    fitted_models=summary['models'][gas])
            rows.append(sr)
            for i,c in enumerate(s['center_index']):
                checks={k:bool(v[i]) for k,v in s['checks'].items()}
                rec=dict(heat_id=heat,gas=gas,center_index=int(c),center_time_s=t[c],flow_m3_s=q[c],
                         oxygen_progress=(clean['oxygen_cumulative'][c]-clean['oxygen_cumulative'][0])/total,
                         stage=min(5,int(5*(clean['oxygen_cumulative'][c]-clean['oxygen_cumulative'][0])/total)+1),
                         best_delay_s=s['best_delay_s'][i],accepted=bool(s['reliable'][i]),
                         max_corr=s['max_corr'][i],peak_separation=s['peak_separation'][i],
                         peak_prominence=s['peak_prominence'][i],peak_width_samples=s['peak_width_samples'][i],
                         boundary_hit=bool(s['boundary_hit'][i]),
                         source_left_index=int(c-p['window']//2),source_right_index=int(c+p['window']//2),
                         search_left_index=int(s['support_left_index'][i]),search_right_index=int(s['support_right_index'][i]),
                         consecutive_group=int(groups[0,i]),source_overlap_group=int(groups[1,i]),search_overlap_group=int(groups[2,i]),
                         failure_reasons=';'.join(k for k,v in checks.items() if not v),
                         **{'pass_'+k:v for k,v in checks.items()})
                windows.append(rec)
                if rec['accepted']:
                    rec=dict(rec,event_identity='uncertain',identity_reason='correlation peak alone does not establish a unique input-response event')
                    accepted.append(rec)
            half=p['window']//2; L=p['max_lag']
            common=s['center_index'][(s['center_index']>=half+L)&(s['center_index']<len(t)-half-L)]
            for signed in [False,True]:
                ss=scan(q,y,**p,signed=signed,centers=common)
                for i,c in enumerate(common):
                    comparisons.append(dict(heat_id=heat,gas=gas,mode='signed' if signed else 'nonnegative',center_index=int(c),
                                            accepted=bool(ss['reliable'][i]),best_delay_s=ss['best_delay_s'][i],max_corr=ss['max_corr'][i]))
            gg=global_scan(q,y,dt,L); globals_.append(dict(heat_id=heat,gas=gas,**gg))
            for win in audit_cfg['synthetic']['windows']:
                pp=settings(cfg,win); half=win//2
                cc=np.arange(half+L,len(t)-half-L,pp['stride'])
                for signed in [False,True]:
                    ss=scan(q,y,**pp,signed=signed,centers=cc)
                    v=ss['delay_observed'][ss['reliable']]
                    sensitivity.append(dict(heat_id=heat,gas=gas,window=win,mode='signed' if signed else 'nonnegative',
                                            candidate_windows=len(cc),accepted_windows=len(v),
                                            median_delay_s=float(np.median(v)) if len(v) else None))
        comparisons_hashes={k:hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest() for k,v in data.items()}
        write_json(baselines/'results'/f'{heat}_array_hashes.json',comparisons_hashes)
        print('baseline verified',heat,'records',len(t),flush=True)
        # Numerical tolerance accommodates SciPy/BLAS platform variation; masks and centers are exact.
        write_json(baselines/'results'/f'{heat}_comparison.json',dict(summary_matches_saved=True,array_checks=array_checks,
            historical_npz_available=False,comparison_basis='configured replay versus cea8248-identical production code; historical JSON summary',
            historical_array_equality_claim=False))
    for row in accepted:
        others=[r for r in accepted if r['heat_id']==row['heat_id'] and r['gas']!=row['gas'] and
                r['source_left_index']<=row['source_right_index'] and r['source_right_index']>=row['source_left_index']]
        row['other_gas_overlapping_source_windows']=len(others)
        row['other_gas_same_center_count']=sum(r['center_index']==row['center_index'] for r in others)
    write_json(out/'input_audit.json',sources); write_json(out/'real_summary.json',rows)
    write_csv(out/'window_diagnostics.csv',windows); write_csv(out/'accepted_windows.csv',accepted)
    write_csv(out/'signed_comparison.csv',comparisons); write_csv(out/'window_sensitivity.csv',sensitivity)
    write_json(out/'global_correlations.json',globals_)
    return rows,accepted,contexts


def synthetic_cases(n,dt,seed,cfg):
    rng=np.random.default_rng(seed); time=np.arange(n)*dt
    source=np.convolve(rng.normal(size=n),np.ones(11)/11,mode='same')
    scale=np.std(source); base_noise=rng.normal(size=n)
    for delay in cfg['fixed_delays_s']:
        truth=np.full(n,delay)
        target=np.interp(time-delay,time,source)+base_noise*scale*.05
        yield f'fixed_{delay:+g}s',source,target,dict(onset=truth,peak=truth,centroid=truth),np.ones(n,bool),'fixed translation'
    for ratio in cfg['noise_ratios'][1:]:
        target=np.interp(time-15,time,source)+base_noise*scale*ratio; truth=np.full(n,15.)
        yield f'noise_{ratio:g}',source,target,dict(onset=truth,peak=truth,centroid=truth),np.ones(n,bool),'fixed translation'
    delay=5+20*time/time[-1]; mapping=time+delay
    target=np.interp(time,mapping,source)+base_noise*scale*.05
    yield 'drifting_5_to_25s',source,target,dict(onset=delay,peak=delay,centroid=delay),np.ones(n,bool),'delay defined on source-time axis; inverse mapping used'
    for theta in [2.,5.]:
        kt=np.arange(0,6*theta+dt/2,dt); kernel=kt*np.exp(-kt/theta); kernel/=kernel.sum()
        target=np.convolve(np.interp(time-5,time,source),kernel,mode='full')[:n]+base_noise*scale*.05
        # Causal kernel k(0)=0; onset is infimum of support; peak/centroid are distinct known estimands.
        truth={name:np.full(n,value) for name,value in [('onset',5.),('peak',5+kt[np.argmax(kernel)]),('centroid',5+np.sum(kt*kernel))]}
        yield f'broadened_theta_{theta:g}s',source,target,truth,np.ones(n,bool),'causal gamma-shaped impulse response; no single rigid-shift truth'
    event_times=np.array([100,118,230,248,350.]); delays=np.array([10,20,10,20,15.]); amps=rng.uniform(.8,1.2,len(event_times))
    x=sum(a*np.exp(-.5*((time-t)/2)**2) for a,t in zip(amps,event_times))
    y=sum(a*np.exp(-.5*((time-t-d)/4)**2) for a,t,d in zip(amps,event_times,delays))+base_noise*np.std(x)*.05
    nearest=np.argmin(abs(time[:,None]-event_times[None,:]),axis=1)
    truth=delays[nearest]; eligible=np.min(abs(time[:,None]-event_times[None,:]),axis=1)<=5
    yield 'adjacent_events',x,y,dict(onset=np.full(n,np.nan),peak=truth,centroid=truth),eligible,'Gaussian event center delays; Gaussian onset is undefined; only source centers within 5s of an event scored'
    unrelated=np.convolve(rng.normal(size=n),np.ones(11)/11,mode='same')
    yield 'unrelated',source,unrelated,dict(onset=np.full(n,np.nan),peak=np.full(n,np.nan),centroid=np.full(n,np.nan)),np.ones(n,bool),'null pair; no true delay'


def synthetic_audit(out,cfg,acfg):
    p=settings(cfg); scfg=acfg['synthetic']; dt=p['dt']; n=scfg['samples']; half=p['window']//2; L=p['max_lag']
    centers=np.arange(half+L,n-half-L,p['stride']); rows=[]; window_rows=[]; calibration=[]
    for seed in scfg['seeds']:
        for name,x,y,truth,eligible,note in synthetic_cases(n,dt,seed,scfg):
            for signed in [False,True]:
                s=scan(x,y,**p,signed=signed,centers=centers); mode='signed' if signed else 'nonnegative'
                for definition,actual in truth.items():
                    if name=='unrelated':
                        metrics=error_metrics(s['delay_observed'],np.zeros(len(centers)))
                        metrics['null_false_acceptance_rate']=float(s['reliable'].mean())
                        for k in ['mae_s','bias_s','median_absolute_error_s','p90_absolute_error_s','wrong_direction_rate_accepted_nonzero','within_1s_rate_eligible']: metrics[k]=None
                        metrics['wrong_direction_count']=None
                    else: metrics=error_metrics(s['delay_observed'],actual[centers],eligible[centers])
                    rows.append(dict(seed=seed,scenario=name,mode=mode,truth_definition=definition,note=note,**metrics))
                for i,c in enumerate(centers):
                    window_rows.append(dict(seed=seed,scenario=name,mode=mode,center_index=int(c),eligible=bool(eligible[c]),
                                            accepted=bool(s['reliable'][i]),estimate_s=float(s['delay_observed'][i]),
                                            onset_s=truth['onset'][c],peak_s=truth['peak'][c],centroid_s=truth['centroid'][c]))
                for i in np.flatnonzero(s['reliable']):
                    calibration.append(dict(seed=seed,scenario=name,mode=mode,center_index=int(centers[i]),
                                            peak_prominence=s['peak_prominence'][i],peak_width_samples=s['peak_width_samples'][i],
                                            estimate_s=s['delay_observed'][i],truth_peak_s=truth['peak'][centers[i]]))
        print('synthetic seed complete',seed,flush=True)
    write_csv(out/'synthetic_runs.csv',rows);write_csv(out/'synthetic_windows.csv',window_rows)
    write_csv(out/'synthetic_peak_calibration.csv',calibration)
    aggregate=[]
    for name in sorted({r['scenario'] for r in window_rows}):
        for mode in ['nonnegative','signed']:
            subset=[r for r in window_rows if r['scenario']==name and r['mode']==mode]
            estimates=np.array([r['estimate_s'] for r in subset]); elig=np.array([r['eligible'] for r in subset])
            for definition in ['onset','peak','centroid']:
                truth=np.array([r[definition+'_s'] for r in subset])
                if name=='unrelated':
                    result=dict(eligible_windows=len(subset),accepted_windows=int(np.isfinite(estimates).sum()),
                                null_false_acceptance_rate=float(np.isfinite(estimates).mean()),
                                acceptance_rate=float(np.isfinite(estimates).mean()),rejection_rate=float(1-np.isfinite(estimates).mean()))
                else: result=error_metrics(estimates,truth,elig)
                aggregate.append(dict(scenario=name,mode=mode,truth_definition=definition,**result))
    write_json(out/'synthetic_summary.json',aggregate)
    return aggregate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,default=ROOT/'Data')
    parser.add_argument('--config',type=Path,default=ROOT/'problem2/configs/delay_audit_config.json')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'problem2/final/delay_audit')
    parser.add_argument('--stage',choices=['all','real','synthetic','figures','report'],default='all')
    args=parser.parse_args(); out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    acfg=json.loads(args.config.read_text()); cfgpath=ROOT/acfg['pipeline_config'];cfg=json.loads(cfgpath.read_text())
    if cfg['reconstruction']!='PCHIP' or cfg['flow']['method']!='Savitzky-Golay derivative' or not cfg['flow']['clip_negative']: raise ValueError('unsupported baseline configuration')
    production=['q2_pipeline.py','q2_pipeline_v2.py','delay_v2.py','data_validation.py']
    for file in production:
        path=ROOT/'problem2/src'/file
        frozen=subprocess.check_output(['git','show',acfg['baseline_commit']+':problem2/src/'+file],cwd=ROOT)
        if path.read_bytes()!=frozen: raise RuntimeError('production differs from baseline: '+file)
    config_path=out/'effective_config.json'
    effective={'audit':acfg,'pipeline':cfg}
    if config_path.exists() and json.loads(config_path.read_text())!=effective and args.stage!='all':
        raise RuntimeError('configuration changed; rerun all stages to avoid mixed experiments')
    write_json(config_path,effective)
    signature=hashlib.sha256(json.dumps(effective,sort_keys=True).encode()+b''.join(p.read_bytes() for p in sorted((ROOT/'problem2/src').glob('*.py')))).hexdigest()
    status_path=out/'stage_status.json'
    state=json.loads(status_path.read_text()) if status_path.exists() else {}
    if state.get('signature')!=signature or args.stage=='all': state={'signature':signature,'completed':[]}
    def require(stages):
        if not set(stages)<=set(state['completed']):
            raise RuntimeError('missing current-code prerequisite stages: '+repr(stages)+'; rerun all stages')
    def completed(stage):
        if stage not in state['completed']: state['completed'].append(stage)
        write_json(status_path,state)
    if args.stage in ['all','real']:
        state['completed']=[s for s in state['completed'] if s not in ['real','figures','report']]
        write_json(status_path,state)
        real_audit(args.data_dir,out,cfg,acfg)
    if args.stage in ['all','synthetic']:
        state['completed']=[s for s in state['completed'] if s not in ['synthetic','figures','report']]
        write_json(status_path,state)
        synthetic_audit(out,cfg,acfg); completed('synthetic')
    if args.stage=='figures': require(['real','synthetic'])
    if args.stage in ['all','real','figures']:
        from .delay_audit_models import assess
        assess(out,cfg)
        if args.stage in ['all','real']: completed('real')
    if args.stage in ['all','figures']:
        require(['real','synthetic'])
        from .delay_audit_figures import render
        render(out,cfg);completed('figures')
    if args.stage in ['all','report']:
        require(['real','synthetic','figures'])
        from .delay_audit_report import update
        update(out);completed('report')
    provenance=dict(baseline_commit=subprocess.check_output(['git','rev-parse',acfg['baseline_commit']],cwd=ROOT,text=True).strip(),
                    execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
                    python=platform.python_version(),packages={p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','matplotlib','pyxlsb']},
                    source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'problem2/src').glob('*.py'))},
                    config_hashes={str(p):digest(p) for p in [args.config.resolve(),cfgpath]},
                    reference_hashes={str(p.relative_to(ROOT)):digest(p) for p in [ROOT/'problem2/audit/data_manifest.csv', ROOT/'problem2/final/results/q2_summary.json',ROOT/'problem2/final/results/process_validation.json',ROOT/'problem2/final/results/cross_validation.json']},
                    report_hash=digest(ROOT/'problem2/final/reports/独立数学审查与修复报告.md'),
                    data_source='XLSB process columns; first-question spectral arrays are not inputs',
                    units=cfg['units'],historical_npz_available=False,
                    reproduction='historical JSON summaries and input hashes checked; configured arrays checked against unchanged baseline code',
                    outputs=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='provenance.json'])
    write_json(out/'provenance.json',provenance)
    print('audit complete',out,flush=True)


if __name__=='__main__': main()
