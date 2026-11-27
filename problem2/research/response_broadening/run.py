"""Independent broadening study. Run with python -m problem2.research.response_broadening.run."""
from __future__ import annotations
import argparse,copy,csv,hashlib,importlib.metadata,json,platform,subprocess
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits

from .core import KernelDesign,changes,descriptors,detect_events,candidate_match,evaluate_matches,adjust_bh
from .benchmark import event_case
from ...src.delay_audit_core import scan
from ...src.delay_audit import serial,write_json,write_csv,synthetic_cases
from ...src.q2_pipeline import oxygen_flow

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
PRIOR=ROOT/'problem2/final/delay_audit'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def short_fit(fit):return {k:v for k,v in fit.items() if k not in ['prediction','baseline_prediction']}


def phase_surrogate(y,rng):
    freq=np.fft.rfft(y);phase=rng.uniform(0,2*np.pi,len(freq));phase[0]=0
    if len(y)%2==0:phase[-1]=0
    return np.fft.irfft(freq*np.exp(1j*phase),n=len(y))


def review_prior(out):
    refs=['synthetic_summary.json','real_summary.json','accepted_windows.csv','provenance.json']
    history=json.loads((PRIOR/'synthetic_summary.json').read_text())
    cases=sorted({r['scenario'] for r in history});pure=[c for c in cases if c.startswith('fixed_')]
    review=dict(commit='9b8513a',scenarios=cases,scenario_count=len(cases),
                pure_translation_reused_without_rerun=pure,
                nonideal_replay=[c for c in cases if c not in pure],
                limitation='Most prior traces have no sparse event ground truth. Event PR/F1 are not defined for them.',
                references={n:sha(PRIOR/n) for n in refs},
                previous_report_sha256=sha(ROOT/'problem2/final/reports/独立数学审查与修复报告.md'))
    write_json(out/'prior_review.json',review)
    return review


def reuse_prior(out,cfg):
    acfg=json.loads((PRIOR/'effective_config.json').read_text())['audit']['synthetic']
    acfg={**acfg,'fixed_delays_s':[]}  # No reexecution of ideal translation checks.
    rows=[]
    for seed in cfg['synthetic']['prior_reuse_seeds']:
        for name,q,g,truth,eligible,note in synthetic_cases(1000,cfg['dt_s'],seed,acfg):
            u,y=changes(q,cfg),changes(g,cfg);fit=KernelDesign(u,cfg).fit(y)
            true_kernel=None
            if name.startswith('broadened_theta_'):
                theta=float(name.split('_')[-1][:-1]);kt=np.arange(0,6*theta+cfg['dt_s']/2,cfg['dt_s'])
                k=kt*np.exp(-kt/theta);k/=k.sum();true_kernel=np.zeros(121);start=round(5/cfg['dt_s']);true_kernel[start:start+len(k)]=k
            elif name.startswith('noise_'):
                true_kernel=np.zeros(121);true_kernel[round(15/cfg['dt_s'])]=1
            td=descriptors(true_kernel,cfg['dt_s']) if true_kernel is not None else None
            sources=detect_events(u,cfg,KernelDesign(u,cfg).train);targets=detect_events(y,cfg,KernelDesign(u,cfg).train)
            er=candidate_match(sources,targets,cfg)
            row=dict(seed=seed,scenario=name,validation_skill=fit['validation_skill'],test_skill=fit['test_skill'],
                     kernel_usable=fit['usable_for_matching'],descriptors=fit['descriptors'],known_kernel=td,
                     event_truth_available=False,event_precision=None,event_recall=None,event_f1=None,
                     descriptive_event_matches=sum(r['status']=='matched' for r in er),
                     event_metric_reason='no compatible labelled change-event truth in the prior audit')
            for key in ['onset','peak','centroid','width']:
                row[key+'_absolute_error_s']=abs(fit['descriptors'][key]-td[key]) if td and fit['descriptors'][key] is not None else None
            rows.append(row)
    write_json(out/'prior_nonideal_kernel.json',rows)
    print('reused 7 nonideal scenarios, no ideal translation scans',flush=True)


def process_case(case,cfg):
    design=KernelDesign(case['source'],cfg);fit=design.fit(case['response'])
    sources=detect_events(case['source'],cfg,design.train);targets=detect_events(case['response'],cfg,design.train)
    v2=scan(case['flow'],case['gas'],dt=cfg['dt_s'])
    predictions={method:candidate_match(sources,targets,cfg,method,fit,v2) for method in ['v2','events','kernel']}
    return design,fit,sources,targets,predictions


def synthetic(out,cfg,development=False):
    seeds=cfg['synthetic']['development_seeds'] if development else cfg['synthetic']['evaluation_seeds']
    rows=[];eventrows=[];truthrows=[];kernelrows=[]
    for seed in seeds:
        for name in cfg['synthetic']['scenarios']:
            case=event_case(name,seed,cfg);design,fit,sources,targets,predictions=process_case(case,cfg)
            for method,records in predictions.items():
                metrics=evaluate_matches(records,sources,targets,case['truth'],fit['test_start_s'],cfg)
                rows.append(dict(seed=seed,scenario=name,method=method,**metrics))
                for r in records:
                    lo=max(0,int(r['source_start']/cfg['dt_s']));hi=min(len(case['flow']),int(r['source_end']/cfg['dt_s'])+1)
                    eventrows.append(dict(seed=seed,scenario=name,method=method,heat_id=f'synthetic_{seed}_{name}',
                                          source_flow=float(np.median(case['flow'][lo:hi])),independence_status='not_established',
                                          evaluation_partition='test' if r['source_peak']>=fit['test_start_s'] else 'training_or_guard',**r))
            td=descriptors(case['kernel_truth'],cfg['dt_s']) if case['kernel_truth'] is not None else None
            kernelrows.append(dict(seed=seed,scenario=name,validation_skill=fit['validation_skill'],test_skill=fit['test_skill'],
                                   usable=fit['usable_for_matching'],estimated=fit['descriptors'],truth=td,
                                   stationary_kernel=case['stationary_kernel'],
                                   errors={k:abs(fit['descriptors'][k]-td[k]) if td and fit['descriptors'][k] is not None else None for k in ['onset','peak','centroid','width']}))
            truthrows.extend(dict(seed=seed,scenario=name,**r) for r in case['truth'])
            if not development and seed==seeds[0]:
                examples=out/'examples';examples.mkdir(exist_ok=True)
                np.savez_compressed(examples/(name+'.npz'),time=case['time'],flow=case['flow'],gas=case['gas'],
                                    source=case['source'],response=case['response'],prediction=fit['prediction'],
                                    estimated_kernel=fit['kernel_weights'],true_kernel=case['kernel_truth'] if td else np.array([]))
        print('development' if development else 'evaluation','seed',seed,flush=True)
    sub=out/('development' if development else 'synthetic');sub.mkdir(exist_ok=True)
    write_csv(sub/'metrics.csv',rows);write_csv(sub/'event_matches.csv',eventrows);write_csv(sub/'truth_events.csv',truthrows)
    write_json(sub/'kernel_metrics.json',kernelrows)
    summary=[]
    for name in cfg['synthetic']['scenarios']:
        for method in ['v2','events','kernel']:
            subset=[r for r in rows if r['scenario']==name and r['method']==method]
            tp,fp,fn=(sum(r[k] for r in subset) for k in ['tp','fp','fn'])
            entry=dict(scenario=name,method=method,seeds=len(seeds),tp=tp,fp=fp,fn=fn,
                       precision=tp/(tp+fp) if tp+fp else None,recall=tp/(tp+fn) if tp+fn else None,
                       f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
                       false_match_rate=fp/(tp+fp) if tp+fp else None,
                       ambiguous_sources=sum(r['ambiguous_sources'] for r in subset),
                       matched_sources=sum(r['matched_sources'] for r in subset),
                       detected_sources=sum(r['detected_sources'] for r in subset),
                       truth_events=sum(r['truth_events'] for r in subset),
                       seed_f1_values=[r['f1'] for r in subset])
            for key in ['true_match_onset_mae_s','true_match_peak_mae_s','true_match_centroid_mae_s']:
                available=[r for r in subset if r[key] is not None]
                entry[key]=float(np.average([r[key] for r in available],weights=[r['tp'] for r in available])) if available else None
            entry['null_match_fraction']=entry['matched_sources']/entry['detected_sources'] if entry['detected_sources'] and not entry['truth_events'] else None
            summary.append(entry)
    write_json(sub/'summary.json',summary)


def real(out,cfg):
    provenance=json.loads((PRIOR/'provenance.json').read_text());hashes={r['path']:r['sha256'] for r in provenance['outputs']}
    prior_inputs=json.loads((PRIOR/'input_audit.json').read_text());rows=[];events=[];detections=[];sensitivity=[];nuisance=[];surrogates=[]
    for info in prior_inputs:
        heat=info['heat_id'];raw=Path(info['path'])
        if not raw.is_file() or sha(raw)!=info['sha256']:raise ValueError('raw file identity no longer matches prior audit: '+heat)
        npz=PRIOR/'baseline/data'/f'{heat}_q2_processed.npz'
        if not npz.is_file():raise FileNotFoundError('Reproduce prior audit NPZ first: '+str(npz))
        if sha(npz)!=hashes[str(npz.relative_to(PRIOR))]:raise ValueError('baseline NPZ hash mismatch')
        with np.load(npz,allow_pickle=False) as data:
            time=data['time'];flow=data['oxygen_flow'];oxygen=data['oxygen_cumulative_clean']
            for gas in ['co','co2']:
                series=data[gas+'_clean'];u,y=changes(flow,cfg),changes(series,cfg)
                design=KernelDesign(u,cfg);fit=design.fit(y)
                sr=detect_events(u,cfg,design.train);gr=detect_events(y,cfg,design.train)
                # Read the formal baseline windows; do not rerun or modify V2.
                v2=dict(center_time=data[f'delay_{gas}_window_time'],reliable=data[f'delay_{gas}_valid_mask'],delay_observed=data[f'delay_{gas}_observed'])
                for method in ['v2','events','kernel']:
                    matches=candidate_match(sr,gr,cfg,method,fit,v2)
                    for r in matches:
                        lo=max(0,int(r['source_start']/cfg['dt_s']));hi=min(len(flow),int(r['source_end']/cfg['dt_s'])+1)
                        events.append(dict(heat_id=heat,gas=gas,method=method,source_flow=float(np.median(flow[lo:hi])),
                                           independence_status='not_established',evidence_status='candidate_only',**r))
                detections.extend(dict(heat_id=heat,gas=gas,role=role,event_index=i,**r) for role,rr in [('source',sr),('gas',gr)] for i,r in enumerate(rr))
                rng=np.random.default_rng(cfg['real']['seed']+len(rows))
                null=[]
                for i in range(cfg['real']['surrogates']):
                    f=design.fit(phase_surrogate(y,rng),details=False);null.append(f['test_skill'])
                    surrogates.append(dict(heat_id=heat,gas=gas,replicate=i,test_skill=f['test_skill'],alpha=f['alpha']))
                p=(1+sum(v>=fit['test_skill'] for v in null))/(1+len(null)) if fit['test_skill'] is not None else 1.
                row=dict(heat_id=heat,gas=gas,source_candidates=len(sr),gas_candidates=len(gr),
                         phase_surrogate_p=p,**short_fit(fit));rows.append(row)
                for window in [21,31,41]:
                    alt=oxygen_flow(time,oxygen,window=window);af=KernelDesign(changes(alt,cfg),cfg).fit(y)
                    sensitivity.append(dict(heat_id=heat,gas=gas,sg_window=window,validation_skill=af['validation_skill'],test_skill=af['test_skill'],
                                            peak_delay_s=af['descriptors']['peak'],centroid_delay_s=af['descriptors']['centroid'],usable=af['usable_for_matching']))
                for order in cfg['real']['nuisance_order_sensitivity']:
                    alternate=copy.deepcopy(cfg);alternate['kernel']['nuisance_order']=order
                    af=KernelDesign(u,alternate).fit(y)
                    nuisance.append(dict(heat_id=heat,gas=gas,nuisance_order=order,test_r2=af['test_r2'],
                                         test_skill=af['test_skill'],validation_skill=af['validation_skill'],peak_delay_s=af['descriptors']['peak'],
                                         centroid_delay_s=af['descriptors']['centroid'],
                                         role='posthoc diagnostic; not used to select the reported main method'))
                example=out/'real_examples';example.mkdir(exist_ok=True)
                np.savez_compressed(example/f'{heat}_{gas}.npz',time=time,flow=flow,gas=series,source=u,response=y,
                                    prediction=fit['prediction'],baseline=fit['baseline_prediction'],kernel=fit['kernel_weights'])
                print('real validation',heat,gas,'skill',round(fit['test_skill'],3),'phase p',p,flush=True)
    qs=adjust_bh([r['phase_surrogate_p'] for r in rows])
    for r,q in zip(rows,qs):
        r['phase_surrogate_bh_q']=float(q)
        ss=[s for s in sensitivity if s['heat_id']==r['heat_id'] and s['gas']==r['gas']]
        r['sg_peak_range_s']=float(np.ptp([s['peak_delay_s'] for s in ss]))
        r['sg_centroid_range_s']=float(np.ptp([s['centroid_delay_s'] for s in ss]))
        r['screen_pass']=bool(r['usable_for_matching'] and r['test_skill']>0 and q<=cfg['real']['fdr_q'] and max(r['sg_peak_range_s'],r['sg_centroid_range_s'])<=5)
    for r in events:
        model=next(s for s in rows if s['heat_id']==r['heat_id'] and s['gas']==r['gas'])
        r['kernel_validation_screen_pass']=model['screen_pass']
        r['validated_physical_match']=False  # No independent response labels on real data.
    sub=out/'real';sub.mkdir(exist_ok=True)
    write_json(sub/'kernel_validation.json',rows);write_csv(sub/'event_matches.csv',events);write_csv(sub/'event_detections.csv',detections)
    write_csv(sub/'sg_sensitivity.csv',sensitivity);write_csv(sub/'phase_controls.csv',surrogates)
    write_csv(sub/'nuisance_sensitivity.csv',nuisance)


def identity(out):
    candidates=[]
    for folder in [ROOT/'Data',Path('/home/user/.codex/attachments')]:
        if folder.exists():
            candidates += [p for p in folder.rglob('*.xlsb*') if any(x in p.name for x in ['6B03205','附件1','Attachment 1.'])]
    files=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in candidates]
    prior=next(r for r in json.loads((PRIOR/'input_audit.json').read_text()) if r['heat_id']=='6B03205')
    write_json(out/'attachment1_identity.json',dict(status='pending_designated_new_file',existing_manifest_sha256=prior['sha256'],
        observed_files=files,existing_file_matches_prior=all(r['sha256']==prior['sha256'] for r in files),
        new_file_identity_verified=False,freeze_gate_pass=False,
        conclusion='The designated new Chinese-named attachment has not been located; equality to the old manifest does not resolve this requirement.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=['all','development','synthetic','real','finish'],default='all')
    args=parser.parse_args();cfg=json.loads((HERE/'config.json').read_text());out=HERE/'results';out.mkdir(exist_ok=True)
    frozen=['problem2/src/q2_pipeline.py','problem2/src/q2_pipeline_v2.py','problem2/src/delay_v2.py',
            'problem2/final/reports/独立数学审查与修复报告.md','problem2/final/delay_audit/synthetic_summary.json']
    for rel in frozen:
        original=subprocess.check_output(['git','show',cfg['baseline_commit']+':'+rel],cwd=ROOT)
        if (ROOT/rel).read_bytes()!=original:raise RuntimeError('protected baseline changed: '+rel)
    write_json(out/'effective_config.json',cfg);review_prior(out);identity(out)
    signature=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()+b''.join((HERE/p).read_bytes() for p in ['core.py','benchmark.py','run.py'])).hexdigest()
    status=out/'stage_status.json';state=json.loads(status.read_text()) if status.exists() else {}
    if state.get('signature')!=signature or args.stage=='all':state={'signature':signature,'completed':[]}
    def done(stage):
        if stage not in state['completed']:state['completed'].append(stage)
        write_json(status,state)
    with threadpool_limits(limits=1):
        if args.stage in ['all','development']:
            synthetic(out,cfg,development=True);done('development')
        if args.stage in ['all','synthetic']:
            reuse_prior(out,cfg);synthetic(out,cfg);done('synthetic')
        if args.stage in ['all','real']:
            real(out,cfg);done('real')
        if args.stage in ['all','finish']:
            if not {'development','synthetic','real'}<=set(state['completed']):raise RuntimeError('rerun numerical stages with current configuration/code before finish')
            from .release_audit import inspect
            inspect(out)
            from .report import finish
            finish(out,cfg);done('finish')
    write_json(out/'provenance.json',dict(baseline_commit=cfg['baseline_commit'],execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),python=platform.python_version(),
        packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','matplotlib','threadpoolctl']},
        config_sha256=sha(HERE/'config.json'),source_hashes={str(p.relative_to(ROOT)):sha(p) for p in HERE.glob('*.py')},
        prior_npz_manifest_sha256=sha(PRIOR/'provenance.json'),raw_input_audit_sha256=sha(PRIOR/'input_audit.json'),
        protected_hashes={rel:sha(ROOT/rel) for rel in frozen},
        outputs=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='provenance.json']))
    print('study stage complete:',args.stage,flush=True)


if __name__=='__main__':main()
