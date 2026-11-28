"""Run the locked paired study. Entry: OPENBLAS_NUM_THREADS=1 python -m ...run."""
import copy,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.signal import lfilter
from ...src.delay_audit import write_csv,write_json
from ..signal_fidelity.integral import generate,fit as original_fit
from ..signal_fidelity.operators import pulse,acquire,clean,fast_hampel
from ..signal_fidelity.run import load_heats
from .model import Design,bootstrap

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OUT=HERE/'results'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def configuration():
    cfg=json.loads((HERE.parent/'signal_fidelity/config.json').read_text())
    cfg.update(json.loads((HERE/'config.json').read_text()))
    return cfg


def signals(cfg,tau,T,kind,seed):
    if kind not in ['short','slow']:return generate(cfg,tau,T,kind,seed)
    p=cfg['integral_probe'];dt=p['fine_dt_s'];t=np.arange(round(p['duration_s']/dt)+1)*dt
    w=4 if kind=='short' else 30;q=np.full_like(t,5.)
    for center,amp in [(85,3),(160,-2),(245,3),(330,2)]:q+=pulse(t,center,w,amp)[0]
    o=25+cumulative_trapezoid(q,t,initial=0.);driver=.8*np.interp(t-tau,t,q,left=q[0])+.4
    a=np.exp(-dt/T);y=np.r_[driver[0],lfilter([1-a],[1.,-a],driver[:-1],zi=[a*driver[0]])[0]]
    stride=round(cfg['dt_s']/dt);return t[::stride],o[::stride],y[::stride]


def settings(cfg):
    return [(kind,d,T) for kind in ['matched','short','slow'] for d in [5,15,30] for T in [2,8]]+[(k,15,8) for k in ['constant','unrelated','two_stage']]


def historical(cfg):
    rows=[]
    for tau in [5,15,30]:
        for T in [2,8]:
            for seed in range(8101,8111):
                t,o,y=generate(cfg,tau,T,'matched',seed);y+=np.random.default_rng(seed).normal(0,.05,len(t))
                original,_,_=original_fit(t,o,y,cfg);r,*_=Design(t,o,cfg,'ordinary').fit(y)
                for k in ['estimate_s','T_s','gain','offset','validation_improvement','profile_width_s']:
                    assert np.isclose(original[k],r[k],rtol=1e-7,atol=1e-8),(k,original[k],r[k])
                assert original['accepted']==r['accepted']
                rows.append(dict(seed=seed,true_delay_s=tau,true_T_s=T,**r,correct=bool(r['accepted'] and abs(r['estimate_s']-tau)<=1)))
    assert sum(r['correct'] for r in rows)==41
    return rows


def controlled(cfg):
    rows=[];profiles=[];boots=[];sensitivity=[];examples=[]
    case=0
    for kind,tau,T in settings(cfg):
        for noise in [0,.05]:
            seeds=[cfg['evaluation_seeds'][0]] if noise==0 else cfg['evaluation_seeds']
            for seed in seeds:
                t,o,truth=signals(cfg,tau,T,kind,seed);y=truth+np.random.default_rng(seed).normal(0,noise,len(t))
                for acq in cfg['acquisition']:
                    measured=acquire(o,t,acq)
                    for gas in cfg['gas_inputs']:
                        yy=y if gas=='raw' else clean(t,y,fast_hampel(y))
                        for method in ['ordinary','weak']:
                            design=Design(t,measured,cfg,method);r,coef,loss,condition=design.fit(yy)
                            label=dict(kind=kind,true_delay_s=tau if kind not in ['constant','unrelated'] else None,
                                       true_T_s=T,noise_std=noise,seed=seed,acquisition=acq,gas_input=gas,method=method)
                            known=label['true_delay_s'] is not None
                            r.update(error_s=r['estimate_s']-tau if known else None,
                                     centroid_error_s=r['estimate_s']+r['T_s']-tau-T if known else None,
                                     correct=bool(known and r['accepted'] and abs(r['estimate_s']-tau)<=1))
                            # All matched cases; representative tau15,T8 for short/slow/mismatch.
                            ci_case=(kind=='matched' or (kind in ['short','slow','two_stage'] and tau==15 and T==8))
                            if acq=='ideal' and gas=='raw' and noise>0 and ci_case:
                                ci,samples=bootstrap(design,yy,r,cfg,cfg['bootstrap_seed']+seed+100*tau+T)
                                r.update(ci);r['ci_covers_true_delay']=bool(ci['ci_lower_s']<=tau<=ci['ci_upper_s'])
                                for j,s in enumerate(samples):boots.append(dict(**label,replicate=j,lag_s=s[0],T_s=s[1],gain=s[2],offset=s[3]))
                            rows.append(dict(**label,**r))
                            if acq=='ideal' and gas=='raw':
                                for lag,cc,ll,cn in zip(design.lags,coef,loss,condition):
                                    profiles.append(dict(**label,lag_s=lag,loss=ll,T_s=cc[0],gain=cc[1],offset=cc[2],scaled_design_condition=cn))
                            if acq=='ideal' and gas=='raw' and noise==.05 and kind=='matched' and method=='weak':
                                for span in cfg['sensitivity_supports_s']:
                                    sr,*_=Design(t,measured,cfg,'weak',span).fit(yy)
                                    sensitivity.append(dict(**label,support_s=span,**sr,correct=bool(sr['accepted'] and abs(sr['estimate_s']-tau)<=1)))
                case+=1
            print('controlled',kind,tau,T,noise,'complete',flush=True)
    return rows,profiles,boots,sensitivity


def summary(rows):
    df=pd.DataFrame(rows);rr=[]
    keys=['kind','noise_std','acquisition','gas_input','method']
    for key,g in df.groupby(keys,dropna=False):
        known=g.true_delay_s.notna();accepted=g.accepted;err=g.loc[known & accepted,'error_s'].abs()
        ci=g.dropna(subset=['ci_lower_s']) if 'ci_lower_s' in g else g.iloc[:0]
        rr.append(dict(zip(keys,key),conditions=len(g),accepted=int(accepted.sum()),rejected=int((~accepted).sum()),
                       correct=int(g.correct.sum()),correct_rate=float(g.correct.mean()) if known.any() else None,
                       accepted_MAE_s=float(err.mean()) if len(err) else None,
                       median_bias_s=float(g.loc[known,'error_s'].median()) if known.any() else None,
                       p90_absolute_error_s=float(g.loc[known,'error_s'].abs().quantile(.9)) if known.any() else None,
                       median_T_s=float(g.T_s.median()),median_validation_improvement=float(g.validation_improvement.median()),
                       ci_conditions=len(ci),ci_coverage=float(ci.ci_covers_true_delay.mean()) if len(ci) else None,
                       ci_median_width_s=float(ci.ci_width_s.median()) if len(ci) else None))
    return rr


def decision(cfg,rows):
    df=pd.DataFrame(rows);primary=df[(df.acquisition=='ideal')&(df.gas_input=='raw')]
    noisy=primary[(primary.kind=='matched')&(primary.noise_std==.05)]
    table=noisy.pivot(index=['seed','true_delay_s','true_T_s'],columns='method',values='correct').astype(float)
    paired=table.weak-table.ordinary
    cluster=paired.groupby('seed').mean().values;rng=np.random.default_rng(cfg['bootstrap_seed'])
    gain=np.mean(rng.choice(cluster,size=(10000,len(cluster))),axis=1);lower,upper=np.quantile(gain,[.025,.975])
    crit=cfg['success'];checks=[]
    def check(name,passed,**evidence):checks.append(dict(name=name,passed=bool(passed),**evidence))
    w0=primary[(primary.kind=='matched')&(primary.noise_std==0)&(primary.method=='weak')]
    check('noiseless_non_degradation',w0.correct.sum()==crit['noiseless_matched_correct'],weak_correct=int(w0.correct.sum()),conditions=len(w0))
    nw=noisy[noisy.method=='weak'];na=noisy[noisy.method=='ordinary']
    check('noisy_stable_paired_gain',nw.correct.mean()>=crit['noisy_matched_correct_rate_min'] and lower>crit['paired_seed_bootstrap_gain_lower_min'],
          ordinary_correct=int(na.correct.sum()),weak_correct=int(nw.correct.sum()),conditions_per_method=len(nw),
          gain=float(paired.mean()),paired_seed_bootstrap_lower=float(lower),paired_seed_bootstrap_upper=float(upper))
    for kind in ['short','slow']:
        g=primary[(primary.kind==kind)&(primary.noise_std==.05)]
        rates=g.groupby('method').correct.mean()
        check(kind+'_no_rate_degradation',rates.weak>=rates.ordinary-crit['time_scale_rate_degradation_max'],ordinary_rate=float(rates.ordinary),weak_rate=float(rates.weak))
    for kind in ['constant','unrelated']:
        g=primary[(primary.kind==kind)&(primary.noise_std==.05)];rates=g.groupby('method').accepted.mean()
        passes=(rates.weak<=rates.ordinary+crit['negative_acceptance_increase_max'])
        if kind=='constant':passes=passes and g[g.method=='weak'].accepted.sum()<=crit['constant_accepted_max']
        check(kind+'_false_acceptance',passes,ordinary_rate=float(rates.ordinary),weak_rate=float(rates.weak),conditions_per_method=10)
    check('conditional_interval_calibration',nw.ci_covers_true_delay.mean()>=crit['conditional_ci_coverage_min'] and nw.ci_width_s.median()<=crit['conditional_ci_median_width_max_s'],
          coverage=float(nw.ci_covers_true_delay.mean()),median_width_s=float(nw.ci_width_s.median()),conditions=len(nw),
          note='Coverage over all 60 primary noisy matched trials, not only accepted trials; seeds shared across six settings.')
    return dict(controlled_gate_passed=all(x['passed'] for x in checks),checks=checks,
                real_data_status='pending' if all(x['passed'] for x in checks) else 'not_run_controlled_gate_failed',
                adoption='research_only_pending_real_validation' if all(x['passed'] for x in checks) else 'do_not_adopt',
                automatic_formal_replacement=False)


def main():
    started=time.time();OUT.mkdir(exist_ok=True);cfg=configuration()
    tracked=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
    protected={p:sha(ROOT/p) for p in tracked if p and (ROOT/p).is_file() and not p.startswith('problem2/research/weak_integral/')}
    code={str(p.relative_to(ROOT)):sha(p) for p in HERE.glob('*.py')};code[str((HERE/'config.json').relative_to(ROOT))]=sha(HERE/'config.json')
    lock=OUT/'protocol_lock.json'
    if lock.exists():
        locked=json.loads(lock.read_text());assert locked['effective_config']==cfg, 'configuration changed after evaluation'
        model_path=str((HERE/'model.py').relative_to(ROOT))
        assert locked['code_hashes'][model_path]==code[model_path], 'estimator changed after evaluation'
    else:
        write_json(lock,dict(code_hashes=code,effective_config=cfg,base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                   rule='Written before historical replay and before new evaluation; no post-evaluation tuning. Gate is evaluated once.'))
    write_json(OUT/'execution_manifest.json',dict(code_hashes=code,original_protocol_sha256=sha(lock),
               note='Original pre-evaluation lock retained. Presentation table adapter and interval rendering repaired after first numeric run; model and configuration unchanged.'))
    historical_rows=historical(cfg);write_csv(OUT/'historical_replay.csv',historical_rows)
    rows,profiles,boots,sensitivity=controlled(cfg)
    for name,data in [('trials',rows),('profiles',profiles),('bootstrap_draws',boots),('support_sensitivity',sensitivity),('summary',summary(rows))]:write_csv(OUT/(name+'.csv'),data)
    result=decision(cfg,rows)
    if result['controlled_gate_passed']:
        from .real import execute
        result.update(execute(cfg,OUT))
    write_json(OUT/'decision.json',result)
    from .present import deliver
    deliver(cfg,OUT,result)
    changed=[p for p,h in protected.items() if not (ROOT/p).exists() or sha(ROOT/p)!=h]
    assert not changed,changed
    write_json(OUT/'protection_check.json',dict(preexisting_tracked_files=len(protected),all_unchanged=True,changed=changed))
    write_json(OUT/'provenance.json',dict(base_commit=cfg['research_baseline'],runtime_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
               code_hashes=code,output_hashes={str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='provenance.json'},
               wall_seconds=time.time()-started))
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
