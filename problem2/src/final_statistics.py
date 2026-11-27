"""Final descriptive oxygen relationships; no changes to V2 or NPZ fields."""
from __future__ import annotations
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import pearsonr,spearmanr

RELATIONS=('co_instantaneous','co2_instantaneous','co_cumulative','co2_cumulative')


def correlations(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);keep=np.isfinite(x)&np.isfinite(y);x=x[keep];y=y[keep]
    if len(x)<3 or np.ptp(x)<=1e-12 or np.ptp(y)<=1e-12:
        return dict(n=len(x),pearson_r=None,spearman_r=None,status='constant or too short')
    return dict(n=len(x),pearson_r=float(pearsonr(x,y).statistic),spearman_r=float(spearmanr(x,y).statistic),status='descriptive')


def interval_increments(time,x,y,seconds):
    """Nonoverlapping equal-time increments; incomplete final block excluded."""
    grid=np.arange(time[0],time[-1]+1e-9,seconds)
    return np.diff(np.interp(grid,time,x)),np.diff(np.interp(grid,time,y))


def stage_statistics(d,edges):
    time=d['time'];oxygen=d['oxygen_cumulative_clean'];delta=oxygen[-1]-oxygen[0]
    if delta<=0 or np.any(np.diff(oxygen)<-1e-10):raise ValueError('oxygen progress requires nondecreasing cumulative input')
    p=(oxygen-oxygen[0])/delta
    # np.interp with nondecreasing repeats uses the last observation at a tie.
    boundaries=np.interp(edges,p,time);boundaries[0]=time[0];boundaries[-1]=time[-1]
    rows=[]
    for k,(lo,hi) in enumerate(zip(boundaries[:-1],boundaries[1:])):
        mask=(p>=edges[k])&((p<=edges[k+1]) if k==len(edges)-2 else (p<edges[k+1]))
        for gas in ['co','co2']:
            emission=d['total_'+gas+'_emission'];de=float(np.interp(hi,time,emission)-np.interp(lo,time,emission))
            do=float(delta*(edges[k+1]-edges[k]));duration=float(hi-lo)
            rows.append(dict(heat_id=d['heat_id'],gas=gas,stage=k+1,progress_start=edges[k],progress_end=edges[k+1],
                time_start_s=float(lo),time_end_s=float(hi),duration_s=duration,sample_count=int(mask.sum()),
                mean_rate_m3_s=de/duration if duration>0 else None,median_rate_m3_s=float(np.median(d[gas+'_clean'][mask])) if mask.any() else None,
                cumulative_emission_start_m3=float(np.interp(lo,time,emission)),cumulative_emission_end_m3=float(np.interp(hi,time,emission)),
                stage_emission_m3=de,stage_emission_share=de/float(emission[-1]) if emission[-1]>0 else None,
                stage_oxygen_m3=do,emission_per_oxygen_m3_m3=de/do))
    return rows


def target(d,relation):
    gas=relation.split('_')[0]
    return d[gas+'_clean'] if relation.endswith('instantaneous') else d['total_'+gas+'_emission']


def progress(d):
    o=d['oxygen_cumulative_clean'];return (o-o[0])/(o[-1]-o[0])


def fit_trend(x,y,degree,weights=None):
    design=np.vander(x,degree+1,increasing=True)
    if weights is not None:
        root=np.sqrt(weights);design=design*root[:,None];y=np.asarray(y)*root
    return np.linalg.lstsq(design,y,rcond=None)[0]


def scores(y,pred):
    error=np.asarray(y)-pred;den=np.sum((y-np.mean(y))**2)
    return dict(mae=float(np.mean(abs(error))),rmse=float(np.sqrt(np.mean(error**2))),r2=float(1-np.sum(error**2)/den) if den>0 else None)


def analyze(heats,cfg):
    levels=[];diffs=[];stages=[];trends=[];loo=[];pooled=[]
    for d in heats:
        stages.extend(stage_statistics(d,np.array(cfg['oxygen_progress_edges'])))
        o=d['oxygen_cumulative_clean'];t=d['time'];p=progress(d)
        for relation in RELATIONS:
            y=target(d,relation);row=dict(heat_id=d['heat_id'],relation=relation,**correlations(o,y))
            if relation.endswith('cumulative'):
                constant_rate=(t-t[0])/(t[-1]-t[0])*y[-1]
                row.update(constant_gas_rate_pearson=correlations(o,constant_rate)['pearson_r'],constant_gas_rate_spearman=correlations(o,constant_rate)['spearman_r'])
            levels.append(row)
            for seconds in cfg['increment_intervals_s']:
                dx,dy=interval_increments(t,o,y,seconds)
                diffs.append(dict(heat_id=d['heat_id'],relation=relation,interval_s=seconds,**correlations(dx,dy)))
            for degree in cfg['trend_degrees']:
                coef=fit_trend(p,y,degree);pred=np.polynomial.polynomial.polyval(p,coef)
                trends.append(dict(heat_id=d['heat_id'],relation=relation,degree=degree,c0=coef[0],c1=coef[1],c2=coef[2] if degree==2 else 0.,**scores(y,pred)))
    for relation in RELATIONS:
        for degree in cfg['trend_degrees']:
            for test in heats:
                train=[d for d in heats if d['heat_id']!=test['heat_id']]
                x=np.concatenate([progress(d) for d in train]);y=np.concatenate([target(d,relation) for d in train]);w=np.concatenate([np.full(len(d['time']),1/len(d['time'])) for d in train])
                coef=fit_trend(x,y,degree,w);prediction=np.polynomial.polynomial.polyval(progress(test),coef)
                loo.append(dict(test_heat=test['heat_id'],training_heats=len(train),relation=relation,degree=degree,c0=coef[0],c1=coef[1],c2=coef[2] if degree==2 else 0.,**scores(target(test,relation),prediction)))
            x=np.concatenate([progress(d) for d in heats]);y=np.concatenate([target(d,relation) for d in heats]);w=np.concatenate([np.full(len(d['time']),1/len(d['time'])) for d in heats])
            coef=fit_trend(x,y,degree,w)
            pooled.append(dict(relation=relation,degree=degree,c0=coef[0],c1=coef[1],c2=coef[2] if degree==2 else 0.,role='retrospective descriptive fit; not a physical law'))
    rng=np.random.default_rng(cfg['random_seed']);aggregate=[]
    for gas in ['co','co2']:
        for stage in range(1,6):
            rr=[r for r in stages if r['gas']==gas and r['stage']==stage];row=dict(gas=gas,stage=stage,heats=len(rr))
            # Resample entire heat rows, never autocorrelated time points.
            ids=rng.integers(0,len(rr),(cfg['bootstrap_heat_resamples'],len(rr)))
            for key in ['mean_rate_m3_s','emission_per_oxygen_m3_m3','stage_emission_share']:
                values=np.array([r[key] for r in rr]);boot=np.median(values[ids],axis=1)
                row.update({key+'_median':float(np.median(values)),key+'_q25':float(np.quantile(values,.25)),key+'_q75':float(np.quantile(values,.75)),key+'_ci_low':float(np.quantile(boot,.025)),key+'_ci_high':float(np.quantile(boot,.975))})
            aggregate.append(row)
    return dict(level_correlations=levels,increment_correlations=diffs,stage_statistics=stages,stage_aggregate=aggregate,trend_parameters=trends,trend_loho=loo,pooled_trends=pooled)
