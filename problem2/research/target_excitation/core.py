"""Paired calculations share all correlations; only excitation placement differs."""
import warnings
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import peak_prominences,peak_widths
from ...src.q2_pipeline import _standardized_change


def compare(source,target,cfg,centers=None):
    x,y=np.asarray(source,float),np.asarray(target,float)
    if x.ndim!=1 or x.shape!=y.shape or not np.isfinite([x,y]).all():raise ValueError('finite equal 1D inputs required')
    w=cfg['window'];L=cfg['max_lag_samples'];half=w//2;dt=cfg['dt_s'];stride=cfg['stride']
    if w%2!=1 or len(x)<=w+L:raise ValueError('window/search unsupported')
    c=np.arange(half,len(x)-half-L,stride) if centers is None else np.asarray(centers,int)
    if len(c)==0 or np.any(c<half) or np.any(c>=len(x)-half-L):raise ValueError('invalid shared centers')
    sx,sy=_standardized_change(x),_standardized_change(y);sg,tg=sx.std(),sy.std();sw,tw=sliding_window_view(sx,w),sliding_window_view(sy,w)
    left=c-half;a=sw[left];astd=a.std(axis=1);bstd0=tw[left].std(axis=1)
    source_ok=(sg>0)&(tg>0)&(astd>=cfg['min_strength_ratio']*sg)
    ac=(a-a.mean(axis=1)[:,None])/(astd[:,None]+1e-12)
    raw=np.full((len(c),L+1),np.nan);shifted_std=np.zeros_like(raw)
    for lag in range(L+1):
        b=tw[left+lag];std=b.std(axis=1);shifted_std[:,lag]=std;keep=source_ok&(std>0)
        if keep.any():raw[keep,lag]=np.mean(ac[keep]*(b[keep]-b[keep].mean(axis=1)[:,None])/(std[keep,None]+1e-12),axis=1)
    A_gate=source_ok&(bstd0>=cfg['min_strength_ratio']*tg)
    B_gate=source_ok[:,None]&(shifted_std>=cfg['min_strength_ratio']*tg)&(shifted_std>0)
    results={}
    for method,cc in [('A',np.where(A_gate[:,None],raw,np.nan)),('B',np.where(B_gate,raw,np.nan))]:
        finite=np.isfinite(cc).any(axis=1);winner=np.argmax(np.where(np.isfinite(cc),cc,-np.inf),axis=1)
        score=cc[np.arange(len(c)),winner];best_std=shifted_std[np.arange(len(c)),winner]
        separation=np.full(len(c),np.nan);full_sep=separation.copy();prominence=separation.copy();width=separation.copy();masked_neighbor=np.zeros(len(c),bool)
        for i in np.flatnonzero(finite):
            far=abs(np.arange(L+1)-winner[i])>=cfg['peak_distance_samples'];eligible=far&np.isfinite(cc[i]);allfinite=far&np.isfinite(raw[i])
            separation[i]=score[i]-(np.max(cc[i,eligible]) if eligible.any() else score[i])
            full_sep[i]=score[i]-(np.max(raw[i,allfinite]) if allfinite.any() else score[i])
            if 0<winner[i]<L:
                # Never fill masked gaps with zero or -inf for a shape diagnostic.
                diagnostic=cc[i] if method=='A' else raw[i]
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    prominence[i]=peak_prominences(diagnostic,[winner[i]])[0][0]
                    width[i]=peak_widths(diagnostic,[winner[i]],rel_height=.5)[0][0]
                masked_neighbor[i]=not np.isfinite(cc[i,winner[i]-1:winner[i]+2]).all()
            else:prominence[i]=0.;width[i]=np.inf
        boundary=(winner==0)|(winner==L)
        accepted=finite&(score>=cfg['min_corr'])&(separation>=cfg['min_separation'])&~boundary
        target_strength=bstd0 if method=='A' else best_std
        conf=np.clip((score-cfg['min_corr'])/(1-cfg['min_corr']),0,1)*np.clip(separation/.25,0,1)*np.minimum(np.minimum(astd/(sg+1e-12),target_strength/(tg+1e-12)),1)
        results[method]=dict(center_index=c,center_time=c*dt,winner_lag_samples=winner,best_delay_s=winner*dt,
            observed=np.where(accepted,winner*dt,np.nan),accepted=accepted,score=score,separation=separation,full_competitor_separation=full_sep,
            B_full_competitor_accept=accepted&(full_sep>=cfg['min_separation']),source_std=astd,unshifted_target_std=bstd0,best_shifted_target_std=best_std,
            source_strength_ratio=astd/(sg+1e-12),unshifted_target_ratio=bstd0/(tg+1e-12),best_shifted_target_ratio=best_std/(tg+1e-12),
            source_gate=source_ok,unshifted_target_gate=A_gate,has_eligible_lag=finite,boundary=boundary,
            confidence=np.where(accepted,conf,0.),prominence=prominence,width=width,masked_neighbor=masked_neighbor,
            correlations=cc,raw_correlations=raw,target_global_std=tg,source_global_std=sg,shifted_std=shifted_std)
    return results


def packet(n,index,rng,dt=.5):
    """Nonnegative compact transient; random internal shape, eight-second support."""
    radius=int(4/dt);z=np.zeros(n);xx=np.linspace(-1,1,2*radius+1)
    shape=np.cos(xx*np.pi/2)**2*np.maximum(.1,1+.35*rng.normal(size=len(xx)))
    if index-radius<0 or index+radius>=n:raise ValueError('packet outside record')
    z[index-radius:index+radius+1]=shape/shape.max();return z


def shift(x,delay,dt):
    steps=int(round(delay/dt));out=np.zeros_like(x)
    if steps==0:return x.copy()
    if steps>0:out[steps:]=x[:-steps]
    else:out[:steps]=x[-steps:]
    return out


def trial_row(pair,method,i,truth,cfg,**labels):
    r=pair[method];accepted=bool(r['accepted'][i]);estimate=float(r['observed'][i]) if accepted else None
    correct=bool(accepted and truth is not None and abs(estimate-truth)<=cfg['error_tolerance_s'])
    return dict(**labels,method=method,center_index=int(r['center_index'][i]),accepted=accepted,truth_s=truth,estimate_s=estimate,
        correct=correct,false_delay=bool(accepted and truth is not None and not correct),absolute_error_s=abs(estimate-truth) if accepted and truth is not None else None,
        source_ratio=float(r['source_strength_ratio'][i]),unshifted_target_ratio=float(r['unshifted_target_ratio'][i]),
        selected_target_ratio=float(r['best_shifted_target_ratio'][i]),max_corr=float(r['score'][i]),separation=float(r['separation'][i]),
        full_competitor_separation=float(r['full_competitor_separation'][i]),full_competitor_accept=bool(r['B_full_competitor_accept'][i]),
        source_gate=bool(r['source_gate'][i]),unshifted_target_gate=bool(r['unshifted_target_gate'][i]),has_eligible_lag=bool(r['has_eligible_lag'][i]),
        boundary=bool(r['boundary'][i]),confidence=float(r['confidence'][i]),masked_neighbor=bool(r['masked_neighbor'][i]))
