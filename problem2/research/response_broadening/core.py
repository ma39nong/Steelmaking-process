"""Independent response-kernel and conservative event-matching research."""
from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.linalg import cho_factor, cho_solve
from scipy.signal import savgol_filter, find_peaks, peak_widths
from scipy.optimize import linear_sum_assignment


def changes(values, cfg):
    """Offline SG derivative. No replacement of the formal oxygen-flow estimate."""
    x=np.asarray(values,float)
    if x.ndim!=1 or not np.isfinite(x).all(): raise ValueError('finite one-dimensional series required')
    return savgol_filter(x,cfg['change_smoothing_samples'],2,deriv=1,delta=cfg['dt_s'],mode='interp')


def descriptors(weights, dt):
    """Time descriptors of |response|; centroid is not a chemical yield."""
    w=np.abs(np.asarray(weights,float));total=w.sum()
    if not np.isfinite(w).all() or total<=1e-12:
        return dict(onset=None,peak=None,centroid=None,width=None,polarity=0,minority_mass=None,tail_mass=None)
    c=np.cumsum(w)/total;t=np.arange(len(w))*dt
    onset=float(t[np.searchsorted(c,.05)]);end=float(t[np.searchsorted(c,.95)])
    signed=np.asarray(weights);positive=np.maximum(signed,0).sum();negative=np.maximum(-signed,0).sum()
    return dict(onset=onset,peak=float(t[np.argmax(w)]),centroid=float(t@w/total),width=end-onset,
                polarity=1 if positive>=negative else -1,minority_mass=float(min(positive,negative)/total),
                tail_mass=float(w[-max(1,int(2/dt)):].sum()/total))


class KernelDesign:
    """Signed FIR with second-difference regularization and training-only scaling.

    Train/validation/test are blocked and separated by the full kernel horizon
    plus the SG half-window. Lambda is selected on validation only. Test outputs
    never enter scaling, nuisance fitting, coefficient fitting or lambda choice.
    """
    def __init__(self, u, cfg):
        self.cfg=cfg;self.u=np.asarray(u,float);n=len(u);self.L=int(round(cfg['max_response_s']/cfg['dt_s']))
        gap=self.L+cfg['change_smoothing_samples']//2
        train_end=int(n*cfg['kernel']['train_fraction']);val_end=int(n*cfg['kernel']['validation_end_fraction'])
        self.rows=np.arange(self.L,n)
        self.train=self.rows[self.rows<train_end]
        self.val=self.rows[(self.rows>=train_end+gap)&(self.rows<val_end)]
        self.test=self.rows[self.rows>=val_end+gap]
        if min(map(len,[self.train,self.val,self.test]))<20:raise ValueError('record too short for blocked response validation')
        scale=float(np.std(np.asarray(u)[self.train]))
        self.scale=max(scale,1e-12)
        windows=sliding_window_view(self.u,self.L+1)[:,::-1]/self.scale
        self.X=np.zeros((n,self.L+1));self.X[self.rows]=windows
        t=np.arange(n)/n;self.B=np.column_stack([np.ones(n),t])[:,:cfg['kernel'].get('nuisance_order',1)+1]
        self.nuisance=np.linalg.lstsq(self.B[self.train],self.X[self.train],rcond=None)[0]
        xr=self.X[self.train]-self.B[self.train]@self.nuisance
        self.xr=xr;self.gram=xr.T@xr/len(self.train)
        d2=np.diff(np.eye(self.L+1),n=2,axis=0)
        self.solvers=[cho_factor(self.gram+a*(d2.T@d2)+cfg['kernel']['ridge_floor']*np.eye(self.L+1),lower=True)
                      for a in cfg['kernel']['alphas']]
        ev=np.linalg.eigvalsh(self.gram)
        self.rank=int(np.sum(ev>max(ev.max(),1e-12)*1e-8))
        self.condition=float(ev.max()/max(ev.min(),1e-16))

    def fit(self,y, details=True):
        y=np.asarray(y,float)
        if y.shape!=self.u.shape or not np.isfinite(y).all():raise ValueError('response shape/values invalid')
        baseline_beta=np.linalg.lstsq(self.B[self.train],y[self.train],rcond=None)[0]
        residual=y[self.train]-self.B[self.train]@baseline_beta
        rhs=self.xr.T@residual/len(self.train)
        solutions=[];val_mse=[];betas=[]
        for factor in self.solvers:
            h=cho_solve(factor,rhs);beta=np.linalg.lstsq(self.B[self.train],y[self.train]-self.X[self.train]@h,rcond=None)[0]
            solutions.append(h);betas.append(beta)
            val_mse.append(float(np.mean((y[self.val]-self.X[self.val]@h-self.B[self.val]@beta)**2)))
        best=int(np.argmin(val_mse));h=solutions[best];beta=betas[best]
        prediction=self.X@h+self.B@beta
        def skill(indices):
            base=float(np.mean((y[indices]-self.B[indices]@baseline_beta)**2))
            return float(1-np.mean((y[indices]-prediction[indices])**2)/base) if base>1e-15 else None
        result=dict(alpha=self.cfg['kernel']['alphas'][best],validation_skill=skill(self.val),test_skill=skill(self.test),
                    test_r2=float(1-np.mean((y[self.test]-prediction[self.test])**2)/max(np.var(y[self.test]),1e-15)),
                    test_mae=float(np.mean(abs(y[self.test]-prediction[self.test]))),
                    test_response_std=float(np.std(y[self.test])),kernel_weights=h/self.scale,train_rows=len(self.train),validation_rows=len(self.val),test_rows=len(self.test),
                    test_start_s=float(self.test[0]*self.cfg['dt_s']),effective_rank=self.rank,unregularized_condition=self.condition)
        if not details:return result
        dt=self.cfg['dt_s'];desc=descriptors(h,dt)
        near=[i for i,mse in enumerate(val_mse) if mse<=max(val_mse[best],1e-15)*self.cfg['kernel']['stability_mse_ratio']]
        ds=[descriptors(solutions[i],dt) for i in near]
        peak_range=max(d['peak'] for d in ds)-min(d['peak'] for d in ds) if all(d['peak'] is not None for d in ds) else None
        centroid_range=max(d['centroid'] for d in ds)-min(d['centroid'] for d in ds) if all(d['centroid'] is not None for d in ds) else None
        stable=len(near)>=2 and peak_range is not None and max(peak_range,centroid_range)<=self.cfg['kernel']['stability_tolerance_s'] and len({d['polarity'] for d in ds})==1
        structural=(desc['minority_mass'] is not None and desc['minority_mass']<=self.cfg['kernel']['maximum_minority_mass'] and desc['tail_mass']<=self.cfg['kernel']['maximum_tail_mass'])
        result.update(descriptors=desc,near_optimal_alphas=[self.cfg['kernel']['alphas'][i] for i in near],
                      peak_range_s=peak_range,centroid_range_s=centroid_range,stable=stable,
                      structural_gate=bool(structural),
                      usable_for_matching=bool(stable and structural and result['validation_skill'] is not None and result['validation_skill']>=self.cfg['kernel']['minimum_validation_skill']),
                      candidate_descriptors=[descriptors(h,dt) for h in solutions],validation_mse=val_mse,
                      prediction=prediction,baseline_prediction=self.B@baseline_beta)
        return result


def detect_events(signal, cfg, scale_indices=None):
    """Detect signed extrema; intervals split at midpoint to prevent double use.

    Detection uses no gas series or response labels. Threshold scale is fitted
    on a supplied training interval, then fixed for the complete record.
    """
    x=np.asarray(signal,float);dt=cfg['dt_s'];ec=cfg['events'];n=len(x)
    ref=x if scale_indices is None else x[scale_indices]
    noise=1.4826*np.median(abs(np.diff(ref)-np.median(np.diff(ref))))/np.sqrt(2) if len(ref)>1 else 0
    threshold=max(ec['noise_multiplier']*noise,ec['relative_prominence']*np.max(abs(ref)),1e-10)
    records=[]
    for polarity in [1,-1]:
        oriented=polarity*x
        peaks,properties=find_peaks(oriented,height=threshold,prominence=threshold,distance=max(1,int(ec['minimum_distance_s']/dt)),
                                    width=max(1,ec['minimum_width_s']/dt))
        for j,p in enumerate(peaks):
            # Peak prominence bases may include other extrema. Boundaries are
            # additionally clipped to each signed lobe and neighboring events.
            left=int(p);right=int(p)
            while left>0 and oriented[left-1]>0:left-=1
            while right<n-1 and oriented[right+1]>0:right+=1
            records.append(dict(peak_index=int(p),left_index=left,right_index=right,polarity=polarity,
                                prominence=float(properties['prominences'][j]),threshold=threshold))
    records.sort(key=lambda e:e['peak_index'])
    for j,e in enumerate(records):
        if j:e['left_index']=max(e['left_index'],(e['peak_index']+records[j-1]['peak_index'])//2+1)
        if j+1<len(records):e['right_index']=min(e['right_index'],(e['peak_index']+records[j+1]['peak_index'])//2)
        lo,hi=e['left_index'],e['right_index'];w=np.maximum(e['polarity']*x[lo:hi+1],0)
        d=descriptors(w,dt)
        e.update(start=float(lo*dt+d['onset']),peak=float(e['peak_index']*dt),
                 centroid=float(lo*dt+d['centroid']),width=d['width'],end=float(hi*dt),
                 confidence=float(e['prominence']/(e['prominence']+threshold)))
    return records


def candidate_match(sources, targets, cfg, method='events', kernel=None, v2=None, v2_center_tolerance_s=2.5):
    """Only isolated edges are accepted. Any competing edge preserves ambiguity.

    event_id is a candidate identifier, not a claim of statistical independence.
    A target with two possible sources is never silently assigned to one.
    """
    edges=[];expected={};L=cfg['max_response_s'];tol=cfg['events']['kernel_peak_tolerance_s']
    for i,s in enumerate(sources):
        if method=='kernel':
            if not kernel or not kernel['usable_for_matching']:continue
            expected[i]=kernel['descriptors']['peak']
        elif method=='v2':
            if v2 is None:continue
            keep=v2['reliable']&(abs(v2['center_time']-s['peak'])<=v2_center_tolerance_s)
            values=v2['delay_observed'][keep]
            if not len(values):continue
            if np.ptp(values)>tol:continue
            expected[i]=float(np.median(values))
        for j,g in enumerate(targets):
            polarity=s['polarity']*(kernel['descriptors']['polarity'] if method=='kernel' else 1)
            delay=g['peak']-s['peak']
            if g['polarity']!=polarity or not (0<=delay<=L):continue
            if method in ['kernel','v2'] and abs(delay-expected[i])>tol:continue
            edges.append((i,j))
    source_degree=np.bincount([i for i,j in edges],minlength=len(sources));target_degree=np.bincount([j for i,j in edges],minlength=len(targets))
    records=[]
    for i,s in enumerate(sources):
        candidates=[j for ii,j in edges if ii==i]
        unique=len(candidates)==1 and target_degree[candidates[0]]==1
        state='matched' if unique else ('ambiguous' if candidates else 'unmatched')
        j=candidates[0] if unique else None;g=targets[j] if j is not None else None
        records.append(dict(event_id=f'source_{i:04d}',source_index=i,target_index=j,status=state,
                            source_start=s['start'],source_peak=s['peak'],source_centroid=s['centroid'],source_end=s['end'],
                            source_polarity=s['polarity'],gas_start=g['start'] if g else None,gas_peak=g['peak'] if g else None,
                            delay_onset=g['start']-s['start'] if g else None,delay_peak=g['peak']-s['peak'] if g else None,
                            delay_centroid=g['centroid']-s['centroid'] if g else None,response_width=g['width'] if g else None,
                            confidence=min(s['confidence'],g['confidence']) if g else None,
                            candidate_target_ids=';'.join(str(j) for j in candidates),
                            reason='unique candidate edge' if unique else ('multiple compatible sources or responses' if candidates else 'no supported response'),
                            expected_peak_delay=expected.get(i),
                            estimator_peak=(expected.get(i) if method=='v2' else g['peak']-s['peak']) if g else None,
                            estimator_onset=(g['start']-s['start']) if g and method!='v2' else None,
                            estimator_centroid=(g['centroid']-s['centroid']) if g and method!='v2' else None))
    return records


def truth_assignment(detected, truth_peaks, tolerance):
    if not len(detected) or not len(truth_peaks):return {}
    costs=np.abs(np.array([e['peak'] for e in detected])[:,None]-np.asarray(truth_peaks)[None,:])
    a,b=linear_sum_assignment(costs)
    return {int(i):int(j) for i,j in zip(a,b) if costs[i,j]<=tolerance}


def evaluate_matches(records,sources,targets,truth,test_start,cfg):
    """Truth is used for evaluation only, never candidate generation."""
    srcmap=truth_assignment(sources,[e['source_peak'] for e in truth],cfg['events']['truth_peak_tolerance_s'])
    gasids=[i for i,e in enumerate(truth) if e.get('gas_peak') is not None]
    localgas=truth_assignment(targets,[truth[i]['gas_peak'] for i in gasids],cfg['events']['truth_peak_tolerance_s'])
    gasmap={i:gasids[j] for i,j in localgas.items()}
    eligible={i for i,e in enumerate(truth) if e['source_peak']>=test_start and e.get('gas_peak') is not None}
    # All null truth records have absent gas peaks, so gasmap is empty there.
    matched=[r for r in records if r['status']=='matched' and r['source_peak']>=test_start]
    tp=0;errors={k:[] for k in ['onset','peak','centroid']};all_peak=[]
    for r in matched:
        i=srcmap.get(r['source_index']);j=gasmap.get(r['target_index'])
        if i in eligible:
            all_peak.append(abs(r['estimator_peak']-truth[i]['delay_peak']))
        if i in eligible and i==j:
            tp+=1
            for k in errors:
                if r['estimator_'+k] is not None:errors[k].append(abs(r['estimator_'+k]-truth[i]['delay_'+k]))
    fp=len(matched)-tp;fn=len(eligible)-tp
    precision=tp/(tp+fp) if tp+fp else None;recall=tp/len(eligible) if eligible else None
    f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None
    test_sources=[r for r in records if r['source_peak']>=test_start]
    detected_true=sum(i in eligible for i in srcmap.values())
    return dict(tp=tp,fp=fp,fn=fn,precision=precision,recall=recall,f1=f1,
                false_match_rate=fp/len(matched) if matched else None,
                source_detection_recall=detected_true/len(eligible) if eligible else None,
                matched_sources=len(matched),detected_sources=len(test_sources),truth_events=len(eligible),
                ambiguous_sources=sum(r['status']=='ambiguous' for r in test_sources),
                null_match_fraction=len(matched)/len(test_sources) if test_sources and not eligible else None,
                all_matched_source_peak_mae_s=float(np.mean(all_peak)) if all_peak else None,
                **{f'true_match_{k}_mae_s':float(np.mean(v)) if v else None for k,v in errors.items()})


def adjust_bh(pvalues):
    p=np.asarray(pvalues,float);order=np.argsort(p);ranked=p[order]*len(p)/np.arange(1,len(p)+1)
    values=np.minimum.accumulate(ranked[::-1])[::-1];out=np.empty_like(p);out[order]=np.clip(values,0,1)
    return out
