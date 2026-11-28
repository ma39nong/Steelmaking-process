"""Timing definitions are explicit; none is an industrial ground-truth label."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from ...src.q2_pipeline import hampel_mask, pchip_clean, oxygen_flow


def integrate(t, q):
    return cumulative_trapezoid(q, t, initial=0.)


def reconstruct(t, oxygen, gas=None):
    om = hampel_mask(oxygen)
    oc = pchip_clean(t, oxygen, om)
    q = oxygen_flow(t, oc)
    if gas is None:
        return q, oc, om
    gm = hampel_mask(gas)
    return q, pchip_clean(t, gas, gm), om, gm


def pulse(t, center, duration, shape='triangle'):
    z = np.abs(t-center)/(duration/2)
    if shape == 'triangle':
        return np.maximum(1-z, 0.)
    if shape == 'rectangle':
        return ((t >= center-duration/2) & (t < center+duration/2)).astype(float)
    raise ValueError(shape)


def held(x, period_samples, phase_samples=0):
    """Right-continuous last-update hold on a uniform grid; no wraparound."""
    if period_samples < 1 or not 0 <= phase_samples < period_samples:
        raise ValueError('invalid hold period or phase')
    i = np.arange(len(x))
    update = np.maximum(phase_samples + ((i-phase_samples)//period_samples)*period_samples, 0)
    return np.asarray(x)[update]


def crossing(t, x, threshold, falling=False):
    ix = np.flatnonzero((x[:-1] >= threshold) & (x[1:] < threshold) if falling
                        else (x[:-1] < threshold) & (x[1:] >= threshold))
    if not len(ix):
        return None
    i = ix[-1] if falling else ix[0]
    return float(t[i] + (threshold-x[i])*(t[i+1]-t[i])/(x[i+1]-x[i]))


def features(t, signal):
    """Positive-part, amplitude-relative features; earliest sampled max is explicit."""
    x = np.maximum(np.asarray(signal), 0.)
    area = float(np.trapezoid(x, t))
    if not area > 1e-12:
        return dict(onset10_s=None, rise50_s=None, fall50_s=None, peak_s=None,
                    centroid_s=None, width50_s=None, area=area, peak_amplitude=float(x.max()))
    peak = float(x.max())
    rise, fall = crossing(t,x,.5*peak), crossing(t,x,.5*peak,True)
    return dict(onset10_s=crossing(t,x,.1*peak),rise50_s=rise,fall50_s=fall,
                peak_s=float(t[np.argmax(x)]),centroid_s=float(np.trapezoid(x*t,t)/area),
                width50_s=fall-rise if rise is not None and fall is not None else None,
                area=area,peak_amplitude=peak)


def errors(t, truth, estimate, prefix=''):
    a,b = features(t,truth),features(t,estimate)
    r={prefix+k+'_error': b[k]-a[k] if a[k] is not None and b[k] is not None else None
       for k in ['onset10_s','rise50_s','fall50_s','centroid_s','width50_s']}
    # A rectangular input has a set of equally valid peak times. Do not invent one.
    valid_peak = t[np.isclose(truth, np.max(truth),rtol=1e-10,atol=1e-12)]
    r[prefix+'peak_set_distance_s'] = float(np.min(abs(valid_peak-b['peak_s']))) if b['peak_s'] is not None else None
    r[prefix+'area_relative_error']=(b['area']-a['area'])/a['area']
    r[prefix+'peak_ratio']=b['peak_amplitude']/a['peak_amplitude']
    return r


def delay_features(t, source, target, truth_delay, unique_peak=True):
    a,b=features(t,source),features(t,target)
    result={}
    for name,key in [('onset','onset10_s'),('centroid','centroid_s'),('peak','peak_s')]:
        identified=a[key] is not None and b[key] is not None and (name!='peak' or unique_peak)
        result['feature_delay_'+name+'_s']=b[key]-a[key] if identified else None
        result['feature_delay_'+name+'_error_s']=b[key]-a[key]-truth_delay if identified else None
    result['feature_peak_definition']='unique sampled triangle peak' if unique_peak else 'nonunique rectangular maximum; no single peak-delay truth'
    return result


def runs(x):
    """Constant runs; adjacent NaNs never imply a held valid measurement."""
    x=np.asarray(x)
    cuts=np.r_[0,np.flatnonzero(~(np.isfinite(x[1:]) & np.isfinite(x[:-1]) & (x[1:]==x[:-1])))+1,len(x)]
    return [(int(a),int(b)) for a,b in zip(cuts[:-1],cuts[1:])]


def lag_row(source,target,center,truth,cfg,**labels):
    # The verified exact-A implementation avoids running every center; B is unused.
    from ..target_excitation.core import compare
    gates=dict(dt_s=cfg['dt_s'],window=121,max_lag_samples=120,stride=5,
               min_corr=.35,min_separation=.05,min_strength_ratio=.05,peak_distance_samples=5)
    r=compare(source,target,gates,[center])['A']
    accepted=bool(r['accepted'][0]);estimate=float(r['observed'][0]) if accepted else None
    return dict(**labels,accepted=accepted,truth_s=truth,estimate_s=estimate,
                correct=bool(accepted and truth is not None and abs(estimate-truth)<=cfg['correct_error_s']),
                wrong=bool(accepted and truth is not None and abs(estimate-truth)>cfg['correct_error_s']),
                error_s=estimate-truth if accepted and truth is not None else None,
                max_corr=float(r['score'][0]),separation=float(r['separation'][0]))
