"""Fidelity of known planted components, not labels for unknown industrial events."""
import numpy as np
from ..measurement_chain.core import crossing


def metrics(t,response,center,width,amplitude,shape,cfg):
    from .operators import pulse
    truth,integral=pulse(t,center,width,amplitude,shape)
    area=amplitude*width/(2 if shape=='triangle' else 1)
    padding=cfg['roi_padding_s'];keep=abs(t-center)<=width/2+padding
    tt=t[keep];r=np.asarray(response)[keep];expected=truth[keep];positive=np.maximum(r,0.)
    peak=max(float(np.max(r)),0.);threshold=cfg['onset_true_peak_fraction']*amplitude
    onset=crossing(tt,r,threshold) if peak>=threshold else None
    trueonset=center-width/2+(cfg['onset_true_peak_fraction']*width/2 if shape=='triangle' else 0.)
    ownonset=crossing(tt,positive,.1*peak) if peak>1e-12 else None
    parea=float(np.trapezoid(positive,tt))
    centroid=float(np.trapezoid(tt*positive,tt)/parea) if parea>1e-12 else None
    if peak>1e-12:
        maxima=tt[np.isclose(r,peak,rtol=1e-9,atol=1e-10)]
        if shape=='triangle':peakdist=float(min(abs(maxima-center)))
        else:peakdist=float(np.min(np.maximum(np.maximum(center-width/2-maxima,maxima-(center+width/2)),0)))
        lo,hi=crossing(tt,positive,.5*peak),crossing(tt,positive,.5*peak,True)
        apparent_width=hi-lo if lo is not None and hi is not None else None
    else:maxima=[];peakdist=apparent_width=None
    far=abs(t-center)>width/2+10
    return dict(peak_ratio=peak/amplitude,peak_at_least_half=bool(peak>=.5*amplitude),
                absolute_onset_defined=onset is not None,onset_absolute_error_s=onset-trueonset if onset is not None else None,
                onset_relative_error_s=ownonset-trueonset if ownonset is not None else None,
                centroid_error_s=centroid-center if centroid is not None else None,
                nearest_peak_set_distance_s=peakdist,estimated_peak_first_s=float(maxima[0]) if len(maxima) else None,
                estimated_peak_last_s=float(maxima[-1]) if len(maxima) else None,
                apparent_width50_error_s=apparent_width-width/(2 if shape=='triangle' else 1) if apparent_width is not None else None,
                signed_area_relative_error=float(np.trapezoid(response,t)/area-1),
                positive_area_ratio=float(np.trapezoid(np.maximum(response,0),t)/area),
                negative_area_ratio=float(np.trapezoid(np.maximum(-response,0),t)/area),
                far_leakage_absolute_area_ratio=float(np.trapezoid(np.where(far,abs(response),0),t)/area),
                waveform_nrmse=float(np.sqrt(np.mean((r-expected)**2))/amplitude))


def cumulative_metrics(t,delta,truth_integral,center,width,volume):
    keep=abs(t-center)<=width/2+15
    return dict(final_volume_relative_error=float(delta[-1]/volume-1),
                cumulative_nrmse=float(np.sqrt(np.mean((delta[keep]-truth_integral[keep])**2))/volume),
                maximum_cumulative_error_relative=float(np.max(abs(delta-truth_integral))/volume))
