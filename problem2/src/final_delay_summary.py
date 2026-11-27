"""Describe existing V2 windows by flow; no re-estimation or model fitting."""
import numpy as np


def summarize(heats,quantiles=(0.,1/3,2/3,1.)):
    flows=[];records=[]
    for d in heats:
        centers=d['delay_co_center_index']
        if not np.array_equal(centers,d['delay_co2_center_index']):raise ValueError('gas candidate centers differ; do not silently double-count')
        q=d['oxygen_flow'][centers];flows.extend(q)
        for gas in ['co','co2']:
            for i in np.flatnonzero(d[f'delay_{gas}_valid_mask']):
                records.append(dict(heat_id=d['heat_id'],gas=gas,center_index=int(centers[i]),flow_m3_s=float(q[i]),delay_s=float(d[f'delay_{gas}_observed'][i])))
    bounds=np.quantile(flows,quantiles)
    if not np.all(np.diff(bounds)>0):raise ValueError('flow quantile boundaries not distinct')
    membership=lambda q,k:bounds[k]<=q and (q<bounds[k+1] or k==2 and q<=bounds[k+1])
    rows=[]
    for k,label in enumerate(['low','medium','high']):
        denominator=sum(membership(q,k) for q in flows)
        for gas in ['co','co2']:
            rr=[r for r in records if r['gas']==gas and membership(r['flow_m3_s'],k)];delays=[r['delay_s'] for r in rr]
            rows.append(dict(flow_band=label,flow_min_m3_s=float(bounds[k]),flow_max_m3_s=float(bounds[k+1]),upper_inclusive=k==2,
                gas=gas,candidate_windows=denominator,accepted_windows=len(rr),contributing_heats=len({r['heat_id'] for r in rr}),
                heat_ids=';'.join(sorted({r['heat_id'] for r in rr})),delay_median_s=float(np.median(delays)) if rr else None,
                delay_min_s=float(min(delays)) if rr else None,delay_max_s=float(max(delays)) if rr else None,
                evidence='algorithm-accepted overlapping windows; physical response identity unconfirmed' if rr else 'unidentifiable: no accepted windows',
                confirmed_physical_events=0,independent_sample_size=None))
    if sum(r['accepted_windows'] for r in rows)!=len(records):raise AssertionError('window allocation lost or duplicated')
    return rows,dict(basis='all common V2 candidate centers, counted once across gases; instantaneous center flow, not event mean',
        quantiles=list(quantiles),boundaries_m3_s=bounds.tolist(),candidate_centers=len(flows),accepted_windows=len(records),
        precision='CSV/JSON full precision determines membership; chapter thresholds are rounded for display',
        inference='descriptive only; no independence assumption, confidence interval or flow-delay regression')
