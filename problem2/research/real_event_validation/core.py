"""Evidence accounting with explicit overlap and no physical-label inference."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def graph(sources, targets, horizon, same_sign=True):
    edges=[]
    for i,s in enumerate(sources):
        for j,g in enumerate(targets):
            lag=g['peak']-s['peak']
            if 0 <= lag <= horizon and (not same_sign or s['polarity']==g['polarity']):
                edges.append((i,j))
    sd=np.bincount([i for i,j in edges],minlength=len(sources))
    td=np.bincount([j for i,j in edges],minlength=len(targets))
    return edges,sd,td


def intervals(sources, horizon):
    """Connected unions of [start, end+horizon], not independent events."""
    ids=[];group=-1;right=-np.inf
    for s in sources:
        if s['start']>right:group+=1
        right=max(right,s['end']+horizon)
        ids.append(group)
    return ids


def temporal_recurrence(reference, other, tolerance):
    """One-to-one same-polarity temporal agreement, not detection accuracy."""
    if not reference or not other:return set()
    cost=np.array([[abs(a['peak']-b['peak']) if a['polarity']==b['polarity'] and abs(a['peak']-b['peak'])<=tolerance else 1e9 for b in other] for a in reference])
    i,j=linear_sum_assignment(cost)
    return {int(a) for a,b in zip(i,j) if cost[a,b]<=tolerance}


def shifted_targets(targets, shift, duration):
    """Circularly shift complete event lobes; drop lobes crossing record seam."""
    out=[]
    for t in targets:
        peak=(t['peak']+shift)%duration;delta=peak-t['peak']
        s={**t,**{key:t[key]+delta for key in ['start','peak','centroid','end']}}
        if 0<=s['start'] and s['end']<duration:out.append(s)
    return sorted(out,key=lambda t:t['peak'])


def unique_count(sources,targets,horizon):
    edges,sd,td=graph(sources,targets,horizon)
    return sum(sd[i]==1 and td[j]==1 for i,j in edges)
