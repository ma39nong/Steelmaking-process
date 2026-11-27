"""Event-labelled additions to the prior audit, not another pure-shift suite."""
import numpy as np
from .core import changes,descriptors


def response_kernel(kind,dt,L):
    t=np.arange(int(round(L/dt))+1)*dt
    if kind=='gaussian_broadening':h=np.exp(-.5*((t-18)/5)**2)
    elif kind=='slow_exponential':h=np.where(t>=7,np.exp(-np.maximum(t-7,0)/15),0)
    else:
        z=np.maximum(t-5,0);h=z*np.exp(-z/5)
    h=h/h.sum()
    if kind=='biphasic_response':
        second=np.exp(-.5*((t-38)/4)**2);second/=second.sum();h=h-.8*second
    return h


def event_case(kind,seed,cfg):
    rng=np.random.default_rng(seed);dt=cfg['dt_s'];n=cfg['synthetic']['samples'];t=np.arange(n)*dt
    times=np.array([80,170,260,350,440,530,620,710,800.])
    if kind=='adjacent_events':times=np.array([80,95,260,275,440,455,710,725,800,815.])
    signs=np.ones(len(times)) if kind=='adjacent_events' else np.array([1 if i%3 else -1 for i in range(len(times))]);amps=rng.uniform(.8,1.2,len(times))
    if kind=='mixed_strength':amps*=np.resize(np.array([1,.15,.35]),len(times))
    h=response_kernel(kind,dt,cfg['max_response_s']);source_parts=[];gas_parts=[];truth=[]
    for i,(center,sign,a) in enumerate(zip(times,signs,amps)):
        u=sign*a*np.exp(-.5*((t-center)/2)**2)
        hi=h
        if kind=='adjacent_events':hi=response_kernel('gaussian_broadening' if i%2 else 'slow_exponential',dt,cfg['max_response_s'])
        g=np.convolve(u,hi)[:n]
        qpart=np.cumsum(u)*dt;gpart=np.cumsum(g)*dt
        # Exact noiseless descriptors AFTER the same offline measurement derivative.
        sp=changes(qpart,cfg);gp=changes(gpart,cfg)
        sd=descriptors(sp,dt);gd=descriptors(gp,dt)
        source_parts.append(u);gas_parts.append(g)
        truth.append(dict(true_event_id=i,source_peak=sd['peak'],source_start=sd['onset'],
                          gas_peak=gd['peak'] if kind!='null_response' else None,
                          gas_start=gd['onset'] if kind!='null_response' else None,
                          delay_onset=gd['onset']-sd['onset'],delay_peak=gd['peak']-sd['peak'],
                          delay_centroid=gd['centroid']-sd['centroid'],response_width=gd['width'],strength=a))
    q=20+np.cumsum(np.sum(source_parts,axis=0))*dt
    gas=20+np.cumsum(np.sum(gas_parts,axis=0))*dt
    q+=rng.normal(0,.005,n);gas+=rng.normal(0,.005,n)
    if kind=='noise_trend':
        q+=rng.normal(0,.08,n)+.002*t+.2*np.sin(t/55)
        gas+=rng.normal(0,.08,n)-.001*t+.4*np.sin(t/43)
    if kind=='null_response':
        independent=np.convolve(rng.normal(size=n),np.ones(41)/41,mode='same')
        gas=20+np.cumsum(independent)*dt+rng.normal(0,.005,n)
    return dict(time=t,flow=q,gas=gas,source=changes(q,cfg),response=changes(gas,cfg),truth=truth,
                kernel_truth=h if kind not in ['adjacent_events','null_response'] else None,
                stationary_kernel=kind not in ['adjacent_events','null_response'])
