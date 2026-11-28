"""Exact baseline operators and one-component experimental ablations."""
import warnings
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter
from ...src.q2_pipeline import hampel_mask,pchip_clean


def fast_hampel(x,window=9,threshold=6.):
    """Vectorized evaluation of the unchanged production rule, including edges."""
    x=np.asarray(x,float);h=window//2
    finite=np.where(np.isfinite(x),x,np.nan)
    a=sliding_window_view(np.pad(finite,(h,h),constant_values=np.nan),window)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        med=np.nanmedian(a,axis=1);scale=1.4826*np.nanmedian(abs(a-med[:,None]),axis=1)
        count=np.sum(np.isfinite(a),axis=1);mask=~np.isfinite(x)|(x<0)
        mask|=(count>=5)&(scale>0)&(abs(x-med)>threshold*scale)
        left,right=a[:,:h],a[:,h+1:]
        lm,rm=np.nanmedian(left,axis=1),np.nanmedian(right,axis=1)
        nb=np.concatenate([left,right],axis=1);nm=np.nanmedian(nb,axis=1)
        ns=1.4826*np.nanmedian(abs(nb-nm[:,None]),axis=1)
        candidate=(count>=5)&(scale==0)&(np.sum(np.isfinite(left),axis=1)>=2)&(np.sum(np.isfinite(right),axis=1)>=2)
        mask|=candidate&(abs(lm-rm)<=np.maximum(1e-12,3*ns))&(abs(x-nm)>threshold*np.maximum(ns,1e-12))
    return mask


def pulse(t,center,width,amplitude=1.,shape='triangle'):
    s=np.asarray(t)-center+width/2
    if shape=='triangle':
        q=amplitude*np.maximum(1-abs(np.asarray(t)-center)/(width/2),0)
        z=np.clip(s,0,width)
        cumulative=np.where(z<=width/2,amplitude*z*z/width,
                            amplitude*width/2-amplitude*(width-z)**2/width)
    elif shape=='rectangle':
        q=amplitude*((s>=0)&(s<width)).astype(float)
        cumulative=amplitude*np.clip(s,0,width)
    else:raise ValueError(shape)
    return q,cumulative


def acquire(x,t,name,schedule=None):
    """A hypothetical acquisition operator applied to whole baseline and augmented signals."""
    x=np.asarray(x,float)
    if name=='ideal':return x.copy()
    if name=='grid_0.01':return np.rint(x/.01)*.01
    if name=='empirical_hold':
        if schedule is None:raise ValueError('empirical schedule required')
        changed=np.r_[True,np.diff(np.asarray(schedule))!=0]
        idx=np.maximum.accumulate(np.where(changed,np.arange(len(x)),0))
    elif name.startswith('hold'):
        period,phase=name[4:].split('_phase');period=float(period);phase=float(phase)
        dt=float(np.median(np.diff(t)));p=round(period/dt);ph=round(phase/dt)
        if p<1 or not 0<=ph<p:raise ValueError(name)
        idx=np.maximum(ph+((np.arange(len(x))-ph)//p)*p,0)
    else:raise ValueError(name)
    return np.rint(x[idx]/.01)*.01


def clean(t,x,mask,linear=False):
    if linear:
        good=(~mask)&np.isfinite(x)
        if sum(good)<2:raise ValueError('not enough visible samples')
        return np.maximum(np.interp(t,t[good],x[good]),0.)
    return pchip_clean(t,x,mask)


def source_methods(t,x,names):
    m=fast_hampel(x);y=clean(t,x,m);dt=float(np.median(np.diff(t)))
    sg31=savgol_filter(y,31,3,deriv=1,delta=dt,mode='interp');result={}
    for name in names:
        yy=y;mm=m
        if name=='formal':q=sg31
        elif name=='omit_oxygen_hampel':
            mm=~np.isfinite(x)|(x<0);yy=clean(t,x,mm);q=savgol_filter(yy,31,3,deriv=1,delta=dt,mode='interp')
        elif name=='linear_oxygen_interpolation':
            yy=clean(t,x,m,True);q=savgol_filter(yy,31,3,deriv=1,delta=dt,mode='interp')
        elif name=='SG11':q=savgol_filter(y,11,3,deriv=1,delta=dt,mode='interp')
        elif name=='analytic_PCHIP_derivative':
            good=(~m)&np.isfinite(x);p=PchipInterpolator(t[good],x[good],extrapolate=False)
            q=p.derivative()(t);q=np.where(np.isfinite(q),q,0.)
        elif name=='secant_2s':
            lo=np.maximum(t[0],t-1);hi=np.minimum(t[-1],t+1)
            q=(np.interp(hi,t,y)-np.interp(lo,t,y))/(hi-lo)
        elif name=='omit_flow_clipping':q=sg31
        else:raise ValueError(name)
        if name!='omit_flow_clipping':q=np.maximum(q,0.)
        result[name]=dict(flow=q,clean=yy,mask=mm)
    return result


def gas_methods(t,x,names):
    m=fast_hampel(x);result={}
    for name in names:
        if name=='formal':mm=m
        elif name=='omit_gas_hampel':mm=~np.isfinite(x)|(x<0)
        elif name=='hampel_threshold12':mm=fast_hampel(x,threshold=12.)
        elif name=='hampel_window5':mm=fast_hampel(x,window=5)
        elif name=='linear_gas_interpolation':mm=m
        else:raise ValueError(name)
        result[name]=dict(clean=clean(t,x,mm,name=='linear_gas_interpolation'),mask=mm)
    return result


def changes(x,name):
    d=np.diff(x,prepend=x[0]);series=pd.Series(d)
    if name=='formal_median21':return d-series.rolling(21,center=True,min_periods=1).median().to_numpy()
    if name=='difference_only':return d
    if name=='mean21':return d-series.rolling(21,center=True,min_periods=1).mean().to_numpy()
    raise ValueError(name)
