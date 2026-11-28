"""Controlled integral-form first-order delay probe, not a new formal model."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.signal import lfilter
from .operators import pulse,acquire,fast_hampel,clean


def forward(t,oxygen,start_y,T,K,c,tau):
    dt=t[1]-t[0]
    shift=np.interp(t-tau,t,oxygen)
    q=np.diff(shift)/dt
    if T<=1e-8:return np.r_[start_y,K*q+c]
    a=np.exp(-dt/T);u=(1-a)*(K*q+c)
    out=np.empty(len(t));out[0]=start_y
    out[1:]=lfilter([1.],[1.,-a],u,zi=[a*start_y])[0]
    return out


def bounded_T_fit(X,y):
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    if beta[0]<0:
        beta=np.r_[0.,np.linalg.lstsq(X[:,1:],y,rcond=None)[0]]
    return beta


def fit(t,oxygen,y,cfg):
    p=cfg['integral_probe'];w=p['integration_interval_s'];a=np.arange(*p['training_bounds_s'],w);b=a+w
    e=cumulative_trapezoid(y,t,initial=0.)
    dy=(np.interp(b,t,y)-np.interp(a,t,y))/w
    avg=(np.interp(b,t,e)-np.interp(a,t,e))/w
    l=p['lag_grid_s'];lags=np.arange(l['min'],l['max']+l['step']/2,l['step'])
    coefficients=[];loss=[];excitation=[]
    for tau in lags:
        q=(np.interp(b-tau,t,oxygen)-np.interp(a-tau,t,oxygen))/w
        X=np.column_stack([dy,-q,-np.ones(len(q))]);coef=bounded_T_fit(X,-avg)
        coefficients.append(coef);loss.append(np.mean((X@coef+avg)**2));excitation.append(np.std(q))
    loss=np.array(loss);best=int(np.argmin(loss));T,K,c=coefficients[best];tau=lags[best]
    tolerance=p['profile_loss_tolerance_fraction_of_output_variance']*max(float(np.var(avg)),1e-12)
    near=lags[loss<=loss[best]+tolerance];width=float(near[-1]-near[0])
    # Fit the no-input nested model independently using training data only.
    X0=np.column_stack([dy,-np.ones(len(dy))]);nullcoef=np.linalg.lstsq(X0,-avg,rcond=None)[0]
    if nullcoef[0]<0:nullcoef=np.array([0.,np.mean(avg)])
    va,vb=p['validation_bounds_s'];keep=(t>=va)&(t<=vb);tt=t[keep]
    # Use full input history; validation only initializes from the first observed output.
    shifted=np.interp(tt-tau,t,oxygen);dt=t[1]-t[0];qv=np.diff(shifted)/dt
    def predict(T,K,c,q):
        out=np.empty(len(tt));out[0]=y[keep][0]
        if T<=1e-8:out[1:]=K*q+c
        else:
            alpha=np.exp(-dt/T);out[1:]=lfilter([1.],[1.,-alpha],(1-alpha)*(K*q+c),zi=[alpha*out[0]])[0]
        return out
    pred=predict(T,K,c,qv);pred0=predict(nullcoef[0],0.,nullcoef[1],np.zeros_like(qv))
    rmse=float(np.sqrt(np.mean((pred[1:]-y[keep][1:])**2)));nullrmse=float(np.sqrt(np.mean((pred0[1:]-y[keep][1:])**2)))
    improvement=1-rmse/nullrmse if nullrmse>1e-8 else 0.
    reasons=[]
    if max(excitation)<p['source_mean_flow_std_min']:reasons.append('insufficient_input_variation')
    if width>p['profile_width_max_s']:reasons.append('broad_profile')
    if abs(K)<p['absolute_gain_min']:reasons.append('near_zero_gain')
    if improvement<p['validation_improvement_min']:reasons.append('no_heldout_improvement')
    if best in [0,len(lags)-1]:reasons.append('search_boundary')
    return dict(estimate_s=float(tau),T_s=float(T),gain=float(K),offset=float(c),accepted=not reasons,
                rejection_reasons=';'.join(reasons),profile_lower_s=float(near[0]),profile_upper_s=float(near[-1]),
                profile_width_s=width,train_loss=float(loss[best]),validation_rmse=rmse,null_validation_rmse=nullrmse,
                validation_improvement=float(improvement),max_source_variation=float(max(excitation))),lags,loss


def generate(cfg,tau,T,kind,seed):
    p=cfg['integral_probe'];fine=p['fine_dt_s'];tf=np.arange(round(p['duration_s']/fine)+1)*fine
    q=np.full_like(tf,5.)
    if kind!='constant':
        rng=np.random.default_rng(seed)
        for center,width,amp in [(85,4,3),(160,16,-2),(245,30,3),(330,16,2)]:q+=pulse(tf,center,width,amp)[0]
    o=25+cumulative_trapezoid(q,tf,initial=0.)
    driver=q
    if kind=='unrelated':
        rng=np.random.default_rng(seed);driver=np.full_like(tf,5.)
        for center,width,amp in [(90,25,3),(180,40,-2),(235,8,2),(330,10,-2)]:
            driver+=pulse(tf,center+rng.uniform(-20,20),width*rng.uniform(.7,1.3),amp)[0]
    delayed=np.interp(tf-tau,tf,driver,left=driver[0]);K=p['true_gain'];c=p['true_offset']
    stages=2 if kind=='two_stage' else 1;ys=K*delayed+c
    for _ in range(stages):
        a=np.exp(-fine/(T/stages));ys=np.r_[ys[0],lfilter([1-a],[1.,-a],ys[:-1],zi=[a*ys[0]])[0]]
    stride=round(cfg['dt_s']/fine)
    return tf[::stride],o[::stride],ys[::stride]


def execute(cfg):
    p=cfg['integral_probe'];rows=[];profiles=[];examples=[]
    settings=[('matched',d,T) for d in p['true_delays_s'] for T in p['true_time_constants_s']]
    settings += [('constant',15,8),('unrelated',15,8),('two_stage',15,8)]
    for kind,tau,T in settings:
        for noise in p['noise_std']:
            seeds=[cfg['seeds'][0]] if noise==0 else cfg['seeds']
            for seed in seeds:
                t,o,ytrue=generate(cfg,tau,T,kind,seed);rng=np.random.default_rng(seed)
                noisy=ytrue+rng.normal(0,noise,len(t))
                for mechanism in p['source_acquisition']:
                    measured=acquire(o,t,mechanism)
                    for gas_input in p['gas_inputs']:
                        yy=noisy if gas_input=='raw' else clean(t,noisy,fast_hampel(noisy))
                        r,grid,loss=fit(t,measured,yy,cfg)
                        truth=tau if kind not in ['constant','unrelated'] else None
                        labels=dict(kind=kind,true_delay_s=truth,true_T_s=T,noise_std=noise,seed=seed,
                                    acquisition=mechanism,gas_input=gas_input)
                        rows.append(dict(**labels,**r,error_s=r['estimate_s']-truth if truth is not None else None,
                                         model_kernel_centroid_error_s=r['estimate_s']+r['T_s']-truth-T if truth is not None else None,
                                         correct=bool(r['accepted'] and truth is not None and abs(r['estimate_s']-truth)<=p['correct_delay_tolerance_s'])))
                        if noise==0 and mechanism in ['ideal','hold1_phase0'] and gas_input=='raw':
                            for g,l in zip(grid,loss):profiles.append(dict(**labels,lag_s=float(g),loss=float(l)))
            print('integral probe',kind,tau,T,noise,flush=True)
    return rows,profiles
