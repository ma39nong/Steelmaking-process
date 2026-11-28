"""Same first-order model; ordinary versus compact smooth weighted moments.

A = -integral(phi' y)/integral(phi), B = integral(phi y)/integral(phi),
Q(tau) = -integral(phi' O(t-tau))/integral(phi).
T*A + B = K*Q + c. phi=sin(pi*(t-a)/(b-a))**2 vanishes at endpoints.
All moments stay inside training bounds, for every lag. Overlapping rows are
never treated as independent noise samples in the full-record bootstrap.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.signal import lfilter


class Design:
    def __init__(self,t,o,cfg,method='weak',support=None):
        self.t=np.asarray(t);self.o=np.asarray(o);self.cfg=cfg;self.method=method
        p=cfg['integral_probe'];lo,hi=p['training_bounds_s'];dt=t[1]-t[0]
        self.lags=np.arange(p['lag_grid_s']['min'],p['lag_grid_s']['max']+.25,p['lag_grid_s']['step'])
        n=len(t)
        if lo-self.lags.max()<t[0] or hi>t[-1]:raise ValueError('insufficient training history')
        if method=='ordinary':
            w=p['integration_interval_s'];a=np.arange(lo,hi,w);b=a+w
            # Matches the original endpoint/interpolated cumulative-trapezoid definition.
            self.D=np.zeros((len(a),n));self.W=np.zeros_like(self.D)
            for j,(aa,bb) in enumerate(zip(a,b)):
                ia=round((aa-t[0])/dt);ib=round((bb-t[0])/dt)
                if not np.isclose(t[ia],aa) or not np.isclose(t[ib],bb):raise ValueError('unaligned ordinary bounds')
                self.D[j,ia]=-1/w;self.D[j,ib]=1/w
                self.W[j,ia:ib+1]=dt/w;self.W[j,[ia,ib]]=dt/(2*w)
            self.Q=np.array([(np.interp(b-l,t,o)-np.interp(a-l,t,o))/w for l in self.lags])
        else:
            w=support or cfg['weak_support_s'];a=np.arange(lo,hi-w+1e-8,cfg['weak_stride_s'])
            self.D=np.zeros((len(a),n));self.W=np.zeros_like(self.D)
            for j,aa in enumerate(a):
                idx=np.flatnonzero((t>=aa-1e-8)&(t<=aa+w+1e-8));s=(t[idx]-aa)/w
                phi=np.sin(np.pi*s)**2;der=np.pi/w*np.sin(2*np.pi*s)
                quad=np.full(len(idx),dt);quad[[0,-1]]*=.5;norm=np.sum(phi*quad)
                self.W[j,idx]=phi*quad/norm;self.D[j,idx]=-der*quad/norm
                # Remove roundoff-only constant leakage; essential for large cumulative offsets.
                self.D[j,idx]-=self.D[j,idx].sum()/len(idx)
            self.Q=np.array([self.D@(np.interp(t-l,t,o)-o[0]) for l in self.lags])
        self.max_excitation=float(np.std(self.Q,axis=1).max())

    def profile(self,y):
        d=self.D@y;b=self.W@y
        X=np.empty((len(self.lags),len(d),3));X[:,:,0]=d;X[:,:,1]=-self.Q;X[:,:,2]=-1.
        # Scale columns, use SVD pseudoinverses: rank deficiency is not regularized away.
        scales=np.sqrt(np.mean(X**2,axis=1));scales=np.maximum(scales,1e-12)
        z=X/scales[:,None,:];sv=np.linalg.svd(z,compute_uv=False)
        co=np.einsum('lij,j->li',np.linalg.pinv(z,rcond=1e-12),-b)/scales
        neg=co[:,0]<0
        if np.any(neg):
            co[neg,0]=0
            co[neg,1:]=np.einsum('lij,j->li',np.linalg.pinv(X[neg,:,1:],rcond=1e-12),-b)
        res=np.einsum('lij,lj->li',X,co)+b
        loss=np.mean(res**2,axis=1)
        condition=sv[:,0]/np.maximum(sv[:,-1],1e-15)
        return co,loss,condition,d,b

    def fit(self,y):
        co,loss,cond,d,b=self.profile(y);j=int(np.argmin(loss));T,K,c=co[j];tau=self.lags[j]
        p=self.cfg['integral_probe'];tol=p['profile_loss_tolerance_fraction_of_output_variance']*max(float(np.var(b)),1e-12)
        near=self.lags[loss<=loss[j]+tol];width=float(np.ptp(near))
        x0=np.column_stack([d,-np.ones(len(d))]);null=np.linalg.lstsq(x0,-b,rcond=None)[0]
        if null[0]<0:null=np.array([0.,float(np.mean(b))])
        va,vb=p['validation_bounds_s'];idx=np.flatnonzero((self.t>=va)&(self.t<=vb));tt=self.t[idx]
        pred=predict(tt,self.t,self.o,y[idx[0]],T,K,c,tau)
        pred0=predict(tt,self.t,self.o,y[idx[0]],null[0],0.,null[1],0.)
        rmse=float(np.sqrt(np.mean((pred[1:]-y[idx[1:]])**2)));nr=float(np.sqrt(np.mean((pred0[1:]-y[idx[1:]])**2)))
        improvement=1-rmse/nr if nr>1e-8 else 0.
        reasons=[]
        if self.max_excitation<p['source_mean_flow_std_min']:reasons.append('insufficient_input_variation')
        if width>p['profile_width_max_s']:reasons.append('broad_profile')
        if abs(K)<p['absolute_gain_min']:reasons.append('near_zero_gain')
        if improvement<p['validation_improvement_min']:reasons.append('no_heldout_improvement')
        if j in (0,len(loss)-1):reasons.append('search_boundary')
        if cond[j]>1e8:reasons.append('rank_deficient_moments')
        return dict(estimate_s=float(tau),T_s=float(T),gain=float(K),offset=float(c),
                    accepted=not reasons,rejection_reasons=';'.join(reasons),train_loss=float(loss[j]),
                    profile_lower_s=float(near[0]),profile_upper_s=float(near[-1]),profile_width_s=width,
                    validation_rmse=rmse,null_validation_rmse=nr,validation_improvement=improvement,
                    max_source_variation=self.max_excitation,scaled_design_condition=float(cond[j]),
                    null_T_s=float(null[0]),null_offset=float(null[1])),co,loss,cond


def predict(tt,history_t,oxygen,initial,T,K,c,tau):
    dt=tt[1]-tt[0];q=np.diff(np.interp(tt-tau,history_t,oxygen))/dt
    out=np.empty(len(tt));out[0]=initial
    if T<=1e-8:out[1:]=K*q+c
    else:
        a=np.exp(-dt/T);out[1:]=lfilter([1.],[1.,-a],(1-a)*(K*q+c),zi=[a*initial])[0]
    return out


def bootstrap(design,y,fit,cfg,seed,block=False):
    """Re-fit whole noisy series, NOT independent overlapped moment rows.

    Synthetic: conditional Gaussian model bootstrap with sigma estimated from
    training residual. Real: moving-block centered residual bootstrap. Neither
    can validate a misspecified pure transport delay; coverage is tested outside.
    """
    t=design.t;p=cfg['integral_probe'];a,b=p['training_bounds_s'];train=(t>=a)&(t<=b)
    mean=predict(t,t,design.o,y[0],fit['T_s'],fit['gain'],fit['offset'],fit['estimate_s'])
    resid=y[train]-mean[train];resid=resid-np.mean(resid)
    sigma=float(np.std(resid,ddof=3));rng=np.random.default_rng(seed);samples=[]
    for _ in range(cfg['bootstrap_replicates']):
        if block:
            length=round(cfg['real_data_if_pass']['moving_block_s']/(t[1]-t[0]));starts=rng.integers(0,len(resid),int(np.ceil(len(t)/length)))
            noise=np.concatenate([resid[(j+np.arange(length))%len(resid)] for j in starts])[:len(t)]
        else:noise=rng.normal(0,sigma,len(t))
        co,loss,_,_,_=design.profile(mean+noise);j=int(np.argmin(loss))
        samples.append([design.lags[j],*co[j]])
    samples=np.asarray(samples);alpha=(1-cfg['bootstrap_confidence'])/2
    lower,upper=np.quantile(samples[:,0],[alpha,1-alpha])
    corr=float(np.corrcoef(samples[:,:2].T)[0,1]) if np.std(samples[:,0])>0 and np.std(samples[:,1])>0 else None
    return dict(ci_lower_s=float(lower),ci_upper_s=float(upper),ci_width_s=float(upper-lower),
                bootstrap_noise_sigma=sigma,bootstrap_tau_T_correlation=corr,
                bootstrap_search_boundary_fraction=float(np.mean((samples[:,0]==design.lags[0])|(samples[:,0]==design.lags[-1])))),samples
