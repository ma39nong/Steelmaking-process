"""Exploratory sensitivity on full grids of three predeclared contrasting heats; no retuning final selection."""
import json,time
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from audit import P,HEATS
from models import fit
from anomalies import detect
from experiment import metrics
from eda import save

def run():
 rows=[]
 for h in ['6B03205','6C02058','6C02108']:
  d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];valid=np.load(P/f'data/masks/{h}.npz')['valid'];rng=np.random.default_rng(44551+int(h[-3:]));q=(rng.random(y.shape)<.03)&valid;qr=np.zeros(y.shape,bool);i=int(len(t)*.48);qr[i:i+4]=valid[i:i+4];hidden=q|qr;vis=y.copy();vis[hidden]=np.nan;train,_,_=detect(vis,t,w,window=5,threshold=9)
  base=dict(kt=96,kw=384,alpha_t=.01,alpha_w=.005,alpha_cross=0.,delta=1.345,irls=10)
  sweep={'alpha_t':[.001,.01,.1],'alpha_w':[.0005,.005,.05],'kt':[48,96,192],'kw':[192,384,768],'delta':[1.345,2.5,5.0],'alpha_cross':[0.,.001,.01]}
  for name,values in sweep.items():
   for value in values:
    p=dict(base);p[name]=value;start=time.perf_counter();model=fit('C',vis,t,w,train,p);pred=model.grid
    r=metrics(y,pred,q,w);rp=metrics(y,pred,qr,w,True)
    rows.append(dict(heat_id=h,parameter=name,value=value,seconds=time.perf_counter()-start,rmse_random=r['rmse'],nrmse_random=r['nrmse'],time_block_peak_height_relative=rp['peak_height_relative'],time_block_peak_position_nm=rp['peak_position_nm'],solver_converged=model.diagnostics['converged'],params=json.dumps(p)))
   print(h,name,flush=True)
 df=pd.DataFrame(rows);df.to_csv(P/'results/validation/spline_sensitivity.csv',index=False)
 for name in sweep:
  fig,axs=plt.subplots(1,2,figsize=(10,4))
  for h in df.heat_id.unique():
   sub=df[(df.heat_id==h)&(df.parameter==name)].sort_values('value');axs[0].plot(sub.value,sub.rmse_random,'o-',label=h);axs[1].plot(sub.value,1-sub.time_block_peak_height_relative,'o-',label=h)
  axs[0].set(xlabel=name,ylabel='Held random RMSE');axs[1].set(xlabel=name,ylabel='1 - peak height relative error (hidden frames)');axs[0].legend(fontsize=8);save(fig,P/f'figures/comparison/sensitivity_{name}')
 inj=pd.read_csv(P/'results/validation/anomaly_injection.csv');dev=inj[inj.split=='development'].groupby(['window','threshold'],as_index=False).f1.mean();fig,ax=plt.subplots(figsize=(7,4))
 for win,g in dev.groupby('window'):ax.plot(g.threshold,g.f1,'o-',label=f'Window {win}')
 ax.legend();ax.set(xlabel='Threshold',ylabel='Injected development macro F1');save(fig,P/'figures/comparison/anomaly_threshold_sensitivity')

def artificial_clipping():
 selection=json.loads((P/'results/final/model_selection.json').read_text());rows=[]
 for h in HEATS:
  d=np.load(P/f'data/processed/{h}.npz');yy,t0,w0=d['Y'],d['time'],d['wavelength'];n,m=yy.shape;i=n//2;j=int(np.argmax(np.median(yy,axis=0)));a=max(0,min(i-60,n-120));b=max(0,min(j-160,m-320));y=yy[a:a+120,b:b+320];t=t0[a:a+120];w=w0[b:b+320]
  original_valid=(y<65535)&np.isfinite(y);cap=float(np.quantile(y[original_valid],.7));block=np.zeros(y.shape,bool);block[45:65,130:190]=True;q=block&original_valid&(y>cap);vis=y.copy();vis[q]=cap;mask=original_valid&~q
  if not q.any():continue
  cp=dict(kt=32,kw=128,alpha_t=.01,alpha_w=.001,alpha_cross=0.,delta=2.5,irls=20,censor_level=cap)
  models={}
  for label,family,p in [('A_mask','A',dict(time_window=3,wave_window=9)),('B_mask','B',dict(rank=4,eta=.003,max_iter=40)),('C_mask','C',cp),('C_one_sided','C',dict(cp,censor_weight=.5))]:
   model=fit(family,vis,t,w,mask,p);pred=model.grid;r=metrics(y,pred,q,w)
   rows.append(dict(heat_id=h,model=label,cap=cap,artificial_count=int(q.sum()),known_reference_max=float(y[q].max()),below_cap_fraction=float(np.mean(pred[q]<cap)),**r));models[label]=pred
  np.savez_compressed(P/f'results/validation/{h}_artificial_clipping.npz',time=t,wavelength=w,reference=y,observed=vis,artificial_mask=q,**models)
  print('artificial clipping',h,int(q.sum()),flush=True)
 pd.DataFrame(rows).to_csv(P/'results/validation/artificial_clipping.csv',index=False)
if __name__=='__main__':
 if not (P/'results/final/model_selection.json').exists():raise RuntimeError('Lock family selection before independent sensitivity')
 run();artificial_clipping()
