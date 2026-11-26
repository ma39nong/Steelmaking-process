"""Nested visibility splits. All model inputs hide development and outer test observations."""
from pathlib import Path
import sys,json,time,logging,resource,platform
import numpy as np,pandas as pd
from audit import ROOT,P,HEATS
from anomalies import detect,inject,detection_metrics,fill_training,spectrum_candidates,channel_candidates
from models import fit,load_surface
from eda import save,heatmap
import matplotlib.pyplot as plt
CONFIG=json.loads((P/'config/experiment.json').read_text())

def metrics(y,pred,q,w,shape_metrics=False):
 err=pred[q]-y[q]
 if not len(err):return dict(n=0,mae=np.nan,rmse=np.nan,nrmse=np.nan,sam_deg=np.nan,peak_position_nm=np.nan,peak_height_relative=np.nan,area_relative=np.nan)
 spread=np.quantile(y[q],.95)-np.quantile(y[q],.05)
 dot=np.sum(np.where(q,y*pred,0),axis=1);n1=np.sqrt(np.sum(np.where(q,y*y,0),axis=1));n2=np.sqrt(np.sum(np.where(q,pred*pred,0),axis=1));eligible=q.sum(1)>=3
 sam=np.degrees(np.arccos(np.clip(dot[eligible]/np.maximum(n1[eligible]*n2[eligible],1e-12),-1,1)))
 result=dict(n=int(q.sum()),mae=float(np.mean(abs(err))),rmse=float(np.sqrt(np.mean(err**2))),nrmse=float(np.sqrt(np.mean(err**2))/max(spread,1)),sam_deg=float(np.mean(sam)) if len(sam) else np.nan,peak_position_nm=np.nan,peak_height_relative=np.nan,area_relative=np.nan,negative_prediction_fraction=float(np.mean(pred[q]<0)))
 if shape_metrics:
  pe=[];he=[];ar=[]
  for i in np.flatnonzero(q.sum(1)>=5):
   ix=np.flatnonzero(q[i]);a=y[i,ix];b=pred[i,ix];ww=w[ix];pe.append(abs(ww[np.argmax(a)]-ww[np.argmax(b)]));he.append(abs(b.max()-a.max())/max(abs(a.max()),1));ar.append(abs(np.trapezoid(b,ww)-np.trapezoid(a,ww))/max(abs(np.trapezoid(a,ww)),1))
  if pe:result.update(peak_position_nm=float(np.mean(pe)),peak_height_relative=float(np.mean(he)),area_relative=float(np.mean(ar)))
 return result

def split_masks(y,eligible,rng):
 n,m=y.shape;tests={}
 # The three random tests are separate experiments and can overlap each other.
 for frac in CONFIG['test_random_fractions']:tests[f'random_{frac:g}']=(rng.random(y.shape)<frac)&eligible
 tb=np.zeros(y.shape,bool);i=int(.68*n);tb[i:i+CONFIG['test_time_block_rows']]=True;tests['time_block']=tb&eligible
 wb=np.zeros(y.shape,bool);j=int(.68*m);wb[:,j:j+CONFIG['test_wavelength_block_channels']]=True;tests['wavelength_block']=wb&eligible
 union=np.logical_or.reduce(list(tests.values()))|tb|wb
 dev={};dev['random']=(rng.random(y.shape)<.02)&eligible&~union
 a=np.zeros(y.shape,bool);i=int(.31*n);a[i:i+4]=True;dev['time_block']=a&eligible&~union
 a=np.zeros(y.shape,bool);j=int(.31*m);a[:,j:j+12]=True;dev['wavelength_block']=a&eligible&~union
 all_dev=np.logical_or.reduce(list(dev.values()))
 a=np.zeros(y.shape,bool);i=int(.31*n);a[i:i+4]=True;all_dev|=a
 a=np.zeros(y.shape,bool);j=int(.31*m);a[:,j:j+12]=True;all_dev|=a
 return dev,tests,all_dev,union

def tune_detector(y,t,w,h):
 # Separate synthetic tuning/test seeds and coordinates, retaining actual physical axes.
 ix=np.arange(len(t)//3,len(t)//3+160);jx=np.arange(len(w)//4,len(w)//4+512)
 patch=y[np.ix_(ix,jx)];tt=t[ix];ww=w[jx]
 rows=[]
 for ci,cfg in enumerate(CONFIG['detector_candidates']):
  rng=np.random.default_rng(CONFIG['seed']+int(h[-3:])+31)
  clean,_,_=detect(patch,tt,ww,**cfg)
  for kind in ['spike','drop','short_block','row','column','saturation']:
   corrupt,truth=inject(patch,kind,rng);mask,types,_=detect(corrupt,tt,ww,**cfg)
   score=detection_metrics(~mask,truth,clean|truth)
   rows.append(dict(heat_id=h,split='development',candidate=ci,kind=kind,**cfg,**score))
 score=pd.DataFrame(rows).groupby('candidate').f1.mean();best=int(score.idxmax());cfg=CONFIG['detector_candidates'][best]
 rng=np.random.default_rng(CONFIG['seed']+int(h[-3:])+371)
 clean,_,_=detect(patch,tt,ww,**cfg)
 for kind in ['spike','drop','short_block','row','column','saturation']:
  corrupt,truth=inject(patch,kind,rng);mask,types,_=detect(corrupt,tt,ww,**cfg);rows.append(dict(heat_id=h,split='test',candidate=best,kind=kind,**cfg,**detection_metrics(~mask,truth,clean|truth)))
 return cfg,rows

def main():
 if not (ROOT/'references/metadata/LITERATURE_REVIEW_COMPLETED').exists():raise RuntimeError('Literature review required before models')
 if not (P/'results/eda/COMPLETE').exists():raise RuntimeError('Audit and EDA required before models')
 from provenance import check_or_stamp
 check_or_stamp()
 logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s',handlers=[logging.StreamHandler(),logging.FileHandler(P/'results/validation/experiment.log')])
 versions={m:__import__(m).__version__ for m in ['numpy','scipy','pandas','matplotlib','sklearn','pyxlsb','pybaselines']};versions['python']=sys.version;versions['platform']=platform.platform();(P/'results/validation/environment.json').write_text(json.dumps(versions,indent=2))
 tuned={};devrows=[];injection=[];anomalyrows=[]
 # Cache per-heat finished work so interrupted runs can resume without silently rerunning results.
 for h in HEATS:
  cache=P/f'results/validation/{h}_development.json'
  if cache.exists():
   d=json.loads(cache.read_text());tuned[h]=d['tuned'];devrows.extend(d['development']);injection.extend(d['injection']);anomalyrows.append(d['anomaly']);logging.info('Resume development %s',h);continue
  d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];logging.info('Development %s',h)
  detector,inj=tune_detector(y,t,w,h);injection.extend(inj)
  valid,types,scores=detect(y,t,w,**detector)
  distance,frames=spectrum_candidates(fill_training(np.where(valid,y,np.nan),t,w));cs,cols=channel_candidates(y)
  types[frames]|=8;types[:,cols]|=16
  np.savez_compressed(P/f'data/masks/{h}.npz',valid=valid,types=types,frame_candidate=frames,channel_candidate=cols,time=t,wavelength=w)
  np.savez_compressed(P/f'data/processed/{h}_cleaned.npz',Y=np.where(valid,y,np.nan),time=t,wavelength=w)
  pd.DataFrame(dict(time_s=t,robust_mahalanobis_squared=distance,candidate=frames)).to_csv(P/f'results/anomaly_detection/{h}_frames.csv',index=False)
  pd.DataFrame(dict(wavelength_nm=w,score=cs,candidate=cols)).to_csv(P/f'results/anomaly_detection/{h}_channels.csv',index=False)
  ii,jj=np.where(~valid);pd.DataFrame(dict(time_s=t[ii],wavelength_nm=w[jj],type=types[ii,jj],value=y[ii,jj])).to_csv(P/f'results/anomaly_detection/{h}_positions.csv',index=False)
  heatmap((~valid).astype(float),t,w,P/f'figures/comparison/{h}_anomaly_mask',h+' reconstruction mask',vmax=1)
  heatmap(np.where(valid,y,np.nan),t,w,P/f'figures/cleaned/{h}_heatmap',h+' cleaned (masked)')
  ar=dict(heat_id=h,local_count=int(((types&4)>0).sum()),saturated_count=int(((types&2)>0).sum()),masked_count=int((~valid).sum()),masked_fraction=float((~valid).mean()),frame_candidates=int(frames.sum()),channel_candidates=int(cols.sum()),noise_floor=scores['noise_floor']);anomalyrows.append(ar)
  # Eligibility is a fixed reference-quality definition, not fitting data.
  eligible=valid & (y<64000)
  rng=np.random.default_rng(CONFIG['seed']+int(h[-3:]));dev,tests,du,tu=split_masks(y,eligible,rng)
  np.savez_compressed(P/f'results/validation/{h}_splits.npz',**{'dev_'+k:v for k,v in dev.items()},**{'test_'+k:v for k,v in tests.items()})
  visible=y.copy();visible[du|tu]=np.nan
  train,_,_=detect(visible,t,w,**detector)
  info=dict(detector=detector,models={});localrows=[]
  for family,candidates in CONFIG['models'].items():
   scores_cfg=[];quantiles=[]
   for ci,params in enumerate(candidates):
    start=time.perf_counter();model=fit(family,visible,t,w,train,params);pred=model.grid;seconds=time.perf_counter()-start;part=[];qs={}
    for scenario,q in dev.items():
     r=metrics(y,pred,q,w,scenario!='random');row=dict(heat_id=h,model=family,candidate=ci,scenario=scenario,seconds=seconds,params=json.dumps(params),**r);localrows.append(row);part.append(r['nrmse']);qs[scenario]=float(np.quantile(abs(y[q]-pred[q]),.95))
    scores_cfg.append(float(np.mean(part)));quantiles.append(qs);logging.info('%s model %s cfg %d dev NRMSE %.5f seconds %.2f',h,family,ci,scores_cfg[-1],seconds)
   best=int(np.argmin(scores_cfg));info['models'][family]=dict(candidate=best,params=candidates[best],development_score=scores_cfg[best],interval_q95=quantiles[best])
  tuned[h]=info;devrows.extend(localrows)
  cache.write_text(json.dumps(dict(tuned=info,development=localrows,injection=inj,anomaly=ar),ensure_ascii=False,indent=2))
 pd.DataFrame(devrows).to_csv(P/'results/validation/development_scores.csv',index=False);pd.DataFrame(injection).to_csv(P/'results/validation/anomaly_injection.csv',index=False);pd.DataFrame(anomalyrows).to_csv(P/'results/anomaly_detection/summary.csv',index=False)
 # Decide model family using development results only, before opening any outer test score.
 family_scores={f:float(np.mean([tuned[h]['models'][f]['development_score'] for h in HEATS])) for f in ['A','B','C']};chosen=min(family_scores,key=family_scores.get)
 selection=dict(family=chosen,development_macro_nrmse=family_scores,per_heat=tuned,selection_basis='only development masks; outer tests not used for selection')
 (P/'results/final/model_selection.json').write_text(json.dumps(selection,ensure_ascii=False,indent=2));logging.info('Locked family selection %s %s',chosen,family_scores)
 results=[]
 for h in HEATS:
  cache=P/f'results/validation/{h}_outer.json'
  if cache.exists():results.extend(json.loads(cache.read_text()));logging.info('Resume test %s',h);continue
  d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];splits=np.load(P/f'results/validation/{h}_splits.npz');info=tuned[h];local=[]
  for scenario in [k[5:] for k in splits.files if k.startswith('test_')]:
   q=splits['test_'+scenario];hidden=q.copy()
   if scenario=='time_block':
    i=int(.68*len(t));hidden[i:i+CONFIG['test_time_block_rows']]=True
   if scenario=='wavelength_block':
    j=int(.68*len(w));hidden[:,j:j+CONFIG['test_wavelength_block_channels']]=True
   visible=y.copy();visible[hidden]=np.nan;mask,_,_=detect(visible,t,w,**info['detector'])
   for family in ['A','B','C']:
    start=time.perf_counter();model=fit(family,visible,t,w,mask,info['models'][family]['params']);pred=model.grid;seconds=time.perf_counter()-start
    result=metrics(y,pred,q,w,not scenario.startswith('random'))
    cal='random' if scenario.startswith('random') else scenario;radius=info['models'][family]['interval_q95'][cal]
    f1=float(np.mean([x['f1'] for x in injection if x['heat_id']==h and x['split']=='test']))
    row=dict(heat_id=h,model=family,scenario=scenario,seconds=seconds,process_peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,anomaly_f1=f1,empirical_band_radius=radius,empirical_band_coverage=float(np.mean(abs(pred[q]-y[q])<=radius)),solver_converged=model.diagnostics['converged'],**result)
    (P/f'results/validation/{h}_{family}_{scenario}_solver.json').write_text(json.dumps(model.diagnostics,indent=2));local.append(row);logging.info('Test %s %s %s RMSE %.2f',h,scenario,family,result['rmse'])
  results.extend(local);cache.write_text(json.dumps(local,ensure_ascii=False,indent=2))
 pd.DataFrame(results).to_csv(P/'results/model_comparison/model_comparison.csv',index=False)
 pd.DataFrame(results).groupby(['model','scenario']).mean(numeric_only=True).to_csv(P/'results/model_comparison/aggregate_by_scenario.csv')
 (P/'results/validation/COMPLETE').write_text('Ten heats; all three models; five outer mask scenarios.\n')
if __name__=='__main__':main()
