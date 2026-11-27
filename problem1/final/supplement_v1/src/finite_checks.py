"""Finite supplementary channel recovery and diagnostic replay; frozen release untouched."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.ndimage import maximum_filter1d
S=Path(__file__).resolve().parents[1];F=S.parent
sys.path.insert(0,str(F/'src'))
import core as C
from detectors import adaptive,local_mask
from analyze import channel_reference


def peaks(ref,corrupt,D):
 env=maximum_filter1d(ref,size=5,axis=1,mode='nearest');threshold=np.maximum(500,.1*ref.max(1)[:,None]);loc=(D>=np.roll(D,1,1))&(D>=np.roll(D,-1,1));loc[:,:2]=False;loc[:,-2:]=False
 oldloc=(corrupt>=np.roll(corrupt,1,1))&(corrupt>=np.roll(corrupt,-1,1));return loc&(D>env+threshold)&~(oldloc&(corrupt>env+threshold))


def main():
 rows=[];panels=[];guard=[]
 for hi,h in enumerate(C.HEATS):
  y,t,w=C.load(h);starts=np.linspace(0,len(t)-96,6).astype(int)
  for split,gap in [('development',3),('confirmation',4)]:
   i=int((starts[gap]+96+starts[gap+1]-32)//2)
   for col in [800,1600]:
    ref=y[i:i+32,col:col+96].copy();tt=t[i:i+32];ww=w[col:col+96];_,d,scale=channel_reference(ref,ww);zz=np.median(abs(d)/scale,axis=0);choices=np.arange(32,64);j=int(choices[np.argmin(zz[choices])]);panels.append(dict(heat_id=h,split=split,row_start=i,column_start=col,height=32,width=96,injected_column=col+j))
    for kind in ['offset','gain','freeze']:
     for strength in [.3,.7]:
      v=ref.copy();rng=np.random.default_rng(27017+hi+col);amplitude=strength*max(5000,float(np.quantile(ref,.95)-np.quantile(ref,.05)))
      if kind=='offset':v[:,j]=np.minimum(64000,ref[:,j]+amplitude)
      elif kind=='gain':v[:,j]=np.minimum(64000,ref[:,j]*(1+strength))
      else:v[:,j]=(1-strength)*ref[:,j]+strength*np.median(ref[:,j])
      truth=v!=ref;rr,dd,ss=channel_reference(v,ww);md=np.median(dd,axis=0);pers=np.maximum((dd>0).mean(0),(dd<0).mean(0));zz=np.median(abs(dd)/ss,axis=0);pred=(pers>=.95)&(abs(md)>=250)&(zz>=4);pred[:4]=False;pred[-4:]=False
      tp=int(pred[j]);fp=int(pred.sum()-tp);fn=1-tp
      oracle=np.zeros(v.shape,bool);oracle[:,j]=True;valid=(v!=65535)&~oracle
      baseline=np.sqrt(np.mean((v[truth]-ref[truth])**2))
      for method in ['bidirectional_linear','rank4_huber']:
       if method=='bidirectional_linear':latent=C.fill_training(np.where(valid|(v==65535),v,np.nan),tt,ww);iterations=0;conv=True
       else:r=C.solve(v,tt,ww,valid,max_iter=100);latent=r['latent'];iterations=len(r['history']);conv=r['converged']
       D=C.project_record(v,valid,latent,tt,ww,True,True)[0];ar=np.trapezoid(ref,ww,axis=1);ad=np.trapezoid(D,ww,axis=1)
       rows.append(dict(heat_id=h,split=split,column_start=col,kind=kind,strength=strength,method=method,tp=tp,fp=fp,fn=fn,before_rmse=baseline,repair_rmse=float(np.sqrt(np.mean((D[truth]-ref[truth])**2))),unpolluted_rmse=float(np.sqrt(np.mean((D[~truth]-ref[~truth])**2))),peak_relative=float(np.mean(abs(D.max(1)-ref.max(1))/np.maximum(ref.max(1),1))),integral_relative=float(np.mean(abs(ad-ar)/np.maximum(abs(ar),1))),introduced_peaks=int(peaks(ref,v,D).sum()),oracle_mask=True,iterations=iterations,converged=conv))
    # One predeclared guard check per type at col800, using the unchanged released pipeline.
    if col==800:
     for kind in ['spike','drop']:
      v,truth=C.inject(ref,np.random.default_rng(27019+hi),kind,.5);e=adaptive(v,tt,ww);mask=local_mask(e,4);valid=(v!=65535)&~mask
      if mask.any():r=C.complete_with_stability(v,tt,ww,valid);D=r['grid'];fall=r['numerical_fallback']
      else:D=v.copy();fall=np.zeros(v.shape,bool)
      protect=mask&fall&(v>=.9*65535)&(e['wave_reference']>=.9*65535) if 'wave_reference' in e else mask&fall&(v>=.9*65535)&(e['wavelength_reference']>=.9*65535)
      wr=e['wave_reference'] if 'wave_reference' in e else e['wavelength_reference'];tr=e['time_reference'];protect&=abs(tr-wr)>.5*np.maximum(abs(wr),25)
      alt=D.copy();alt[protect]=wr[protect];alt=C.project_record(v,valid,alt,tt,ww,True,True)[0]
      for method,out in [('frozen',D),('guard_candidate',alt)]:guard.append(dict(heat_id=h,split=split,kind=kind,method=method,protected_count=int(protect.sum()),pollution_rmse=float(np.sqrt(np.mean((out[truth]-ref[truth])**2))),overall_rmse=float(np.sqrt(np.mean((out-ref)**2))),unpolluted_rmse=float(np.sqrt(np.mean((out[~truth]-ref[~truth])**2))),introduced_peaks=int(peaks(ref,v,out).sum())))
  pd.DataFrame(rows).to_csv(S/'results/channel_validation.csv',index=False);pd.DataFrame(guard).to_csv(S/'results/guard_validation.csv',index=False);print('finite checks',h,flush=True)
 pd.DataFrame(panels).to_csv(S/'results/new_panel_index.csv',index=False)
 heat=pd.DataFrame(rows).groupby(['split','heat_id','kind','method']).mean(numeric_only=True).reset_index();heat.groupby(['split','kind','method']).mean(numeric_only=True).reset_index().to_csv(S/'results/channel_validation_macro.csv',index=False)
 # Replay only historically identified introduced-peak cases; inspect locations, not select parameters.
 old=pd.read_csv(F/'results/exact_final_pipeline_confirmation.csv');old=old[old.introduced_unsupported_peak_count>0];lookup={(h,b):(ref,t,w,i,j) for h,b,ref,t,w,i,j in C.panels('confirmation')};locations=[];cases=[]
 for row in old.itertuples():
  ref,t,w,i,j=lookup[(row.heat_id,row.band)];v,truth=C.inject(ref,np.random.default_rng(row.seed+sum(map(ord,row.heat_id))+j),row.kind,row.strength);e=adaptive(v,t,w);mask=local_mask(e,4);r=C.complete_with_stability(v,t,w,(v!=65535)&~mask);D=r['grid'];new=peaks(ref,v,D);count=int(new.sum());cases.append(dict(heat_id=row.heat_id,band=row.band,seed=row.seed,kind=row.kind,strength=row.strength,expected=row.introduced_unsupported_peak_count,replayed=count))
  assert count==row.introduced_unsupported_peak_count,(row,count)
  for a,b in zip(*np.where(new)):
   locations.append(dict(heat_id=row.heat_id,band=row.band,seed=row.seed,kind=row.kind,strength=row.strength,row=int(a+i),column=int(b+j),time_s=t[a],wavelength_nm=w[b],reference=ref[a,b],corrupted=v[a,b],repaired=D[a,b],was_injected=bool(truth[a,b]),was_masked=bool(mask[a,b]),fallback=bool(r['numerical_fallback'][a,b]),score=float(e['score'][a,b])))
  print('peak replay',row.heat_id,row.band,count,flush=True)
 assert len(locations)==17
 pd.DataFrame(locations).to_csv(S/'results/introduced_peaks_17_locations.csv',index=False);pd.DataFrame(cases).to_csv(S/'results/introduced_peak_replay_cases.csv',index=False)
 print('ALL_FINITE_CHECKS_COMPLETE',flush=True)
if __name__=='__main__':main()
