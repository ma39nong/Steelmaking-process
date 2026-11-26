"""Final refits, continuous evaluator exports, ten surfaces, BEADS and censor sensitivity."""
import json,time,logging,sys
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from pybaselines import Baseline
from audit import ROOT,P,HEATS
from models import fit,load_surface
from anomalies import fill_training,detect
from eda import save,heatmap

def surface_plot(model,path,title):
 t=np.linspace(model.time[0],model.time[-1],160);w=np.linspace(model.wave[0],model.wave[-1],360);z=model.evaluate(t,w);T,W=np.meshgrid(t,w,indexing='ij')
 fig=plt.figure(figsize=(10,6));ax=fig.add_subplot(111,projection='3d');ax.plot_surface(T,W,z,cmap='viridis',rstride=1,cstride=1,linewidth=0,antialiased=True,rasterized=True);ax.set(xlabel='Time (s)',ylabel='Wavelength (nm)',zlabel='Intensity (recorded units)',title=title);ax.view_init(28,-130);save(fig,path)

def beads_experiments():
 out=[]
 for h in HEATS:
  d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];mask=np.load(P/f'data/masks/{h}.npz')['valid'];clean=fill_training(np.where(mask,y,np.nan),t,w)
  for stage,ix in [('early',len(t)//10),('middle',len(t)//2),('late',len(t)*9//10)]:
   uniform=np.linspace(w[0],w[-1],len(w));curve=np.interp(uniform,w,clean[ix]);sc=max(curve.max(),1);z=curve/sc
   fig,ax=plt.subplots(figsize=(10,4));ax.plot(uniform,curve,color='black',lw=.8,label='Observed (training-cleaned)')
   for fc,lam in [(.08,[.1,1,.8]),(.008,[.001,.01,.008])]:
    baseline,info=Baseline().beads(z,freq_cutoff=fc,lam_0=lam[0],lam_1=lam[1],lam_2=lam[2],asymmetry=6,filter_type=1,max_iter=100,tol=1e-4)
    peak=info['signal'];recovered=peak+baseline
    out.append(dict(heat_id=h,stage=stage,time_s=t[ix],freq_cutoff=fc,lam0=lam[0],lam1=lam[1],lam2=lam[2],iterations=len(info['tol_history']),converged=bool(info['tol_history'][-1]<1e-4),last_tolerance=float(info['tol_history'][-1]),baseline_area_fraction=float(np.trapezoid(baseline,uniform)/np.trapezoid(z,uniform)),peak_only_area_relative_error=float(abs(np.trapezoid(peak-z,uniform))/max(abs(np.trapezoid(z,uniform)),1e-12)),combined_rmse=float(np.sqrt(np.mean((recovered-z)**2))*sc),peak_only_max_relative_error=float(abs(peak.max()-z.max())/max(z.max(),1e-12))))
    np.savez_compressed(P/f'results/validation/{h}_BEADS_{stage}_{fc}.npz',wavelength=uniform,observed=curve,baseline=baseline*sc,peaks=peak*sc,recovered=recovered*sc,residual=(z-recovered)*sc)
    ax.plot(uniform,baseline*sc,lw=.7,label=f'Baseline fc={fc}');ax.plot(uniform,peak*sc,lw=.7,label=f'Peaks only fc={fc}')
   ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (recorded units)',title=f'{h} BEADS {stage}; uniform resampling explicit');ax.legend(fontsize=8);save(fig,P/f'figures/comparison/{h}_BEADS_{stage}')
 pd.DataFrame(out).to_csv(P/'results/validation/beads_experiment.csv',index=False)

def main():
 if not (P/'results/validation/COMPLETE').exists():raise RuntimeError('Independent tests must finish first')
 from provenance import check_or_stamp
 check_or_stamp()
 logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s')
 selection=json.loads((P/'results/final/model_selection.json').read_text());family=selection['family'];out=[];satrows=[]
 for h in HEATS:
  d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];m=np.load(P/f'data/masks/{h}.npz');mask=m['valid'];types=m['types'];params=selection['per_heat'][h]['models'][family]['params'].copy()
  # Final refit uses precisely the validation-selected solver/parameters. No test-driven changes.
  start=time.perf_counter();model=fit(family,y,t,w,mask,params);pred=model.grid;seconds=time.perf_counter()-start
  model.save(P/f'data/reconstructed/{h}_model.npz')
  td=np.linspace(t[0],t[-1],2*len(t)-1);wd=np.linspace(w[0],w[-1],2*len(w)-1);dense=model.evaluate(td,wd).astype(np.float32)
  np.savez_compressed(P/f'data/reconstructed/{h}_dense.npz',time=td,wavelength=wd,intensity=dense)
  np.savez_compressed(P/f'data/reconstructed/{h}_on_observed_grid.npz',time=t,wavelength=w,intensity=pred)
  radius=selection['per_heat'][h]['models'][family]['interval_q95'];sat=(types&2)>0
  np.savez_compressed(P/f'data/reconstructed/{h}_uncertainty.npz',saturated_unidentifiable=sat,valid_reference=mask,random_empirical_radius=radius['random'],time_block_empirical_radius=radius['time_block'],wavelength_block_empirical_radius=radius['wavelength_block'])
  surface_plot(model,P/f'figures/reconstructed/{h}_surface',h+' continuous reconstruction (model '+family+')')
  heatmap(pred,t,w,P/f'figures/reconstructed/{h}_heatmap',h+' reconstruction (model '+family+')')
  fig,axs=plt.subplots(1,2,figsize=(12,4))
  for i in [len(t)//10,len(t)//2,len(t)*9//10]:
   axs[0].plot(w,y[i],lw=.6,alpha=.6,label=f'Raw {t[i]:g}s');axs[0].plot(w,pred[i],lw=.8,label=f'Fit {t[i]:g}s')
  axs[0].set(xlabel='Wavelength (nm)',ylabel='Intensity');axs[0].legend(fontsize=7)
  ii,jj=np.where(~mask)
  if len(ii):
   k=len(ii)//2;i,j=ii[k],jj[k];ix=np.arange(max(0,j-25),min(len(w),j+26));axs[1].plot(w[ix],y[i,ix],'.-',lw=.7,label='Raw');axs[1].plot(w[ix],pred[i,ix],lw=1,label='Reconstruction');bad=~mask[i,ix];axs[1].scatter(w[ix][bad],y[i,ix][bad],c='red',s=12,label='Masked');axs[1].set_title(f't={t[i]:g}s')
  axs[1].set(xlabel='Wavelength (nm)',ylabel='Intensity');axs[1].legend(fontsize=8);fig.suptitle(h);save(fig,P/f'figures/comparison/{h}_before_after')
  # Independent censor sensitivity is not silently substituted into the selected final family.
  if sat.any():
   cp=selection['per_heat'][h]['models']['C']['params'].copy();masked_c=fit('C',y,t,w,mask,cp).grid;cp['censor_weight']=.5;censored=fit('C',y,t,w,mask,cp).grid
   ii,jj=np.where(sat);pd.DataFrame(dict(time_s=t[ii],wavelength_nm=w[jj],recorded=y[sat],selected_mask_reconstruction=pred[sat],spline_mask=masked_c[sat],spline_one_sided=censored[sat],true_upper_bound_known=False)).to_csv(P/f'results/final/{h}_saturation_sensitivity.csv',index=False)
   satrows.append(dict(heat_id=h,count=int(sat.sum()),selected_below_upper_code_fraction=float(np.mean(pred[sat]<65535)),mask_c_shortfall_mean=float(np.mean(np.maximum(65535-masked_c[sat],0))),censor_c_shortfall_mean=float(np.mean(np.maximum(65535-censored[sat],0))),model_disagreement_mean=float(np.mean(abs(masked_c[sat]-censored[sat]))),true_error_identifiable=False))
  out.append(dict(heat_id=h,model=family,seconds=seconds,params=json.dumps(params),solver_converged=model.diagnostics['converged'],negative_fraction=float(np.mean(pred<0)),dense_negative_fraction=float(np.mean(dense<0)),dense_minimum=float(dense.min()),dense_maximum=float(dense.max()),observed_grid_rmse_not_validation=float(np.sqrt(np.mean((pred[mask]-y[mask])**2))),saturated_count=int(sat.sum()),dense_time=len(td),dense_wavelength=len(wd),model_file=f'data/reconstructed/{h}_model.npz'))
  logging.info('Final %s model %s seconds %.2f',h,family,seconds)
 pd.DataFrame(out).to_csv(P/'results/final/reconstruction_summary.csv',index=False);pd.DataFrame(satrows).to_csv(P/'results/final/saturation_sensitivity.csv',index=False)
 fig=plt.figure(figsize=(18,15))
 for k,h in enumerate(HEATS):
  model=load_surface(P/f'data/reconstructed/{h}_model.npz');tt=np.linspace(model.time[0],model.time[-1],70);ww=np.linspace(model.wave[0],model.wave[-1],120);T,W=np.meshgrid(tt,ww,indexing='ij');ax=fig.add_subplot(5,2,k+1,projection='3d');ax.plot_surface(T,W,model.evaluate(tt,ww),cmap='viridis',linewidth=0,rasterized=True);ax.set(title=h,xlabel='Time (s)',ylabel='Wavelength (nm)',zlabel='Intensity');ax.view_init(25,-130)
 save(fig,P/'figures/reconstructed/all_heats_surfaces')
 if '--skip-beads' not in sys.argv:beads_experiments()
 elif not (P/'results/validation/beads_experiment.csv').exists():raise RuntimeError('Cannot skip missing BEADS results')
 (P/'results/final/COMPLETE').write_text('Ten real-data continuous models, dense grids and individual surfaces generated.\n')
if __name__=='__main__':main()
