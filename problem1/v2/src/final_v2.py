"""Final full-native refits under the locked development choice; V1 output untouched."""
import json,time
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *
from solver import fit
from detector import diagnostics
from interpolation import Surface
from models import fit_c

def main():
 if not (P/'results/TEST_COMPLETE.json').exists():raise RuntimeError('Independent completion test must finish first')
 lock=json.loads((P/'config/LOCKED.json').read_text());final=json.loads((P/'config/FINAL_SELECTION_LOCKED.json').read_text());interp=json.loads((P/'config/INTERPOLATION_LOCKED.json').read_text());variant=final['model'];out=[];curves=[]
 for h in HEATS:
  y,t,w=load(h);dc=lock['detector'];mask,types,diag=diagnostics(y,t,w,**dc)
  if lock['detector_method']=='MAD':mask,types,_=detect(y,t,w,window=5,threshold=dc['threshold'])
  sat=y==65535;start=time.perf_counter()
  if variant=='C':
   cp=json.loads((P/'config/C_EXPANDED_LOCKED.json').read_text())['params'].copy();ratio_t=cp['kt']/96;ratio_w=cp['kw']/192;cp['kt']=min(len(t),max(4,int(np.ceil(len(t)*ratio_t))));cp['kw']=min(len(w),max(4,int(np.ceil(len(w)*ratio_w))));cp|=dict(irls=30,cg_iter=250,cg_tol=1e-6,irls_tol=1e-4,objective_tol=1e-6)
   model=fit_c(y,t,w,mask,cp);grid=model.grid;latent=grid;hist=model.diagnostics['history'];converged=model.diagnostics['converged'];model.save(P/f'data/{h}_native_tensor_model.npz');params=cp
   # Exact TensorSurface remains evaluate interface; also provide gridded interpolation export for consistency.
  else:
   params=lock['models'][variant]['params'];r=fit(y,t,w,mask,variant,**params,max_iter=CFG['solver_max_iter'],tol=CFG['solver_tol']);grid=r['grid'];latent=r['latent'];hist=r['history'];converged=r['converged']
  model=Surface(t,w,grid,interp['method']);model.save(P/f'data/{h}_model.npz',latent_grid=latent,model=variant,params_json=json.dumps(params))
  np.savez_compressed(P/f'data/{h}_masks.npz',valid=mask,types=types,saturated=sat,saturated_unidentifiable=sat,time=t,wavelength=w)
  np.savez_compressed(P/f'data/{h}_diagnostics.npz',frame_angle_deg=diag['frame_angle_deg'],frame_brightness_relative=diag['frame_brightness_relative'],frame_residual=diag['frame_residual'],channel_residual=diag['channel_residual'])
  pd.DataFrame(dict(time_s=t,brightness_relative=diag['frame_brightness_relative'],angle_deg=diag['frame_angle_deg'],lowrank_residual=diag['frame_residual'])).to_csv(P/f'results/{h}_frame_diagnostics.csv',index=False)
  pd.DataFrame(dict(wavelength_nm=w,lowrank_persistent_residual=diag['channel_residual'])).to_csv(P/f'results/{h}_channel_diagnostics.csv',index=False)
  dump(P/f'results/{h}_final_solver.json',dict(model=variant,converged=converged,history=hist,params=params));curves.extend(dict(heat_id=h,**x) for x in hist)
  tt=np.linspace(t[0],t[-1],140);ww=np.linspace(w[0],w[-1],300);T,W=np.meshgrid(tt,ww,indexing='ij');pred=model.evaluate(tt,ww)
  fig=plt.figure(figsize=(10,6));ax=fig.add_subplot(111,projection='3d');ax.plot_surface(T,W,pred,cmap='viridis',rasterized=True,linewidth=0);ax.set(title=h+' V2 '+variant+' + '+interp['method'],xlabel='Time (s)',ylabel='Wavelength (nm)',zlabel='Intensity');ax.view_init(28,-130);fig.tight_layout();fig.savefig(P/f'figures/{h}_surface.png',dpi=170);plt.close(fig)
  fig,axs=plt.subplots(2,1,figsize=(10,6));axs[0].imshow((~mask).T,aspect='auto',origin='lower',extent=[t[0],t[-1],w[0],w[-1]],interpolation='nearest');axs[0].set(title=h+' anomaly/missing mask',ylabel='Wavelength (nm)');axs[1].imshow(sat.T,aspect='auto',origin='lower',extent=[t[0],t[-1],w[0],w[-1]],interpolation='nearest');axs[1].set(title='Conditional saturation: unidentifiable upper intensity',xlabel='Time (s)',ylabel='Wavelength (nm)');fig.tight_layout();fig.savefig(P/f'figures/{h}_masks.png',dpi=150);plt.close(fig)
  out.append(dict(heat_id=h,model=variant,interpolation=interp['method'],seconds=time.perf_counter()-start,converged=converged,valid_observation_rmse_not_prediction=float(np.sqrt(np.mean((grid[mask]-y[mask])**2))),masked_fraction=float((~mask).mean()),saturated_count=int(sat.sum()),saturated_lower_violation=float(np.mean(grid[sat]<65535)) if sat.any() else 0,grid_negative_fraction=float(np.mean(grid<0)),function_query_negative_fraction=float(np.mean(pred<0)),params=json.dumps(params),upper_bound_identifiable=False))
  pd.DataFrame(out).to_csv(P/'results/final_summary.csv',index=False);print('final',h,variant,converged,flush=True)
 pd.DataFrame(curves).to_csv(P/'results/final_solver_curves.csv',index=False)
 dump(P/'results/FINAL_COMPLETE.json',dict(heats=HEATS,model=variant,interpolation=interp['method'],selection_sha256=digest(P/'config/FINAL_SELECTION_LOCKED.json')))
if __name__=='__main__':main()
