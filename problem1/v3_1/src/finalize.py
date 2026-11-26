"""Export new native/dense spectra using unchanged V3 completion and interpolation."""
from common31 import *
from stable_completion import complete_with_stability

def one(h):
 frozen=json.loads((P/'config/FINAL_FROZEN.json').read_text());name=frozen['strategy'];y,t,w=load(h);e=dict(np.load(P/f'data/{h}_point_evidence.npz'));remove=np.zeros(y.shape,bool) if name=='V3-P0' else e['high'] if name=='conservative' else e['candidate'];valid=np.isfinite(y)&(y>=0)&(y!=65535)&~remove
 old=dict(np.load(V3/f'data/{h}_interface.npz'));q=old['quality'].copy();grade=old['quality_grade'].copy();sat=y==65535;fallback=np.zeros(y.shape,bool)
 if remove.any():
  r=complete_with_stability(y,t,w,valid);D=r['grid'];fallback=r['numerical_fallback'];latent=r['latent_result']['latent'];s=r['latent_result'];solver={k:v for k,v in s.items() if k not in ['latent','checkpoints']};dump(P/f'results/{h}_solver.json',solver);q[remove]|=4|16384;grade[remove]=2;q[fallback]|=32768
 else:
  D=y.copy();latent=old['latent'];solver=dict(converged=True,stop_reason='no_new_missing_cells_V3_record_grid_reused',history=[],seconds=0.)
 conditional=D.copy();conditional[sat]=old['conditional'][sat];s=Spectrum(t,w,D,y,q,conditional);diff=D-y
 out={**old,'time':t,'wavelength':w,'raw':y,'cleaned':np.where(remove,np.nan,y),'reconstructed':D,'latent':latent,'conditional':conditional,'valid':valid,'invalid':remove,'quality':q,'quality_grade':grade,'statistical_repair':remove,'numerical_fallback':fallback,'point_candidate':e['candidate'],'high_statistical_after_review':e['high'],'suspect_points':e['suspect'],'difference_from_raw':diff,'preprocessing_strategy':name,'hardware_fault_confirmed':False}
 np.savez_compressed(P/f'data/{h}_interface.npz',**out);np.savez_compressed(P/f'data/{h}_difference.npz',difference=diff,row=np.where(remove)[0],column=np.where(remove)[1],raw=y[remove],final=D[remove])
 td=np.linspace(t[0],t[-1],2*len(t)-1);wd=np.linspace(w[0],w[-1],2*len(w)-1);dense=s.evaluate(td,wd).astype(np.float32);np.savez_compressed(P/f'data/{h}_dense.npz',time=td,wavelength=wd,recorded=dense,offline=True)
 frame=pd.read_csv(V3/f'results/{h}_frame_diagnostics.csv');area=np.trapezoid(y,w,axis=1);af=np.trapezoid(D,w,axis=1);peaks=y.max(1);summary=dict(heat_id=h,observations=y.size,candidate_points=int(e['candidate'].sum()),high_statistical_after_review=int(e['high'].sum()),repaired_points=int(remove.sum()),actually_changed_points=int((abs(diff)>1e-8).sum()),repair_fraction=float(remove.mean()),suspect_frames_retained=int(frame.candidate.sum()),whole_frames_deleted=int(remove.all(1).sum()),unrepaired_max_change=float(np.max(abs(diff[~remove]))),repaired_raw_mean=float(y[remove].mean()) if remove.any() else np.nan,repaired_final_mean=float(D[remove].mean()) if remove.any() else np.nan,negative_nodes=int((D<0).sum()),negative_dense=int((dense<0).sum()),sat_record_max_error=float(np.max(abs(D[sat]-65535))) if sat.any() else 0.,conditional_lower_violations=int((conditional[sat]<65535).sum()),peak_p95_relative_change=float(np.quantile(abs(D.max(1)-peaks)/np.maximum(peaks,1),.95)),peak_max_relative_change=float(np.max(abs(D.max(1)-peaks)/np.maximum(peaks,1))),integral_relative_mae=float(np.mean(abs(area-af)/np.maximum(abs(area),1))),numerical_fallback_points=int((fallback&remove).sum()),lowrank_converged=solver['converged'],fit_iterations=len(solver['history']),stop_reason=solver['stop_reason'],fit_seconds=solver['seconds'],interface_sha256=sha(P/f'data/{h}_interface.npz'),dense_sha256=sha(P/f'data/{h}_dense.npz'),source_V3_interface_sha256=sha(V3/f'data/{h}_interface.npz'))
 dump(P/f'results/{h}_final.json',summary);print('final',h,summary['repaired_points'],summary['lowrank_converged'],summary['numerical_fallback_points'],flush=True)
def aggregate():
 summary=[json.loads((P/f'results/{h}_final.json').read_text()) for h in HEATS];pd.DataFrame(summary).to_csv(P/'results/final_summary.csv',index=False);rows=[]
 for h in HEATS:
  d=np.load(P/f'data/{h}_interface.npz');ix,jx=np.where(d['statistical_repair'])
  for i,j in zip(ix,jx):rows.append(dict(heat_id=h,row=int(i),column=int(j),time_s=d['time'][i],wavelength_nm=d['wavelength'][j],raw=d['raw'][i,j],repaired=d['reconstructed'][i,j],difference=d['difference_from_raw'][i,j],numerical_fallback=bool(d['numerical_fallback'][i,j]),statistical_anomaly_not_confirmed_hardware_fault=True))
 pd.DataFrame(rows,columns=['heat_id','row','column','time_s','wavelength_nm','raw','repaired','difference','numerical_fallback','statistical_anomaly_not_confirmed_hardware_fault']).to_csv(P/'results/actual_point_repairs.csv',index=False)
 dump(P/'results/FINAL_COMPLETE.json',dict(heats=10,strategy=json.loads((P/'config/FINAL_FROZEN.json').read_text())['strategy'],repaired_points=sum(r['repaired_points'] for r in summary),whole_frames_deleted=0,model_frozen=True))
if __name__=='__main__':aggregate() if sys.argv[1]=='aggregate' else one(sys.argv[1])
