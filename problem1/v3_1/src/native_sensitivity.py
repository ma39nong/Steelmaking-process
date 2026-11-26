"""Native four-point conservative counterfactual only; never replaces frozen formal P0."""
from common31 import *
from stable_completion import complete_with_stability

def main(h):
 assert json.loads((P/'config/FINAL_FROZEN.json').read_text())['strategy']=='V3-P0'
 y,t,w=load(h);e=np.load(P/f'data/{h}_point_evidence.npz');mask=e['high'];assert mask.any();valid=(y!=65535)&~mask;r=complete_with_stability(y,t,w,valid);D=r['grid'];s=r['latent_result'];dump(P/f'results/{h}_sensitivity_solver.json',{k:v for k,v in s.items() if k not in ['latent','checkpoints']})
 rows=[];ii,jj=np.where(mask)
 for i,j in zip(ii,jj):
  neighbor=np.r_[y[i-2:i,j],y[i+1:i+3,j],y[i,j-2:j],y[i,j+1:j+3]]
  rows.append(dict(heat_id=h,row=int(i),column=int(j),time_s=t[i],wavelength_nm=w[j],raw=y[i,j],counterfactual_repaired=D[i,j],difference=D[i,j]-y[i,j],neighbor_min=float(neighbor.min()),neighbor_max=float(neighbor.max()),within_neighbor_support=bool(neighbor.min()<=D[i,j]<=neighbor.max()),numerical_fallback=bool(r['numerical_fallback'][i,j]),lowrank_converged=s['converged'],iterations=len(s['history']),seconds=s['seconds'],formal_repair_applied=False))
 pd.DataFrame(rows).to_csv(P/f'results/{h}_counterfactual_repair.csv',index=False);np.savez_compressed(P/f'data/{h}_counterfactual_patch.npz',raw=y[max(0,i-5):i+6,max(0,j-10):j+11],repaired=D[max(0,i-5):i+6,max(0,j-10):j+11],time=t[max(0,i-5):i+6],wavelength=w[max(0,j-10):j+11]);print('sensitivity',h,rows,flush=True)
if __name__=='__main__':main(sys.argv[1])
