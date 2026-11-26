"""Final fits: unchanged V3 solver, supports and recorded-scale bilinear export."""
from core import *

def one(h):
 y,t,w=load(h);m=dict(np.load(P/f'data/{h}_masks.npz'));mask=m['repair'];valid=m['valid'];r=complete_with_stability(y,t,w,valid);D=r['grid'];fit=r['latent_result'];dump(P/f'results/{h}_solver.json',{k:v for k,v in fit.items() if k not in ['latent','checkpoints']})
 np.savez_compressed(P/f'data/{h}_fit.npz',reconstructed=D,latent=fit['latent'],numerical_fallback=r['numerical_fallback']);print('fit',h,int(mask.sum()),fit['converged'],int(r['numerical_fallback'][mask].sum()),flush=True)
if __name__=='__main__':one(sys.argv[1])
