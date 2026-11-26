"""Post-lock compute cost and solver curve exports; no new accuracy selection."""
from ctx import *
from preprocess import apply
from model3 import project_record
import pandas as pd,time,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Noto Serif CJK SC','axes.unicode_minus':False})
def main():
 rows=[]
 for h,stage,y,t,w,_,_ in panels('model'):
  if stage!=1:continue
  hidden=dict(hides(y.shape,9307))['random_0.05'];v=y.copy();v[hidden]=np.nan;pp=apply(v,t,w,'P0');mask=pp['valid'];par=json.loads((V1/'results/final/model_selection.json').read_text())['per_heat'][h]['models']['B']['params']
  for name in ['V1','V2','V3-A','V3-B','V3-C']:
   start=time.perf_counter()
   if name=='V1':m=fit_b(v,t,w,mask,par);grid=m.grid
   else:
    r=legacy.fit(v,t,w,mask,'B1' if name=='V3-B' else 'B2',rank=4,eta=0 if name=='V3-B' else .003,max_iter=100,tol=1e-4);grid=r['grid']
    if name.startswith('V3'):grid=project_record(v,mask,r['latent'],t,w,name=='V3-C',name=='V3-C')[0]
   Surface(t,w,grid,'pchip' if name=='V1' else 'linear').evaluate(t,w)
   rows.append(dict(heat_id=h,model=name,seconds=time.perf_counter()-start,shape='32x256',common_preprocessing_excluded=True,repetitions=1))
 pd.DataFrame(rows).to_csv(P/'results/paired_compute_cost.csv',index=False)
 for h in CFG['stability_heats']:
  f,axs=plt.subplots(2,2,figsize=(11,7))
  for scope in ['native','held']:
   hist=json.loads((P/f'results/{h}_{scope}_solver.json').read_text())['history'];x=[r['iteration'] for r in hist];axs[0,0].plot(x,[r['objective'] for r in hist],label=scope);axs[0,1].semilogy(x,[r['relative_change'] for r in hist],label=scope);axs[1,1].plot(x,[r['seconds'] for r in hist],label=scope)
  q=pd.read_csv(P/f'results/{h}_stability.csv');z=q[q.scope=='held_panel'];axs[1,0].plot(z.budget,z.held_rmse,'o-',label='遮蔽点RMSE');axs[0,0].set(title='实际目标函数',xlabel='迭代');axs[0,1].set(title='相对变化',xlabel='迭代');axs[0,1].axhline(1e-4,color='gray',ls='--');axs[1,0].set(title='预算与恢复误差',xlabel='迭代预算',ylabel='counts');axs[1,1].set(title='累计计算时间',xlabel='迭代',ylabel='s')
  for a in axs.flat:a.legend()
  f.suptitle(h+'：原始求解状态，不因输出保真掩盖潜变量变化');f.tight_layout();f.savefig(P/f'figures/appendix/{h}_求解器完整曲线.png',dpi=180);plt.close(f)
if __name__=='__main__':main()
