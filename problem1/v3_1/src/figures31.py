from common31 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Noto Serif CJK SC','font.size':10,'axes.unicode_minus':False})
A=P/'figures/appendix';M=P/'figures/main'
def save(fig,name,main=False):
 folder=M if main else A;fig.savefig(folder/(name+'.png'),dpi=160,bbox_inches='tight');path=folder/(name+'.svg');fig.savefig(path,bbox_inches='tight');path.write_text('\n'.join(x.rstrip() for x in path.read_text().splitlines())+'\n');plt.close(fig)
def main():
 assert (P/'results/FINAL_COMPLETE.json').exists();high=pd.read_csv(P/'results/high_confidence_statistical_points.csv');overview=plt.figure(figsize=(14,22));overview.subplots_adjust(hspace=.5,wspace=.05)
 for z,h in enumerate(HEATS):
  d=np.load(P/f'data/{h}_interface.npz');y=d['raw'];D=d['reconstructed'];t=d['time'];w=d['wavelength'];q=high[high.heat_id==h]
  if len(q):i=int(q.iloc[0].row);j=int(q.iloc[0].column)
  else:
   ix,jx=np.where(d['point_candidate']);i=int(ix[len(ix)//2]);j=int(jx[len(jx)//2])
  a=max(0,i-5);b=min(len(t),i+6);c=max(0,j-10);e=min(len(w),j+11);fig,axs=plt.subplots(2,3,figsize=(15,8))
  axs[0,0].plot(w,y[i],color='black',label='原始记录');axs[0,0].plot(w,D[i],ls='--',label='正式输出');axs[0,0].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title=f't={t[i]} s 全谱');axs[0,0].legend()
  axs[0,1].plot(w[c:e],y[i,c:e],label='原始');axs[0,1].plot(w[c:e],D[i,c:e],ls='--',label='正式P0');axs[0,1].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title='候选邻域（正式保留）');axs[0,1].legend()
  axs[0,2].plot(t[a:b],y[a:b,j],label='原始');axs[0,2].plot(t[a:b],D[a:b,j],ls='--',label='正式');axs[0,2].set(xlabel='时间 / s',ylabel='记录强度 / counts',title=f'λ={w[j]} nm');axs[0,2].legend()
  low=y[a:b,c:e].min();up=y[a:b,c:e].max()
  for ax,img,title in zip(axs[1],[y[a:b,c:e],D[a:b,c:e],D[a:b,c:e]-y[a:b,c:e]],['原始局部热图','正式输出局部热图','正式−原始（全零）']):
   diff='全零' in title;im=ax.imshow(img,aspect='auto',origin='lower',extent=[w[c],w[e-1],t[a],t[b-1]],cmap='coolwarm' if diff else 'viridis',vmin=-1 if diff else low,vmax=1 if diff else up);ax.set(xlabel='波长 / nm',ylabel='时间 / s',title=title);fig.colorbar(im,ax=ax,label='counts')
  fig.suptitle(h+'：确认未通过，正式数据未作点修复');fig.tight_layout();save(fig,h+'_处理前后与局部热图')
  fig,ax=plt.subplots(figsize=(9,4));ax.imshow(D-y,origin='lower',aspect='auto',extent=[w[0],w[-1],t[0],t[-1]],cmap='coolwarm',vmin=-1,vmax=1);ax.set(xlabel='波长 / nm',ylabel='时间 / s',title=h+' 正式矩阵−原始矩阵：修改单元数0，最大绝对差0');save(fig,h+'_全矩阵差值')
  ii=np.linspace(0,len(t)-1,60).astype(int);jj=np.linspace(0,len(w)-1,90).astype(int);W,T=np.meshgrid(w[jj],t[ii])
  for fig,ax in [(plt.figure(figsize=(9,6)),None),(overview,None)]:
   ax=fig.add_subplot(1,1,1,projection='3d') if fig is not overview else fig.add_subplot(5,2,z+1,projection='3d');ax.plot_surface(W,T,D[np.ix_(ii,jj)],cmap='viridis',linewidth=0,rasterized=True);ax.set(title=h,xlabel='波长 / nm',ylabel='时间 / s',zlabel='记录强度 / counts');ax.view_init(25,-125)
   if fig is not overview:save(fig,h+'_最终连续光谱')
 save(overview,'十炉最终连续光谱总览',True)
 fig,ax=plt.subplots(figsize=(13,4));ax.axis('off');items=[(.02,'原始记录\n27,816,683个观测'),(.27,'统计证据分级\n5强候选→复核保留1'),(.53,'端到端独立确认\n保守策略未通过'),(.78,'正式P0保留\n诊断＋连续重构')]
 # Count is generated from actual summary, never hand-assumed.
 total=int(pd.read_csv(P/'results/raw_point_audit.csv').observations.sum());items[0]=(.02,f'原始记录\n{total:,}个观测')
 for x,text in items:ax.text(x,.5,text,transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.7',fc='#e6f1ed',ec='#277f8e'),va='center')
 for x in [.22,.48,.73]:ax.annotate('',xy=(x+.035,.5),xytext=(x,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
 ax.text(.02,.1,'确认Recall 0.071，污染RMSE降低3.93%；实际修复0点，不宣称完成异常消除。',transform=ax.transAxes);save(fig,'最终预处理技术路线',True)
 fig,axs=plt.subplots(1,2,figsize=(12,4))
 for ax,split in zip(axs,['development','confirmation']):
  q=pd.read_csv(P/f'results/contamination_{split}.csv');g=q.groupby('strategy')[['pollution_rmse','unpolluted_rmse','overall_rmse']].mean();g.plot.bar(ax=ax,color=['#277f8e','#db863c','#68aa71']);ax.set(title='开发' if split=='development' else '新区域确认',ylabel='RMSE / counts',xlabel='策略');ax.tick_params(axis='x',rotation=0)
 fig.tight_layout();save(fig,'异常污染端到端恢复比较',True)
 # Exact saved protocol case, for illustration only; no score or tuning changes.
 from experiments import inject
 from points import evidence,repair
 for split,col in [('development',64),('confirmation',1720)]:
  for h,y,t,w,_,j in panels(split):
   if h!='6B03205' or j!=col:continue
   seed=CFG[split+'_seeds'][0];rng=np.random.default_rng(seed+sum(map(ord,h))+j);corrupt,truth=inject(y,rng,'spike',.7);ev=evidence(corrupt,t,w);C=repair(corrupt,t,w,ev,'conservative')[0];B=repair(corrupt,t,w,ev,'V3-P1')[0];ix,jx=np.where(truth);i=int(ix[0]);k=int(jx[0]);c=max(0,k-8);e=min(len(w),k+9);fig,axs=plt.subplots(1,2,figsize=(12,4))
   for a,yy,title in zip(axs,[slice(None),slice(c,e)],['全谱','注入点局部']):
    for data,label,style in [(y,'注入前带噪参考','-'),(corrupt,'污染输入 / P0',':'),(C,'保守处理','--'),(B,'历史P1','-.')]:a.plot(w[yy],data[i,yy],style,label=label)
    a.set(xlabel='波长 / nm',ylabel='记录强度 / counts',title=f'{split}，{title}，t={t[i]} s');a.legend(fontsize=8)
   save(fig,'固定注入案例_'+split,True)

def sensitivity():
 q=pd.concat([pd.read_csv(P/f'results/{h}_counterfactual_repair.csv') for h in ['6B03205','6B03270','6B03271','6C02059']]);q.to_csv(P/'results/counterfactual_repairs.csv',index=False)
 for r in q.itertuples():
  d=np.load(P/f'data/{r.heat_id}_counterfactual_patch.npz');y=d['raw'];D=d['repaired'];t=d['time'];w=d['wavelength'];i=int(np.argmin(abs(t-r.time_s)));fig,axs=plt.subplots(1,3,figsize=(14,4));axs[0].plot(w,y[i],label='正式P0 / 原始');axs[0].plot(w,D[i],ls='--',label='保守修复敏感性（非正式）');axs[0].legend(fontsize=8);axs[0].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title='局部谱对照')
  for ax,img,title in zip(axs[1:],[y,D],['正式P0原始局部','敏感性修复局部']):im=ax.imshow(img,origin='lower',aspect='auto',extent=[w[0],w[-1],t[0],t[-1]],vmin=y.min(),vmax=y.max(),cmap='viridis');ax.set(xlabel='波长 / nm',ylabel='时间 / s',title=title);fig.colorbar(im,ax=ax)
  fig.suptitle(r.heat_id+'：仅检查修复后支撑结构，不宣称得到真实值');fig.tight_layout();save(fig,r.heat_id+'_高证据点修复敏感性')
def review():
 r=json.loads((P/'config/native_review_vetoes.json').read_text())['vetoes'][0];y,t,w=load(r['heat_id']);i=r['row'];j=r['column'];sl=slice(j-10,j+11);fig,axs=plt.subplots(1,2,figsize=(12,4));rows=[]
 for k in range(i-2,i+3):
  axs[0].plot(w[sl],y[k,sl],label=f't={t[k]} s',linewidth=2 if k==i else 1);ix=j-5+int(np.argmax(y[k,j-5:j+6]));rows.append(dict(heat_id=r['heat_id'],row=k,time_s=t[k],target_wavelength=w[j],target_record=y[k,j],nearby_peak_wavelength=w[ix],nearby_peak_record=y[k,ix],detected_current=k==i))
 axs[0].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title='五个相邻时刻均存在的真实记录谱线');axs[0].legend();im=axs[1].imshow(y[i-2:i+3,sl],aspect='auto',origin='lower',extent=[w[j-10],w[j+10],t[i-2],t[i+2]],cmap='viridis');axs[1].set(xlabel='波长 / nm',ylabel='时间 / s',title='候选位置邻通道共同增强：降级保留');fig.colorbar(im,ax=axs[1],label='counts');fig.tight_layout();save(fig,'持续谱线候选复核保留',True);pd.DataFrame(rows).to_csv(P/'results/persistent_line_veto_neighborhood.csv',index=False)

if __name__=='__main__':
 if len(sys.argv)==1:main()
 elif sys.argv[1]=='sensitivity':sensitivity()
 else:review()
