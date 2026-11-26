"""Paper figures derived only from saved observations and locked experiments."""
from ctx import *
import pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.utils.extmath import randomized_svd
plt.rcParams.update({'font.family':'Noto Serif CJK SC','axes.unicode_minus':False,'font.size':10,'savefig.dpi':180})
M=P/'figures/main';A=P/'figures/appendix'
def save(fig,name,main=True):
 d=M if main else A;fig.savefig(d/(name+'.png'),bbox_inches='tight');fig.savefig(d/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def main():
 stats=[];energy=[];peaks=[];overview=plt.figure(figsize=(14,22));overview.subplots_adjust(hspace=.5,wspace=.05)
 for z,h in enumerate(HEATS):
  y,t,w=load(h);d=np.load(P/f'data/{h}_interface.npz');old=np.load(V2/f'data/{h}_model.npz');g=old['grid'];mask=np.load(V2/f'data/{h}_masks.npz')['valid'];D=d['reconstructed'];s=np.linalg.svd(y[::8],compute_uv=False);e=np.cumsum(s*s)/np.sum(s*s)
  for k in [1,2,4,8,16,32]:energy.append(dict(heat_id=h,rank=k,energy=float(e[k-1]),time_stride=8,source='raw recorded including top codes; descriptive only'))
  stats.append(dict(heat_id=h,frames=len(t),channels=len(w),time_min=t[0],time_max=t[-1],wavelength_min=w[0],wavelength_max=w[-1],sat=int((y==65535).sum()),candidate_frames=int(pd.read_csv(P/f'results/{h}_frame_diagnostics.csv')['candidate'].sum()),intensity_median=float(np.median(y)),intensity_max=float(y.max())))
  diff=abs(g.max(1)-y.max(1))/np.maximum(y.max(1),1);i=int(np.argmax(diff));j=int(np.argmax(g[i]));peaks.append(dict(heat_id=h,time=float(t[i]),row=i,wavelength=float(w[j]),V2_peak_relative_change=float(diff[i]),raw_at_V2_peak=float(y[i,j]),V2_at_peak=float(g[i,j]),old_valid_at_peak=bool(mask[i,j]),V3_at_peak=float(D[i,j]),V3_peak_relative_change=float(np.max(abs(D.max(1)-y.max(1))/np.maximum(y.max(1),1)))))
  ii=np.linspace(0,len(t)-1,60).astype(int);jj=np.linspace(0,len(w)-1,90).astype(int);W,T=np.meshgrid(w[jj],t[ii])
  for fig,ax in [(plt.figure(figsize=(9,6)),None),(overview,None)]:
   ax=fig.add_subplot(1,1,1,projection='3d') if fig is not overview else fig.add_subplot(5,2,z+1,projection='3d');ax.plot_surface(W,T,D[np.ix_(ii,jj)],cmap='viridis',linewidth=0,rasterized=True);ax.set(xlabel='波长 / nm',ylabel='时间 / s',zlabel='记录强度 / counts',title=h);ax.view_init(25,-125)
   if fig is not overview:save(fig,h+'_连续光谱曲面',False)
  fig,ax=plt.subplots(figsize=(10,3));ax.imshow(d['quality'],aspect='auto',origin='lower',extent=[w[0],w[-1],t[0],t[-1]],cmap='viridis');ax.set(xlabel='波长 / nm',ylabel='时间 / s',title=h+' 质量位掩码（多标记可叠加）');save(fig,h+'_质量标记',False)
 save(overview,'十炉连续三维光谱总览');pd.DataFrame(stats).to_csv(P/'results/data_statistics.csv',index=False);pd.DataFrame(energy).to_csv(P/'results/singular_energy.csv',index=False);pd.DataFrame(peaks).to_csv(P/'results/peak_origin_audit.csv',index=False)
 fig,axs=plt.subplots(1,3,figsize=(15,4));sdf=pd.DataFrame(stats);x=np.arange(10)
 for ax,col,title in zip(axs,['frames','sat','candidate_frames'],['采样帧数','疑似截顶记录数','帧级统计候选数']):ax.bar(x,sdf[col],color='#277f8e');ax.set_xticks(x,HEATS,rotation=70);ax.set_title(title)
 save(fig,'十炉数据特征统计')
 fig,ax=plt.subplots(figsize=(8,4))
 for h,q in pd.DataFrame(energy).groupby('heat_id'):ax.plot(q['rank'],q.energy,marker='.',label=h)
 ax.set(xlabel='秩',ylabel='累积奇异值平方占比',ylim=(.8,1.005),title='原始记录矩阵描述性低秩分析（每8帧采样）');ax.legend(ncol=2,fontsize=8);save(fig,'奇异值能量贡献率')
 q=pd.read_csv(P/'results/model_confirmation.csv');means=q.groupby(['scenario','model']).rmse.mean().unstack();names=['V1','V2','V3-C'];fig,ax=plt.subplots(figsize=(9,4))
 for k,n in enumerate(names):ax.bar(np.arange(3)+(k-1)*.24,means[n],.24,label=n)
 ax.set_xticks(np.arange(3),means.index);ax.set(ylabel='逐炉宏平均 RMSE / counts',title='新遮蔽确认：缺失位置预测（观测复制不计分）');ax.legend();save(fig,'V1_V2_V3配对比较')
 q=pd.read_csv(P/'results/continuous_confirmation.csv');print('continuous columns',q.columns.tolist(),flush=True);means=q.groupby(['scenario','model']).rmse.mean().unstack();fig,ax=plt.subplots(figsize=(9,4))
 for k,n in enumerate(names):ax.bar(np.arange(len(means))+(k-1)*.24,means[n],.24,label=n)
 ax.set_xticks(np.arange(len(means)),means.index);ax.set(ylabel='留出真实坐标 RMSE / counts',title='真实坐标降采样连续重构');ax.legend();save(fig,'连续插值精度比较')
 h='6C02058';y,t,w=load(h);old=np.load(V2/f'data/{h}_model.npz')['grid'];d=np.load(P/f'data/{h}_interface.npz')['reconstructed'];i=int(np.argmin(abs(t-445.5)));j=int(np.argmax(old[i]));fig,axs=plt.subplots(1,2,figsize=(13,4))
 for a,ran in zip(axs,[slice(None),slice(max(0,j-30),j+31)]):
  a.plot(w[ran],y[i,ran],label='原始记录',color='black');a.plot(w[ran],np.load(V1/f'data/reconstructed/{h}_on_observed_grid.npz')['intensity'][i,ran],label='V1历史输出',alpha=.6);a.plot(w[ran],old[i,ran],label='V2潜变量输出',alpha=.7);a.plot(w[ran],d[i,ran],label='V3记录尺度',ls='--');a.set(xlabel='波长 / nm',ylabel='记录强度 / counts',title=f'{h}, t={t[i]} s');a.legend()
 save(fig,'典型谱峰保真对比')
 fig,axs=plt.subplots(2,2,figsize=(12,8))
 for ax,h in zip(axs.flat,CFG['stability_heats']):
  for scope in ['native','held']:
   hist=json.loads((P/f'results/{h}_{scope}_solver.json').read_text())['history'];ax.semilogy([r['iteration'] for r in hist],[r['relative_change'] for r in hist],label=scope)
  ax.axhline(1e-4,color='gray',ls='--');ax.set(xlabel='迭代次数',ylabel='相对迭代变化率',title=h);ax.legend()
 save(fig,'代表性求解器收敛')
 fig,ax=plt.subplots(figsize=(12,4));ax.axis('off');boxes=[('原始记录\n10炉 / 13,589帧',.03),('分级诊断\n1,129候选；0确认整帧',.27),('保真与边界投影\n仅补缺失 / 无负光强',.52),('双线性连续曲面\n截顶标记；条件估计另存',.77)]
 for text,x in boxes:ax.text(x,.5,text,transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.8',fc='#e8f2ef',ec='#277f8e'),va='center',fontsize=11)
 for x in [.22,.47,.72]:ax.annotate('',xy=(x+.04,.5),xytext=(x,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
 ax.text(.03,.15,'低秩缺失预测：固定秩4、Huber尺度0.003；真实仪器故障无标签，不自动删统计候选。',transform=ax.transAxes);save(fig,'第一问整体技术路线')
 # One actual fixed held-out panel, used only to illustrate already locked experiment.
 for h,stage,y,t,w,_,_ in panels('model'):
  if h=='6C02058' and stage==1:
   hidden=dict(hides(y.shape,9307))['random_0.05'];v=y.copy();v[hidden]=np.nan;from preprocess import apply;from model3 import project_record
   pp=apply(v,t,w,'P0');r=legacy.fit_variant('B2',v,t,w,pp['valid'],rank=4,eta=.003,max_iter=100) if False else None
   # Visible-only linear is explicitly an illustration baseline, not an additional V3 score.
   filled=fill_training(v,t,w);fig,axs=plt.subplots(1,3,figsize=(14,4))
   for ax,img,title in zip(axs,[y,np.where(hidden,np.nan,y),filled],['原始记录参考','固定随机遮蔽','双向线性补全示例']):ax.imshow(img,origin='lower',aspect='auto',extent=[w[0],w[-1],t[0],t[-1]],cmap='viridis');ax.set(title=title,xlabel='波长 / nm',ylabel='时间 / s')
   save(fig,'典型缺失识别与补全')
if __name__=='__main__':main()
