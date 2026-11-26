from core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Noto Serif CJK SC','font.size':10,'axes.unicode_minus':False})
A=P/'figures/appendix';M=P/'figures/main'
def save(fig,name,main=False):
 d=M if main else A;fig.savefig(d/(name+'.png'),dpi=170,bbox_inches='tight');p=d/(name+'.svg');fig.savefig(p,bbox_inches='tight');p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n');plt.close(fig)

def main():
 assert (P/'results/FINAL_COMPLETE.json').exists();repairs=pd.read_csv(P/'results/actual_repairs.csv');overview=plt.figure(figsize=(14,22));overview.subplots_adjust(hspace=.5,wspace=.05)
 for z,h in enumerate(HEATS):
  d=np.load(P/f'data/{h}_interface.npz');y=d['raw'];D=d['reconstructed'];t=d['time'];w=d['wavelength'];sub=repairs[repairs.heat_id==h];r=sub.iloc[int(np.argmax(abs(sub.difference.to_numpy())))];i=int(r.row);j=int(r.column);a=max(0,i-5);b=min(len(t),i+6);c=max(0,j-10);e=min(len(w),j+11)
  fig,axs=plt.subplots(2,3,figsize=(15,8));axs[0,0].plot(w,y[i],label='原始记录',color='black',lw=1);axs[0,0].plot(w,D[i],label='统计修复',ls='--');axs[0,0].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title=f't={t[i]} s');axs[0,0].legend()
  axs[0,1].plot(w[c:e],y[i,c:e],label='原始');axs[0,1].plot(w[c:e],D[i,c:e],ls='--',label='修复估计');axs[0,1].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title='典型异常局部谱');axs[0,1].legend()
  axs[0,2].plot(t[a:b],y[a:b,j],label='原始');axs[0,2].plot(t[a:b],D[a:b,j],ls='--',label='修复');axs[0,2].set(xlabel='时间 / s',ylabel='记录强度 / counts',title=f'λ={w[j]} nm');axs[0,2].legend()
  mn=y[a:b,c:e].min();mx=y[a:b,c:e].max();lim=max(abs(D[a:b,c:e]-y[a:b,c:e]).max(),1)
  for ax,img,title in zip(axs[1],[y[a:b,c:e],D[a:b,c:e],D[a:b,c:e]-y[a:b,c:e]],['原始局部热图','修复后局部热图','修复−原始']):
   diff=title=='修复−原始';im=ax.pcolormesh(w[c:e],t[a:b],img,shading='nearest',cmap='coolwarm' if diff else 'viridis',vmin=-lim if diff else mn,vmax=lim if diff else mx,rasterized=True);ax.set(xlabel='波长 / nm',ylabel='时间 / s',title=title);fig.colorbar(im,ax=ax,label='counts')
  fig.suptitle(h+'：多证据统计异常的实际修复');fig.tight_layout();save(fig,h+'_异常修复曲线与局部热图')
  fig,axs=plt.subplots(1,3,figsize=(15,4));ii=np.arange(0,len(t),3);jj=np.arange(0,len(w),3)
  for ax,grid,title in zip(axs[:2],[y,D],['原始全谱热图','处理后全谱热图']):im=ax.pcolormesh(w[jj],t[ii],grid[np.ix_(ii,jj)],shading='nearest',cmap='viridis',vmin=0,vmax=65535,rasterized=True);ax.set(xlabel='波长 / nm',ylabel='时间 / s',title=title);fig.colorbar(im,ax=ax)
  pnt=axs[2].scatter(sub.wavelength_nm,sub.time_s,c=sub.difference,cmap='coolwarm',s=8,vmin=-max(abs(sub.difference)),vmax=max(abs(sub.difference)));axs[2].set(xlabel='波长 / nm',ylabel='时间 / s',xlim=(w[0],w[-1]),ylim=(t[0],t[-1]),title=f'全矩阵差值位置（{len(sub)}点）');fig.colorbar(pnt,ax=axs[2],label='修复−原始 / counts');fig.tight_layout();save(fig,h+'_全谱处理与差值')
  ii=np.linspace(0,len(t)-1,60).astype(int);jj=np.linspace(0,len(w)-1,90).astype(int);W,T=np.meshgrid(w[jj],t[ii])
  for fig in [plt.figure(figsize=(9,6)),overview]:
   ax=fig.add_subplot(1,1,1,projection='3d') if fig is not overview else fig.add_subplot(5,2,z+1,projection='3d');ax.plot_surface(W,T,D[np.ix_(ii,jj)],cmap='viridis',linewidth=0,rasterized=True);ax.set(title=h,xlabel='波长 / nm',ylabel='时间 / s',zlabel='记录强度 / counts');ax.view_init(25,-125)
   if fig is not overview:save(fig,h+'_最终连续三维光谱')
  # Raw fast coherent transitions remain; no whole-frame deletion.
  frame=pd.read_csv(V3/f'results/{h}_frame_diagnostics.csv');k=int(np.nanargmax(frame.integral_relative_first_difference.abs().to_numpy()));sel=np.arange(max(0,k-2),min(len(t),k+3));fig,axs=plt.subplots(1,2,figsize=(12,4));axs[0].plot(t,np.trapezoid(y,w,axis=1),label='原始');axs[0].plot(t,np.trapezoid(D,w,axis=1),ls='--',label='修复后');axs[0].set(xlabel='时间 / s',ylabel='积分记录强度 / counts·nm',title='整帧积分及快变保留');axs[0].legend()
  for k in sel:axs[1].plot(w,D[k],label=f'{t[k]}s')
  axs[1].set(xlabel='波长 / nm',ylabel='记录强度 / counts',title='快变邻帧的光谱形状');axs[1].legend(fontsize=8);fig.tight_layout();save(fig,h+'_帧级快变诊断')
 save(overview,'十炉最终连续光谱总览',True)
 # Four specifically reviewed real points, not merely the most visually striking cases.
 review=pd.read_csv(V31/'results/high_confidence_statistical_points.csv');fig,axs=plt.subplots(2,2,figsize=(12,8))
 for ax,r in zip(axs.flat,review.itertuples()):
  d=np.load(P/f'data/{r.heat_id}_interface.npz');j=r.column;i=r.row;ix=slice(j-6,j+7);ax.plot(d['wavelength'][ix],d['raw'][i,ix],label='原始');ax.plot(d['wavelength'][ix],d['reconstructed'][i,ix],ls='--',label='实际修复估计');ax.set(title=f'{r.heat_id} t={r.time_s}s',xlabel='波长 / nm',ylabel='counts');ax.legend()
 fig.tight_layout();save(fig,'四个高证据尖峰实际修复',True)
 # Unified comparison and per-band behavior.
 q=pd.read_csv(P/'results/unified_comparison.csv');q=q[q.split=='confirmation'];names=['A_P0','B_V1_MAD','C_V31','D_adaptive4'];q=q.set_index('strategy').loc[names];fig,axs=plt.subplots(1,2,figsize=(12,4));q[['precision','recall','f1']].plot.bar(ax=axs[0],color=['#277f8e','#6da881','#d99138']);q[['pollution_rmse','unpolluted_rmse']].plot.bar(ax=axs[1],color=['#277f8e','#d99138']);axs[0].set(ylabel='绝对检测指标',xlabel='策略');axs[1].set(ylabel='RMSE / counts',xlabel='策略');fig.tight_layout();save(fig,'统一端到端异常处理比较',True)
 z=pd.read_csv(P/'results/confirmation_per_heat_band.csv');fig,axs=plt.subplots(1,2,figsize=(12,4));z.groupby(['band','strategy']).pollution_rmse.mean().unstack()[names].plot.bar(ax=axs[0]);z.groupby(['band','strategy']).f1.mean().unstack()[names].plot.bar(ax=axs[1]);axs[0].set(ylabel='污染RMSE / counts',xlabel='波长区域');axs[1].set(ylabel='案例宏F1',xlabel='波长区域');fig.tight_layout();save(fig,'波长区域稳定性',True)
 sv=pd.read_csv(V3/'results/singular_energy.csv');fig,ax=plt.subplots(figsize=(8,4))
 for h,d in sv.groupby('heat_id'):ax.plot(d['rank'],d.energy,label=h,marker='.')
 ax.set(xlabel='秩',ylabel='累计奇异值平方占比',ylim=(.8,1.001));ax.legend(ncol=2,fontsize=8);save(fig,'低秩能量贡献',True)
 n=pd.read_csv(P/'results/final_summary.csv');fig,axs=plt.subplots(1,2,figsize=(12,4));axs[0].bar(n.heat_id,n.repair_points,color='#277f8e');axs[1].bar(n.heat_id,n.repair_fraction*100,color='#6da881');axs[0].set(ylabel='实际修复点数');axs[1].set(ylabel='修复比例 / %')
 for ax in axs:ax.tick_params(axis='x',rotation=60)
 fig.tight_layout();save(fig,'十炉实际修复统计',True)
 fig,ax=plt.subplots(figsize=(13,4));ax.axis('off');total=int(n.repair_points.sum());boxes=[(.01,'真实坐标与质量审计\n十炉原始光谱'),(.26,'局部尺度＋多证据评分\n区分瞬态异常与持续谱线'),(.53,f'固定秩4 Huber补全\n保真/非负；实际修复{total}点'),(.80,'双线性连续曲面\n饱和记录与条件估计分开')]
 for x,text in boxes:ax.text(x,.5,text,transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.7',fc='#e5f0ed',ec='#277f8e'),va='center',fontsize=10)
 for x in [.22,.49,.76]:ax.annotate('',xy=(x+.03,.5),xytext=(x,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
 ax.text(.01,.12,'共同协议端到端比较 → 开发锁定 → 新区域确认 → 全炉实际修复 → 连续模型检验',transform=ax.transAxes);save(fig,'第一问技术路线',True)
 fig,axs=plt.subplots(2,2,figsize=(12,8))
 for ax,h in zip(axs.flat,['6C02057','6C02058','6C02107','6C02108']):
  hist=json.loads((P/f'results/{h}_solver.json').read_text())['history'];ax.semilogy([x['iteration'] for x in hist],[x['relative_change'] for x in hist]);ax.axhline(1e-4,ls='--',color='gray');ax.set(title=h,xlabel='迭代',ylabel='相对变化')
 fig.tight_layout();save(fig,'最终求解收敛',True)
if __name__=='__main__':main()
