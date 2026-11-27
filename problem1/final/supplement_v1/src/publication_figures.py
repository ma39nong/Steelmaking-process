"""Seven publication figures and evidence appendices from frozen and supplementary data."""
from pathlib import Path
import sys,json,hashlib,re,argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,cm
from matplotlib.patches import Rectangle
from matplotlib.text import Text
S=Path(__file__).resolve().parents[1];F=S.parent;P=F.parent;O=F/'figures/publication_v3'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.2,'axes.labelsize':8.2,'axes.titlesize':8.5,'legend.fontsize':7.6,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1.1,'svg.fonttype':'none','pdf.fonttype':42,'mathtext.fontset':'dejavusans','savefig.dpi':600})
RAW='#344c5c';REPAIR='#178b82';ANOM='#cc624b';FALL='#867399';GROUP={'6B':'#365f82','6C':'#c28a45'}
HEATS=pd.read_csv(F/'results/final_summary.csv').heat_id.tolist();MANIFEST=[];SAMPLING=[]

def load(h):
 with np.load(F/f'data/{h}_interface.npz') as z:return {k:z[k] for k in ['raw','reconstructed','time','wavelength','statistical_repair','suspected_saturation','latent','numerical_fallback']}

def panel(ax,label,title):ax.set_title(f'({label}) {title}',loc='left',pad=6)

def save(fig,name,section='main'):
 path=O/section;fig.canvas.draw();labels=[x.get_text() for x in fig.findobj(match=Text)]
 assert not any(re.search(r'[\u3400-\u9fff]',x) for x in labels),name
 for fmt in ['png','svg','pdf']:
  q=path/f'{name}.{fmt}';fig.savefig(q,dpi=600,bbox_inches='tight',metadata={'Creator':'Steelmaking-process publication v3'})
  if fmt=='svg':
   text=q.read_text();assert not re.search(r'[\u3400-\u9fff]',text);q.write_text('\n'.join(x.rstrip() for x in text.splitlines())+'\n')
  MANIFEST.append(dict(figure=name,section=section,format=fmt,path=q.relative_to(F).as_posix(),bytes=q.stat().st_size,sha256=hashlib.sha256(q.read_bytes()).hexdigest()))
 plt.close(fig);print('rendered',name,flush=True)

def F01():
 g=pd.read_csv(S/'results/common_grid_group_summaries.csv');dark=pd.read_csv(S/'results/dark_bright_focus_spectra.csv').query('fraction==0.15');sv=pd.read_csv(P/'v3/results/singular_energy.csv')
 fig=plt.figure(figsize=(6.8,5.4),layout='constrained');gs=fig.add_gridspec(2,3,height_ratios=[1.05,1]);ax=fig.add_subplot(gs[0,:])
 for group,d in g.groupby('group'):
  ax.plot(d.wavelength_nm,d['median'],color=GROUP[group],label=f'{group}: equal-heat mean median');ax.fill_between(d.wavelength_nm,d.q10,d.q90,color=GROUP[group],alpha=.12)
 for wave,label in [(589.3,'Na I ref.'),(760,'O2 A ref.'),(768,'K I ref.')]:ax.axvline(wave,color='#a1a6aa',ls=':',lw=.7)
 ax.text(589.3,.95,'Na I ref.',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7);ax.text(764,.84,'O2 / K I refs.',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7);ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)',xlim=(340.535,1023.303));ax.legend(loc='upper left');panel(ax,'a','Group spectra and mean within-heat 10-90% bands')
 ax=fig.add_subplot(gs[1,0]);first={}
 for h,d in sv.groupby('heat_id'):
  col=GROUP[h[:2]];ax.semilogy(d['rank'],1-d.energy,color=col,alpha=.45,marker='.',ms=3,label=h[:2] if h[:2] not in first else None);first[h[:2]]=1
 ax.axvline(4,color=ANOM,ls='--',lw=.8);ax.set(xlabel='Rank r',ylabel='Unexplained energy 1-E(r)',xlim=(1,16),xticks=[1,4,8,16]);ax.legend();panel(ax,'b','Raw sampled SVD')
 for k,(prefix,wave) in enumerate([('6B',835.887),('6C',383.293)],1):
  ax=fig.add_subplot(gs[1,k]);d=dark[dark.heat_id.str.startswith(prefix)].groupby('wavelength_nm')[['dark','bright']].mean().reset_index();ax.plot(d.wavelength_nm,d.dark,color=RAW,label='Dark 15%');ax.plot(d.wavelength_nm,d.bright,color=REPAIR,label='Bright 15%');ax.axvline(wave,color=ANOM,ls=':',lw=.8);ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');panel(ax,chr(98+k),f'{prefix}: {wave:.3f} nm');ax.legend(fontsize=7)
 save(fig,'F01_data_structure')

def case(h,i,j):
 z=load(h);points=pd.read_csv(S/'results/repair_reliability_points.csv');r=points[(points.heat_id==h)&(points.row==i)&(points.column==j)].iloc[0];return z,r

def F02():
 z,r=case('6B03205',29,466);y=z['raw'];D=z['reconstructed'];t=z['time'];w=z['wavelength'];i,j=29,466;a,b=i-4,i+5;c,e=j-7,j+8
 fig=plt.figure(figsize=(6.8,6.1),layout='constrained');gs=fig.add_gridspec(3,3,height_ratios=[1,1,.65]);ax=fig.add_subplot(gs[0,0]);im=ax.pcolormesh(w[c:e],t[a:b],y[a:b,c:e],cmap='viridis',shading='nearest',rasterized=True);ax.scatter(w[j],t[i],facecolors='none',edgecolors=ANOM,s=38);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)');panel(ax,'a','Recorded input');fig.colorbar(im,ax=ax,label='counts',shrink=.7)
 ax=fig.add_subplot(gs[0,1]);ax.plot(t[a:b],y[a:b,j],'.-',color=RAW,label='Raw');ax.scatter(t[i],r.time_reference,color=REPAIR,marker='x',s=32,label='Time reference');ax.set(xlabel='Time (s)',ylabel='Intensity (counts)');panel(ax,'b',fr'Time evidence: $z_t$={r.time_z:.2f}');ax.legend(fontsize=6.9)
 ax=fig.add_subplot(gs[0,2]);ax.plot(w[c:e],y[i,c:e],'.-',color=RAW,label='Raw');ax.scatter(w[j],r.wavelength_reference,color=REPAIR,marker='x',s=32,label='Wave reference');ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');panel(ax,'c',fr'Wave evidence: $z_\lambda$={r.wave_z:.2f}');ax.legend(fontsize=6.9)
 ax=fig.add_subplot(gs[1,0]);vals=[r.time_z,r.wave_z,r.lowrank_z];ax.hlines([0,1,2],0,vals,color='#c8d2d5',lw=2);ax.scatter(vals,[0,1,2],c=[RAW,REPAIR,FALL],s=24);ax.set(yticks=[0,1,2],yticklabels=[r'$z_t$',r'$z_\lambda$',r'$z_L$'],xlabel='Standardized evidence');cons=r.score/np.sqrt(r.time_z*(.5*r.wave_z/(1+r.persistence)+.5*r.lowrank_z));ax.set_ylim(-1.1,2.3);ax.text(.02,.04,f'P={r.persistence}; c={cons:.3f}\nS={r.score:.2f} > 4',transform=ax.transAxes,fontsize=8);panel(ax,'d','Evidence fusion')
 with np.load(F/'data/6B03205_masks.npz') as m:score=m['score'][a:b,c:e]
 ax=fig.add_subplot(gs[1,1]);im=ax.pcolormesh(w[c:e],t[a:b],score,cmap='magma',shading='nearest',rasterized=True);ax.scatter(w[j],t[i],s=40,facecolors='none',edgecolors='white');ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)');panel(ax,'e','Score and selected mask');fig.colorbar(im,ax=ax,label='Score',shrink=.7)
 ax=fig.add_subplot(gs[1,2]);ax.plot(w[c:e],y[i,c:e],color=RAW,label='Raw');ax.plot(w[c:e],D[i,c:e],color=REPAIR,label='Repaired');ax.scatter(w[j],D[i,j],color=REPAIR,s=22);ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');ax.legend(fontsize=7);panel(ax,'f',f'{y[i,j]:.0f} to {D[i,j]:.2f} counts')
 ax=fig.add_subplot(gs[2,:]);ax.axis('off');ax.text(.01,.80,r'$S=c\sqrt{z_t\{z_\lambda/[2(1+P)]+z_L/2\}}$',fontsize=11);ax.text(.01,.40,r'$\mathrm{rank}(L)\leq4,\quad\min\sum_{\Omega_v}h_{0.003}(L-Y/U),\quad D=\mathrm{fidelity+support}(UL)$',fontsize=9);ax.text(.01,.05,'Score >= 4; connected span <= 3 frames x 7 channels; top codes handled separately.',fontsize=8)
 save(fig,'F02_model_mechanism')

def F03():
 review=pd.read_csv(F/'results/four_reviewed_actual_repairs.csv');fig=plt.figure(figsize=(6.8,6.5),layout='constrained');gs=fig.add_gridspec(3,4,height_ratios=[1,1,1]);z,r=case('6B03205',29,466);y=z['raw'];D=z['reconstructed'];w=z['wavelength'];t=z['time']
 ax=fig.add_subplot(gs[0,:2]);ii=np.unique(np.r_[np.arange(0,len(t),4),29]);jj=np.unique(np.r_[np.arange(0,len(w),4),466]);ax.pcolormesh(w[jj],t[ii],y[np.ix_(ii,jj)],cmap='viridis',vmin=0,vmax=65535,rasterized=True,shading='nearest');ax.scatter(w[466],t[29],facecolors='none',edgecolors=ANOM,s=40);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)');panel(ax,'a','6B03205: case location')
 ax=fig.add_subplot(gs[0,2:]);ix=slice(25,34);ax.plot(t[ix],y[ix,466],'.-',color=RAW,label='Raw');ax.plot(t[ix],D[ix,466],'.--',color=REPAIR,label='Repaired');ax.scatter(t[29],r.time_reference,marker='x',color=ANOM,s=35,label='Time reference');ax.set(xlabel='Time (s)',ylabel='Intensity (counts)');ax.legend(fontsize=7);panel(ax,'b',f'Temporal evidence; S={r.score:.2f}')
 for k,row in enumerate(review.itertuples()):
  z=load(row.heat_id);w=z['wavelength'];j=row.column;i=row.row;ix=slice(j-5,j+6);ax=fig.add_subplot(gs[1,k]);ax.plot(w[ix],z['raw'][i,ix],color=RAW,label='Raw');ax.plot(w[ix],z['reconstructed'][i,ix],color=REPAIR,ls='--',label='Repaired');ax.scatter(w[j],row.repaired,color=REPAIR,s=18);ax.set(xlabel='Wavelength (nm)',ylabel='counts');ax.tick_params(axis='x',labelrotation=25);panel(ax,chr(99+k),row.heat_id);ax.text(.03,.93,f'{row.original:.0f} -> {row.repaired:.1f}\nS={row.score:.2f}',transform=ax.transAxes,va='top',fontsize=7.5,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.4))
 z,r=case('6C02108',608,1211);w=z['wavelength'];t=z['time'];y=z['raw'];j=1211;i=608
 ax=fig.add_subplot(gs[2,:2]);ix=slice(j-8,j+9)
 for q in range(i-2,i+3):ax.plot(w[ix],y[q,ix],color=RAW if q==i else '#bbc7cd',lw=1.5 if q==i else .7)
 ax.scatter(w[j],y[i,j],color=REPAIR,marker='D',s=20);ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');panel(ax,'g','Retained persistent feature')
 ax=fig.add_subplot(gs[2,2:]);ix=slice(i-5,i+6);ax.plot(t[ix],y[ix,j],'.-',color=RAW,label='Target channel');ax.plot(t[ix],.5*(y[ix,j-1]+y[ix,j+1]),ls='--',color=REPAIR,label='Neighbor mean');ax.set(xlabel='Time (s)',ylabel='Intensity (counts)');ax.legend(fontsize=7);panel(ax,'h','6C02108: 766.413 nm')
 save(fig,'F03_local_evidence')

def F04():
 d=pd.read_csv(P/'v2/results/interpolation_independent_macro.csv',header=[0,1],index_col=[0,1]);fig=plt.figure(figsize=(6.8,4.9),layout='constrained');gs=fig.add_gridspec(2,2);methods=['linear','pchip','cubic'];cols=[REPAIR,RAW,ANOM];labels=['Bilinear','PCHIP','Cubic / tensor']
 for col,(kind,xvals,scenes) in enumerate([('Time',[1,1.5,2],['time_stride_2','time_stride_3','time_stride_4']),('Wave',[2,3],['wave_stride_2','wave_stride_3'])]):
  ax=fig.add_subplot(gs[0,col])
  for m,c,label in zip(methods,cols,labels):ax.plot(xvals,[d.loc[(m,s),('rmse','mean')] for s in scenes],marker='o',color=c,ms=4,label=label)
  ax.set(xlabel='Retained interval (s)' if kind=='Time' else 'Channel stride',ylabel='Held-out RMSE (counts)',xticks=xvals);ax.legend(fontsize=7);panel(ax,chr(97+col),kind+' downsampling')
 z=load('6B03205');i,j=29,466;t=z['time'];w=z['wavelength'];D=z['reconstructed'];a,b=.35,.6;corner=np.array([D[i,j],D[i+1,j],D[i,j+1],D[i+1,j+1]]);weights=np.array([(1-a)*(1-b),a*(1-b),(1-a)*b,a*b]);query=weights@corner
 ax=fig.add_subplot(gs[1,0]);coords=[(0,0),(1,0),(0,1),(1,1)]
 for (x,y),v,wt in zip(coords,corner,weights):ax.scatter(x,y,color=RAW,s=25);ax.annotate(f'{v:.1f}\nw={wt:.2f}',(x,y),xytext=(4,5 if y==0 else -25),textcoords='offset points',fontsize=7)
 ax.add_patch(Rectangle((0,0),1,1,fill=False,ec='#899ca5',lw=.7));ax.scatter(a,b,color=REPAIR,s=25);ax.set(xlim=(-.1,1.4),ylim=(-.25,1.3),xticks=[0,1],xticklabels=[f'{t[i]:.1f}',f'{t[i+1]:.1f}'],yticks=[0,1],yticklabels=[f'{w[j]:.3f}',f'{w[j+1]:.3f}'],xlabel='Time (s)',ylabel='Wavelength (nm)');panel(ax,'c','Actual-coordinate bilinear cell')
 ax=fig.add_subplot(gs[1,1]);ax.axis('off');ax.text(0,.87,r'$\widehat F=\sum_{k=1}^{4}w_kD_k$',fontsize=12);ax.text(0,.61,r'$w_k\geq0,\quad\sum w_k=1$',fontsize=11);ax.text(0,.36,f'Query value: {query:.2f} counts\nCorner range: {corner.min():.1f} - {corner.max():.1f}\nContinuous, nonnegative, within corner range.',fontsize=8.2);panel(ax,'d','Convex-combination properties')
 pd.DataFrame(dict(corner=['00','10','01','11'],value=corner,weight=weights,time=[t[i],t[i+1],t[i],t[i+1]],wavelength=[w[j],w[j],w[j+1],w[j+1]])).to_csv(O/'bilinear_cell_values.csv',index=False)
 save(fig,'F04_continuous_reconstruction')

def F05():
 u=pd.read_csv(F/'results/unified_comparison.csv').query('split=="confirmation"').set_index('strategy');keys=['A_P0','B_V1_MAD','C_V31','D_adaptive4','D_adaptive6'];labels=['Retain all','MAD','Conservative','Adaptive 4','Adaptive 6'];cols=['#8f999e',RAW,FALL,REPAIR,ANOM];fig,axs=plt.subplots(2,2,figsize=(6.8,5.7),layout='constrained');ax=axs[0,0]
 for f in [.2,.4,.6,.8]:
  rr=np.linspace(max(f/2+.001,.05),1,300);pp=f*rr/(2*rr-f);ok=pp<=1;ax.plot(rr[ok],pp[ok],color='#c9cdd0',ls='--',lw=.65);ax.text(.98,f/(2-f)+.015,f'{f:.1f}',fontsize=6.8,ha='right',color='#737c81')
 for key,label,c in zip(keys[1:],labels[1:],cols[1:]):r=u.loc[key];ax.scatter(r.recall,r.precision,color=c,s=25);ax.annotate(label,(r.recall,r.precision),xytext=(4,-11 if key=='C_V31' else 5),textcoords='offset points',fontsize=7)
 ax.set(xlabel='Recall',ylabel='Precision',xlim=(0,1.03),ylim=(0,1.10));panel(ax,'a','Finite strategies and F1 contours')
 ax=axs[0,1]
 for key,label,c in zip(keys,labels,cols):r=u.loc[key];ax.scatter(r.unpolluted_rmse,r.pollution_rmse,color=c,s=25);ax.annotate(label,(r.unpolluted_rmse,r.pollution_rmse),xytext=(4,3),textcoords='offset points',fontsize=7)
 ax.set(xlabel='Unpolluted RMSE (counts)',ylabel='Polluted RMSE (counts)',xlim=(-8,210),ylim=(0,float(u.loc[keys,'pollution_rmse'].max())*1.12));panel(ax,'b','Recovery benefit and cost (300 steps)')
 (O/'F05_axis_bounds.json').write_text(json.dumps({'x_limits':list(ax.get_xlim()),'y_limits':list(ax.get_ylim()),'all_strategy_points_visible':bool(all(ax.get_xlim()[0]<=u.loc[k,'unpolluted_rmse']<=ax.get_xlim()[1] and ax.get_ylim()[0]<=u.loc[k,'pollution_rmse']<=ax.get_ylim()[1] for k in keys)),'retain_all_pollution_rmse':float(u.loc['A_P0','pollution_rmse'])},indent=2)+'\n')
 b=pd.read_csv(F/'results/confirmation_band_summary.csv');mat=np.array([[b[(b.band==band)&(b.strategy==key)].f1.iloc[0] for key in keys] for band in ['blue','middle','red']]);ax=axs[1,0];im=ax.imshow(mat,cmap='viridis',vmin=0,vmax=1,aspect='auto')
 for (i,j),v in np.ndenumerate(mat):ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7.5,color='white' if v<.5 else 'black')
 ax.set(yticks=[0,1,2],yticklabels=['Blue','Middle','Red'],xticks=range(5),xticklabels=['Retain','MAD','Cons.','S=4','S=6']);ax.tick_params(axis='x',rotation=25);panel(ax,'c','Mean case-wise F1 by region')
 ax=axs[1,1];exact=pd.read_csv(F/'results/exact_final_pipeline_confirmation.csv').groupby('heat_id').mean(numeric_only=True).mean()
 for yy,(before,after,name) in enumerate([(exact.before_pollution_rmse,exact.pollution_rmse,'Polluted'),(exact.before_overall_rmse,exact.overall_rmse,'Overall')]):
  ax.plot([before,after],[yy,yy],color='#bbc8cd',lw=2);ax.scatter(before,yy,color=RAW,s=25);ax.scatter(after,yy,color=REPAIR,s=25);ax.annotate(f'{before:.1f} -> {after:.1f}',(after,yy),xytext=(0,12),textcoords='offset points',fontsize=7.5)
 ax.set(xscale='log',xlim=(100,12000),ylim=(-.5,1.8),yticks=[0,1],yticklabels=['Polluted','Overall'],xlabel='RMSE (counts)');ax.text(.02,.89,'17 introduced peaks; same cases',transform=ax.transAxes,fontsize=7.5);panel(ax,'d','Released 600/1200-step configuration')
 save(fig,'F05_joint_validation')

def F06():
 n=pd.read_csv(F/'results/final_summary.csv');fig=plt.figure(figsize=(6.8,5.2),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1,1]);ax=fig.add_subplot(gs[:,0]);yy=np.arange(len(n));ax.hlines(yy,0,n.repair_points,color='#d2dadd',lw=1);ax.scatter(n.repair_points,yy,color=REPAIR,s=22,label='All repaired');ax.scatter(n.numerical_fallback_points,yy,color=FALL,s=16,marker='D',label='Fallback subset')
 for y,r in zip(yy,n.itertuples()):ax.annotate(str(r.repair_points),(r.repair_points,y),xytext=(5,0),textcoords='offset points',va='center',fontsize=7)
 ax.set(yticks=yy,yticklabels=n.heat_id,xlabel='Number of observations',xlim=(-10,550));ax.invert_yaxis();ax.legend(loc='lower right');panel(ax,'a','Actual repair distribution')
 z,r=case('6B03205',29,466);w=z['wavelength'];t=z['time'];y=z['raw'];D=z['reconstructed'];i,j=29,466;ix=slice(j-7,j+8);ax=fig.add_subplot(gs[0,1]);ax.plot(w[ix],y[i,ix],color=RAW,label='Raw');ax.plot(w[ix],D[i,ix],color=REPAIR,ls='--',label='Repaired');ax.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');ax.legend(fontsize=7);panel(ax,'b','Local change: 14.5 s')
 ax=fig.add_subplot(gs[1,1]);a,b=i-3,i+4;c,e=j-7,j+8;delta=D[a:b,c:e]-y[a:b,c:e];v=max(abs(delta).max(),1);im=ax.pcolormesh(w[c:e],t[a:b],delta,cmap='RdBu_r',norm=colors.TwoSlopeNorm(0,vmin=-v,vmax=v),shading='nearest',rasterized=True);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)');fig.colorbar(im,ax=ax,label='Change (counts)',shrink=.8);panel(ax,'c','Signed repair footprint')
 save(fig,'F06_actual_repair_results')

def surface_indices(z,h):
 D=z['reconstructed'];n,m=D.shape;pi,pj=np.unravel_index(np.argmax(D),D.shape);ii=np.unique(np.r_[np.linspace(0,n-1,75).astype(int),pi]);jj=np.unique(np.r_[np.linspace(0,m-1,120).astype(int),pj]);review=pd.read_csv(F/'results/four_reviewed_actual_repairs.csv');r=review[review.heat_id==h]
 if len(r):ii=np.unique(np.r_[ii,r.row]);jj=np.unique(np.r_[jj,r.column])
 if h=='6C02058':ii=np.unique(np.r_[ii,679]);jj=np.unique(np.r_[jj,936])
 if h=='6C02108':ii=np.unique(np.r_[ii,608]);jj=np.unique(np.r_[jj,1211])
 oldi=np.linspace(0,n-1,60).astype(int);oldj=np.linspace(0,m-1,90).astype(int)
 SAMPLING.append(dict(heat_id=h,native_max=float(D.max()),old_sample_max=float(D[np.ix_(oldi,oldj)].max()),new_sample_max=float(D[np.ix_(ii,jj)].max()),rows=len(ii),columns=len(jj),all_reviewed_points_retained=bool(all((x.row in ii)&(x.column in jj) for x in r.itertuples())),native_arrays_unchanged=True))
 return ii.astype(int),jj.astype(int)

def draw_surface(ax,z,h,small=False):
 ii,jj=surface_indices(z,h);w=z['wavelength'];t=z['time'];D=z['reconstructed'];W,T=np.meshgrid(w[jj],t[ii]);ax.plot_surface(W,T,D[np.ix_(ii,jj)],cmap='viridis',norm=colors.Normalize(0,65535),rasterized=True,linewidth=0,rcount=len(ii),ccount=len(jj),antialiased=False)
 ax.set(xlim=(w[0],w[-1]),ylim=(t[0],t[-1]),zlim=(0,65535),xticks=[400,700,1000],yticks=[0,t[-1]],zticks=[0,65535]);ax.set_zticklabels(['0','65k']);ax.set_xlabel('Wavelength (nm)' if not small else r'$\lambda$ (nm)',labelpad=1,fontsize=8 if not small else 7.5);ax.set_ylabel('Time (s)' if not small else 't (s)',labelpad=1,fontsize=8 if not small else 7.5);ax.tick_params(labelsize=7 if small else 8,pad=0);ax.view_init(26,-125);ax.set_title(h,pad=2,fontweight='bold',fontsize=8.5 if small else 10);ax.set_box_aspect((1.4,1,.65));ax.xaxis.pane.fill=False;ax.yaxis.pane.fill=False;ax.zaxis.pane.fill=False
 if not small:ax.set_zlabel('Intensity (counts)',labelpad=4)


def F07():
 fig=plt.figure(figsize=(6.8,9.3));fig.subplots_adjust(left=.01,right=.97,top=.98,bottom=.08,hspace=.15,wspace=.02)
 for k,h in enumerate(HEATS):ax=fig.add_subplot(5,2,k+1,projection='3d');draw_surface(ax,load(h),h,True)
 cax=fig.add_axes([.23,.035,.54,.012]);fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0,65535),cmap='viridis'),cax=cax,orientation='horizontal',label='Recorded intensity (counts)',ticks=[0,16384,32768,49152,65535]);save(fig,'F07_ten_heat_surfaces')
 for k,h in enumerate(HEATS):
  fig=plt.figure(figsize=(6.8,4.5));ax=fig.add_subplot(111,projection='3d');draw_surface(ax,load(h),h);save(fig,f'A{k+1:02d}_{h}_surface','appendix')
 pd.DataFrame(SAMPLING).drop_duplicates('heat_id').to_csv(O/'surface_sampling_audit.csv',index=False)


def appendices():
 fig,axs=plt.subplots(1,2,figsize=(6.8,2.9),layout='constrained')
 for h,c in zip(['6C02057','6C02058','6C02107','6C02108'],[RAW,REPAIR,ANOM,FALL]):
  a=json.loads((F/f'results/{h}_solver.json').read_text())['history'];x=[r['iteration'] for r in a];axs[0].plot(x,np.array([r['objective'] for r in a])/a[0]['objective'],color=c,label=h);axs[1].semilogy(x,[r['relative_change'] for r in a],color=c)
 for ax in axs:
  for x in [600,1200]:ax.axvline(x,color='#ccd1d4',ls=':',lw=.6)
  ax.set_xlabel('Iteration')
 axs[0].set_ylabel('Objective / first recorded objective');axs[0].legend(fontsize=7);axs[1].set_ylabel('Relative iterate change');axs[1].axhline(1e-4,color='gray',ls='--',lw=.8);panel(axs[0],'a','Normalized objectives');panel(axs[1],'b','Iterate changes');save(fig,'B01_solver_stability','appendix')
 z,r=case('6C02058',679,936);w=z['wavelength'];t=z['time'];y=z['raw'];D=z['reconstructed'];i,j=679,936;fig,axs=plt.subplots(2,2,figsize=(6.8,5),layout='constrained');ix=slice(j-12,j+13);axs[0,0].plot(w[ix],y[i,ix],color=RAW,label='Raw');axs[0,0].plot(w[ix],D[i,ix],color=FALL,ls='--',label='Released fallback');axs[0,0].scatter(w[j],r.latent,color=REPAIR,s=24,label='Latent estimate');axs[0,0].legend(fontsize=7);axs[0,0].set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');panel(axs[0,0],'a','Wide-peak notch: 339.5 s')
 ix=slice(i-4,i+5);axs[0,1].plot(t[ix],y[ix,j],'.-',color=RAW,label='Raw');axs[0,1].plot(t[ix],D[ix,j],'.--',color=FALL,label='Released');axs[0,1].scatter(t[i],r.wavelength_reference,color=REPAIR,label='Wave reference');axs[0,1].scatter(t[i],r.time_reference,color=ANOM,marker='x',label='Time reference');axs[0,1].set(xlabel='Time (s)',ylabel='Intensity (counts)');axs[0,1].legend(fontsize=6.8);panel(axs[0,1],'b','Disagreeing references')
 a,b=i-3,i+4;c,e=j-8,j+9;delta=D[a:b,c:e]-y[a:b,c:e];v=max(abs(delta).max(),1);im=axs[1,0].pcolormesh(w[c:e],t[a:b],delta,shading='nearest',cmap='RdBu_r',norm=colors.TwoSlopeNorm(0,vmin=-v,vmax=v),rasterized=True);fig.colorbar(im,ax=axs[1,0],label='Change (counts)');axs[1,0].set(xlabel='Wavelength (nm)',ylabel='Time (s)');panel(axs[1,0],'c','Signed output difference')
 ax=axs[1,1];ax.axis('off');panel(ax,'d','Audited quantities');ax.text(0,.92,f'Raw: {r.raw:.0f} counts\nLatent: {r.latent:.2f} counts\nReleased: {r.repaired:.2f} counts\nzt={r.time_z:.2f}, zw={r.wave_z:.4f}\nzL={r.lowrank_z:.2f}, S={r.score:.2f}\nNumerical fallback: yes\nTop code at target: no',va='top',fontsize=8.5,linespacing=1.5);save(fig,'B02_wide_peak_tradeoff','appendix')
 for h in ['6B03205','6C02058']:
  for mode in ['raw','centered']:
   with np.load(S/f'data/{h}_svd_{mode}.npz') as z:d=dict(z)
   fig,axs=plt.subplots(3,2,figsize=(6.8,6.8),layout='constrained')
   for k in range(4):axs[0,0].plot(d['wavelength'],d['basis'][k],label=f'v{k+1}');axs[0,1].plot(d['time'],d['coefficients'][:,k],label=f'sigma{k+1} u{k+1}')
   axs[0,0].set(xlabel='Wavelength (nm)',ylabel='Unit-norm loading');axs[0,1].set(xlabel='Time (s)',ylabel='Coefficient (counts)');axs[0,0].legend(fontsize=6.5,ncol=2);axs[0,1].legend(fontsize=6.5,ncol=2);panel(axs[0,0],'a',f'{h}: {mode} SVD');panel(axs[0,1],'b','Temporal coefficients')
   for ax,key,label in [(axs[1,0],'rank1','Rank-1 term'),(axs[1,1],'background','Rank-4 sum'),(axs[2,0],'residual','Residual')]:
    arr=d[key];signed=(key=='residual') or mode=='centered';lim=max(abs(arr).max(),1);im=ax.pcolormesh(d['wavelength'][::4],d['time'],arr[:,::4],shading='nearest',cmap='RdBu_r' if signed else 'viridis',vmin=-lim if signed else 0,vmax=lim,rasterized=True);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)');ax.set_title(label,loc='left');fig.colorbar(im,ax=ax,label='counts',shrink=.7)
   ax=axs[2,1];ax.axis('off');ax.text(.02,.92,('Y - channel mean' if mode=='centered' else 'Y')+' = sum sigma_k u_k v_k^T\n\nSign convention:\nlargest absolute spectral loading > 0\n\nCentered: subtract each channel mean\nRaw and centered energies stay separate.\n\nTime sample: every 8th frame',va='top',fontsize=8.2,linespacing=1.35);save(fig,f'B03_{h}_svd_{mode}','appendix')
 fig,axs=plt.subplots(2,2,figsize=(6.8,5.1),layout='constrained');follow=pd.read_csv(S/'results/channel_followup.csv').query('focus and fraction==0.15')
 for row,group in enumerate(['6B','6C']):
  f=follow[follow.heat_id.str.startswith(group)]
  for h in f.heat_id:
   tr=pd.read_csv(S/f'results/{h}_focus_trace.csv');axs[row,0].plot(tr.time_s,tr.residual,lw=.6,label=h);ii=np.linspace(0,len(tr)-1,100).astype(int);axs[row,1].scatter(tr.brightness.iloc[ii],tr.residual.iloc[ii],s=3,alpha=.4,label=h)
  axs[row,0].set(xlabel='Time (s)',ylabel='Neighbor residual (counts)');axs[row,1].set(xlabel='Excluded-window brightness (counts)',ylabel='Neighbor residual (counts)');axs[row,0].legend(fontsize=6.2,ncol=2);panel(axs[row,0],chr(97+2*row),group+' target trajectories');panel(axs[row,1],chr(98+2*row),group+' residual versus brightness')
 save(fig,'B04_persistent_channel_evidence','appendix')


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['priority','remaining'],required=True);args=parser.parse_args()
 if args.stage=='priority':F01();F02();F03();F05()
 else:F04();F06();F07();appendices()
 old=O/'asset_manifest.csv';d=pd.DataFrame(MANIFEST)
 if old.exists():prior=pd.read_csv(old);d=pd.concat([prior[~prior.figure.isin(d.figure)],d],ignore_index=True)
 d.to_csv(old,index=False)
 (O/'render_config.json').write_text(json.dumps({'base_commit':'4d28373','png_dpi':600,'main_width_inches':6.8,'figure_labels':'English','captions':'Chinese','group_colors':GROUP,'raw':RAW,'repair':REPAIR,'anomaly':ANOM,'fallback':FALL,'frozen_arrays_modified':False,'source_protocol':'../../supplement_v1/protocol.json'},indent=2)+'\n')
if __name__=='__main__':main()
