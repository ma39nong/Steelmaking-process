"""English publication figures from frozen arrays and saved results; no model fitting."""
from pathlib import Path
import csv
import hashlib
import json
import re
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parents[1]
V3=P.parent/'v3'
V31=P.parent/'v3_1'
HEATS=pd.read_csv(P/'results/final_summary.csv').heat_id.tolist()
FIGURE_RECORDS=[]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.unicode_minus':True,'svg.fonttype':'none','pdf.fonttype':42,'mathtext.fontset':'dejavuserif','savefig.dpi':300})
OUT=P/'figures/english'
A=OUT/'appendix';M=OUT/'main'
A.mkdir(parents=True,exist_ok=True);M.mkdir(parents=True,exist_ok=True)
def save(fig,name,main=False):
 d=M if main else A
 fig.canvas.draw()
 from matplotlib.text import Text
 labels=[obj.get_text() for obj in fig.findobj(match=Text) if obj.get_text()]
 if any(re.search(r'[\u3400-\u9fff]',text) for text in labels):
  raise ValueError('Non-English CJK label in '+name)
 for ext in ['png','svg','pdf']:
  path=d/(name+'.'+ext)
  fig.savefig(path,dpi=300,bbox_inches='tight',metadata={'Creator':'Steelmaking-process: frozen-data English figure export'})
  if ext=='svg':
   text=path.read_text()
   if re.search(r'[\u3400-\u9fff]',text):raise ValueError('CJK text in SVG: '+name)
   path.write_text('\n'.join(x.rstrip() for x in text.splitlines())+'\n')
  FIGURE_RECORDS.append(dict(figure=name,group='main' if main else 'appendix',format=ext,path=path.relative_to(P).as_posix(),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),visible_text_items=len(labels),cjk_text_items=0))
 plt.close(fig)

def main():
 assert (P/'results/FINAL_COMPLETE.json').exists();repairs=pd.read_csv(P/'results/actual_repairs.csv');overview=plt.figure(figsize=(14,22));overview.subplots_adjust(hspace=.5,wspace=.05)
 for z,h in enumerate(HEATS):
  d=np.load(P/f'data/{h}_interface.npz');y=d['raw'];D=d['reconstructed'];t=d['time'];w=d['wavelength'];sub=repairs[repairs.heat_id==h];r=sub.iloc[int(np.argmax(abs(sub.difference.to_numpy())))];i=int(r.row);j=int(r.column);a=max(0,i-5);b=min(len(t),i+6);c=max(0,j-10);e=min(len(w),j+11)
  fig,axs=plt.subplots(2,3,figsize=(15,8));axs[0,0].plot(w,y[i],label='Raw record',color='black',lw=1);axs[0,0].plot(w,D[i],label='Point repair',ls='--');axs[0,0].set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)',title=f't={t[i]} s');axs[0,0].legend()
  axs[0,1].plot(w[c:e],y[i,c:e],label='Raw');axs[0,1].plot(w[c:e],D[i,c:e],ls='--',label='Repair estimate');axs[0,1].set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)',title='Local spectrum at a repaired point');axs[0,1].legend()
  axs[0,2].plot(t[a:b],y[a:b,j],label='Raw');axs[0,2].plot(t[a:b],D[a:b,j],ls='--',label='Repaired');axs[0,2].set(xlabel='Time (s)',ylabel='Intensity (counts)',title=f'λ={w[j]} nm');axs[0,2].legend()
  mn=y[a:b,c:e].min();mx=y[a:b,c:e].max();lim=max(abs(D[a:b,c:e]-y[a:b,c:e]).max(),1)
  for ax,img,title in zip(axs[1],[y[a:b,c:e],D[a:b,c:e],D[a:b,c:e]-y[a:b,c:e]],['Raw local spectrum','Repaired local spectrum','Repaired - raw']):
   diff=title=='Repaired - raw';im=ax.pcolormesh(w[c:e],t[a:b],img,shading='nearest',cmap='coolwarm' if diff else 'viridis',vmin=-lim if diff else mn,vmax=lim if diff else mx,rasterized=True);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)',title=title);fig.colorbar(im,ax=ax,label='counts')
  fig.suptitle(h+': Statistical point repair');fig.tight_layout();save(fig,h+'_repair_detail')
  fig,axs=plt.subplots(1,3,figsize=(15,4));ii=np.arange(0,len(t),3);jj=np.arange(0,len(w),3)
  for ax,grid,title in zip(axs[:2],[y,D],['Raw spectral matrix','Repaired spectral matrix']):im=ax.pcolormesh(w[jj],t[ii],grid[np.ix_(ii,jj)],shading='nearest',cmap='viridis',vmin=0,vmax=65535,rasterized=True);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)',title=title);fig.colorbar(im,ax=ax,label='Intensity (counts)')
  pnt=axs[2].scatter(sub.wavelength_nm,sub.time_s,c=sub.difference,cmap='coolwarm',s=8,vmin=-max(abs(sub.difference)),vmax=max(abs(sub.difference)));axs[2].set(xlabel='Wavelength (nm)',ylabel='Time (s)',xlim=(w[0],w[-1]),ylim=(t[0],t[-1]),title=f'Changed cells ({len(sub)} points)');fig.colorbar(pnt,ax=axs[2],label='Repaired - raw (counts)');fig.tight_layout();save(fig,h+'_full_spectrum_difference')
  ii=np.linspace(0,len(t)-1,60).astype(int);jj=np.linspace(0,len(w)-1,90).astype(int);W,T=np.meshgrid(w[jj],t[ii])
  for fig in [plt.figure(figsize=(9,6)),overview]:
   ax=fig.add_subplot(1,1,1,projection='3d') if fig is not overview else fig.add_subplot(5,2,z+1,projection='3d');ax.plot_surface(W,T,D[np.ix_(ii,jj)],cmap='viridis',linewidth=0,rasterized=True,vmin=0,vmax=65535);ax.set(title=h,xlabel='Wavelength (nm)',ylabel='Time (s)',zlabel='Intensity (counts)');ax.view_init(25,-125)
   if fig is not overview:save(fig,h+'_continuous_surface')
  # Raw fast coherent transitions remain; no whole-frame deletion.
  frame=pd.read_csv(V3/f'results/{h}_frame_diagnostics.csv');k=int(np.nanargmax(frame.integral_relative_first_difference.abs().to_numpy()));sel=np.arange(max(0,k-2),min(len(t),k+3));fig,axs=plt.subplots(1,2,figsize=(12,4));axs[0].plot(t,np.trapezoid(y,w,axis=1),label='Raw');axs[0].plot(t,np.trapezoid(D,w,axis=1),ls='--',label='Repaired');axs[0].set(xlabel='Time (s)',ylabel='Integrated intensity (counts nm)',title='Integrated intensity at rapid transitions');axs[0].legend()
  for k in sel:axs[1].plot(w,D[k],label=f'{t[k]}s')
  axs[1].set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)',title='Spectral shapes of neighboring frames');axs[1].legend(fontsize=8);fig.tight_layout();save(fig,h+'_rapid_transition')
 save(overview,'all_heats_continuous_surfaces',True)
 # Four specifically reviewed real points, not merely the most visually striking cases.
 review=pd.read_csv(V31/'results/high_confidence_statistical_points.csv');fig,axs=plt.subplots(2,2,figsize=(12,8))
 for ax,r in zip(axs.flat,review.itertuples()):
  d=np.load(P/f'data/{r.heat_id}_interface.npz');j=r.column;i=r.row;ix=slice(j-6,j+7);ax.plot(d['wavelength'][ix],d['raw'][i,ix],label='Raw');ax.plot(d['wavelength'][ix],d['reconstructed'][i,ix],ls='--',label='Actual repair estimate');ax.set(title=f'{r.heat_id} t={r.time_s}s',xlabel='Wavelength (nm)',ylabel='counts');ax.legend()
 fig.tight_layout();save(fig,'four_reviewed_spikes',True)
 # Unified comparison and per-band behavior.
 q=pd.read_csv(P/'results/unified_comparison.csv');q=q[q.split=='confirmation'];names=['A_P0','B_V1_MAD','C_V31','D_adaptive4'];q=q.set_index('strategy').loc[names];fig,axs=plt.subplots(1,2,figsize=(12,4));q[['precision','recall','f1']].rename(columns={'precision':'Precision','recall':'Recall','f1':'F1'}).plot.bar(ax=axs[0],color=['#277f8e','#6da881','#d99138']);q[['pollution_rmse','unpolluted_rmse']].rename(columns={'pollution_rmse':'Polluted region','unpolluted_rmse':'Unpolluted region'}).plot.bar(ax=axs[1],color=['#277f8e','#d99138']);axs[0].set(ylabel='Detection score',xlabel='Strategy');axs[1].set(ylabel='RMSE (counts)',xlabel='Strategy');
 for ax in axs:ax.set_xticklabels(['Retain all','Bidirectional MAD','Conservative','Adaptive'],rotation=15,ha='right')
 axs[0].set_ylim(0,1.05)
 fig.tight_layout();save(fig,'end_to_end_comparison',True)
 z=pd.read_csv(P/'results/confirmation_per_heat_band.csv');fig,axs=plt.subplots(1,2,figsize=(12,4));z.groupby(['band','strategy']).pollution_rmse.mean().unstack()[names].plot.bar(ax=axs[0]);z.groupby(['band','strategy']).f1.mean().unstack()[names].plot.bar(ax=axs[1]);axs[0].set(ylabel='Polluted-region RMSE (counts)',xlabel='Wavelength region');axs[1].set(ylabel='Mean case-wise F1 score',xlabel='Wavelength region');
 for ax in axs:
  ax.set_xticklabels(['Blue','Middle','Red'],rotation=0)
  ax.legend(['Retain all','Bidirectional MAD','Conservative','Adaptive'],fontsize=8)
 fig.tight_layout();save(fig,'wavelength_region_comparison',True)
 sv=pd.read_csv(V3/'results/singular_energy.csv');fig,ax=plt.subplots(figsize=(8,4))
 for h,d in sv.groupby('heat_id'):ax.plot(d['rank'],d.energy,label=h,marker='.')
 ax.set(xlabel='Rank',ylabel='Cumulative squared singular-value energy',ylim=(.8,1.001));ax.legend(ncol=2,fontsize=8);save(fig,'singular_value_energy',True)
 n=pd.read_csv(P/'results/final_summary.csv');fig,axs=plt.subplots(1,2,figsize=(12,4));axs[0].bar(n.heat_id,n.repair_points,color='#277f8e');axs[1].bar(n.heat_id,n.repair_fraction*100,color='#6da881');axs[0].set(ylabel='Number of repaired observations');axs[1].set(ylabel='Repaired observations (%)')
 for ax in axs:ax.tick_params(axis='x',rotation=60)
 fig.tight_layout();save(fig,'repair_counts',True)
 fig,ax=plt.subplots(figsize=(16,4.2));ax.axis('off');total=int(n.repair_points.sum());boxes=[(.02,'Coordinate and quality audit\nTen raw spectral matrices'),(.28,'Local scales + combined evidence\nTransient / persistent structures'),(.54,f'Rank-4 Huber completion\nFidelity / nonnegative support\n{total:,} repaired observations'),(.82,'Bilinear spectral surfaces\nRecorded / conditional outputs')]
 for x,text in boxes:ax.text(x,.5,text,transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.7',fc='#e5f0ed',ec='#277f8e'),va='center',fontsize=10)
 for x in [.245,.505,.785]:ax.annotate('',xy=(x+.03,.5),xytext=(x,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
 ax.text(.01,.12,'Paired end-to-end comparison -> Fixed development selection -> New-region confirmation\nActual point repair across ten heats -> Continuous reconstruction validation',transform=ax.transAxes);save(fig,'method_workflow',True)
 fig,axs=plt.subplots(4,2,figsize=(12,13),layout='constrained')
 for row,h in enumerate(['6C02057','6C02058','6C02107','6C02108']):
  hist=json.loads((P/f'results/{h}_solver.json').read_text())['history'];it=[x['iteration'] for x in hist]
  axs[row,0].plot(it,[x['objective'] for x in hist],color='#277f8e');axs[row,0].set(title=h,xlabel='Iteration',ylabel='Huber objective')
  axs[row,1].semilogy(it,[x['relative_change'] for x in hist],color='#277f8e');axs[row,1].axhline(1e-4,ls='--',color='gray',label='Stopping tolerance');axs[row,1].set(title=h,xlabel='Iteration',ylabel='Relative iterate change');axs[row,1].legend(fontsize=8)
 save(fig,'solver_convergence',True)
 interpolation_comparison()
 pd.DataFrame(FIGURE_RECORDS).to_csv(OUT/'figure_manifest.csv',index=False)
 (OUT/'language_check.json').write_text(json.dumps({'unique_figures':len({x['figure'] for x in FIGURE_RECORDS}),'files':len(FIGURE_RECORDS),'formats':['png','svg','pdf'],'visible_text_cjk_count':0,'svg_cjk_count':0,'data_or_model_refitted':False},indent=2)+'\n')


def interpolation_comparison():
 d=pd.read_csv(P.parent/'v2/results/interpolation_test.csv')
 g=d.groupby(['heat_id','scenario','method']).rmse.mean().groupby(['scenario','method']).mean().unstack()
 methods=['linear','pchip','cubic','tensor_bspline'];labels=['Bilinear','PCHIP','Cubic spline','Tensor B-spline']
 fig,axs=plt.subplots(1,2,figsize=(13,4.6),layout='constrained')
 for ax,scenarios,ticks,title in [(axs[0],['time_stride_2','time_stride_3','time_stride_4'],['1.0 s','1.5 s','2.0 s'],'Time downsampling'),(axs[1],['wave_stride_2','wave_stride_3'],['Every 2nd channel','Every 3rd channel'],'Wavelength downsampling')]:
  g.loc[scenarios,methods].rename(columns=dict(zip(methods,labels))).plot.bar(ax=ax,color=['#277f8e','#6da881','#d99138','#8b6daf'])
  ax.set(title=title,xlabel='Retained sampling interval',ylabel='Held-out RMSE (counts)');ax.set_xticklabels(ticks,rotation=0);ax.legend(fontsize=8)
 save(fig,'interpolation_comparison',True)


if __name__=='__main__':main()
