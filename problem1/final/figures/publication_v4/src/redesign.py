"""Targeted scientific redesign; frozen data and publication_v3 are read-only."""
from pathlib import Path
import sys,json,re,hashlib,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,cm
from matplotlib.patches import FancyArrowPatch,Rectangle
from matplotlib.text import Text
from mpl_toolkits.mplot3d import proj3d
O=Path(__file__).resolve().parents[1];F=O.parents[1];S=F/'supplement_v1'
sys.path.insert(0,str(S/'src'))
import publication_figures as prior
RAW=prior.RAW;TEAL=prior.REPAIR;ORANGE=prior.ANOM;PURPLE=prior.FALL;GRAY='#697c85';PALE='#edf2f3'
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'svg.hashsalt':'problem1-publication-v4'})
MAN=[];CHECKS={}
import typography

def save(fig,name):
 typography.apply(fig,name);fig.canvas.draw()
 assert not any(re.search(r'[\u3400-\u9fff]',x.get_text()) for x in fig.findobj(match=Text))
 for ext in ['png','svg','pdf']:
  p=O/f'main/{name}.{ext}';fig.savefig(p,dpi=600,bbox_inches='tight',pad_inches=.06,metadata={'Creator':'Steelmaking-process publication v4'})
  if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
  MAN.append(dict(figure=name,format=ext,path=p.relative_to(O).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
 plt.close(fig);print('rendered',name,flush=True)
def txt(ax,x,y,s,**kw):
 if '\n' in s:kw.setdefault('va','top')
 ax.text(x,y,s,transform=ax.transAxes,**kw)
def arrow(ax,start,end,color=GRAY,style='-|>',lw=1.1,rad=0):
 ax.add_patch(FancyArrowPatch(start,end,arrowstyle=style,mutation_scale=10,lw=lw,color=color,connectionstyle=f'arc3,rad={rad}',transform=ax.transAxes))
def heading(ax,x,y,num,title):txt(ax,x,y,f'{num}  {title}',fontsize=9,fontweight='bold',va='top',color=RAW)
def mechanism():
 from simplify import mechanism as current_mechanism
 current_mechanism()


def repair_atlas():
 p=pd.read_csv(F/'results/actual_repairs.csv');s=pd.read_csv(F/'results/final_summary.csv');heats=s.heat_id.tolist();lookup={h:k for k,h in enumerate(heats)};times={}
 for h in heats:
  z=prior.load(h);times[h]=float(z['time'][-1])
 p['phase']=[r.time_s/times[r.heat_id] for r in p.itertuples()];p['atlas_y']=[lookup[r.heat_id]+.12+.76*r.phase for r in p.itertuples()]
 fig=plt.figure(figsize=(6.8,6.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[2.4,1]);ax=fig.add_subplot(gs[0,:]);norm=colors.TwoSlopeNorm(vmin=-max(abs(p.difference)),vcenter=0,vmax=max(abs(p.difference)))
 for k,h in enumerate(heats):
  ax.axhspan(k,k+1,color='#f0f3f4' if k%2==0 else 'white',zorder=0);ax.axhline(k,color='#d6dee1',lw=.4);ax.text(1035,k+.5,str(int(s.set_index('heat_id').loc[h,'repair_points'])),va='center',fontsize=7,color=GRAY)
 for fall,mark,label in [(False,'o','No fallback'),(True,'D','Fallback subset')]:
  d=p[p.numerical_fallback==fall];ax.scatter(d.wavelength_nm,d.atlas_y,c=d.difference,cmap='RdBu_r',norm=norm,s=9 if not fall else 12,marker=mark,lw=.2,edgecolors='#6b7d87',alpha=.85,rasterized=False,label=label)
 ax.set(xlim=(339,1078),ylim=(10,0),yticks=np.arange(10)+.5,yticklabels=heats,xlabel='Wavelength (nm)',xticks=[400,600,800,1000]);ax.set_title('(a) All 1,043 repairs: wavelength, heat and within-heat time',loc='center',fontsize=8.6);ax.text(.99,1.015,'n',ha='right',transform=ax.transAxes,color=GRAY,fontsize=8);ax.legend(loc='upper left',bbox_to_anchor=(0,-.19),ncol=2,fontsize=7.2,frameon=False,borderaxespad=0)
 fig.colorbar(cm.ScalarMappable(norm=norm,cmap='RdBu_r'),ax=ax,location='right',fraction=.026,pad=.015,label='Estimated change D-Y (counts)',ticks=[-30000,0,30000])
 ax.text(.01,-.30,'Within each band: top = start, bottom = end of that heat.  Marker shape identifies fallback.',transform=ax.transAxes,fontsize=7,color=GRAY)
 a=fig.add_subplot(gs[1,0]);b=fig.add_subplot(gs[1,1]);stats=[]
 for fall,col,label in [(False,TEAL,'No fallback (n=578)'),(True,PURPLE,'Fallback (n=465)')]:
  d=p[p.numerical_fallback==fall];x=np.sort(abs(d.difference));a.step(x,np.arange(1,len(x)+1)/len(x),where='post',color=col,label=label)
  edges=np.linspace(0,1,11);hist,_=np.histogram(d.phase,edges);b.stairs(hist/len(d),edges,color=col,label=label)
  stats.append(dict(fallback=fall,n=len(d),median_abs_change=np.median(x),p90_abs_change=np.quantile(x,.9),mean_phase=d.phase.mean()))
 a.set(xscale='log',xlabel='Absolute estimated change (counts)',ylabel='Fraction of repaired nodes',ylim=(0,1.03));a.legend(fontsize=7);a.set_title('(b) Magnitude of the actual changes',loc='center',fontsize=8.5)
 b.set(xlabel='Within-heat time t / T',ylabel='Fraction within each subset',xlim=(0,1));b.set_title('(c) Where repairs occur in time',loc='center',fontsize=8.5)
 p.to_csv(O/'repair_atlas_points.csv',index=False);pd.DataFrame(stats).to_csv(O/'repair_atlas_summary.csv',index=False)
 CHECKS['atlas']={'total':len(p),'fallback':int(p.numerical_fallback.sum()),'nonfallback':int((~p.numerical_fallback).sum()),'each_actual_repair_plotted_once':bool(not p.duplicated(['heat_id','row','column']).any()),'phase_is_display_coordinate_only':True}
 save(fig,'F06_actual_repair_results')

def surfaces():
 fig=plt.figure(figsize=(6.8,9.2));checks=[]
 # Wider projected objects, lower boxes and shared external units remove repeated empty axes.
 for k,h in enumerate(prior.HEATS):
  row,col=divmod(k,2);ax=fig.add_axes([.005+.50*col,.84-.177*row,.49,.175],projection='3d');z=prior.load(h);D=z['reconstructed'];w=z['wavelength'];t=z['time'];ii,jj=prior.surface_indices(z,h);W,T=np.meshgrid(w[jj],t[ii]);Z=D[np.ix_(ii,jj)]
  surface=ax.plot_surface(W,T,Z,cmap='viridis',norm=colors.Normalize(0,65535),rasterized=True,linewidth=0,antialiased=False,rcount=len(ii),ccount=len(jj),clip_on=False)
  ax.set(xlim=(w[0],w[-1]),ylim=(0,t[-1]),zlim=(0,65535),xticks=([700,1000] if row<4 else [400,700,1000]),yticks=([t[-1]] if row<4 else [0,t[-1]]),zticks=[0,65535]);ax.set_zticklabels(['0','65k'] if col==0 else []);ax.set_yticklabels([f'{t[-1]:g}'] if row<4 else ['0',f'{t[-1]:g}']);ax.tick_params(labelsize=7,pad=-1);ax.view_init(23,-125);ax.set_proj_type('ortho');ax.set_box_aspect((1.65,1,.65),zoom=1.38)
  ax.text2D(.50,.90,h,ha='center',transform=ax.transAxes,fontsize=9,fontweight='bold',color=RAW)
  for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis._axinfo['grid'].update(color='#d5dde0',linewidth=.4)
  fig.canvas.draw();px,py,_=proj3d.proj_transform(W.ravel(),T.ravel(),Z.ravel(),ax.get_proj());xy=fig.transFigure.inverted().transform(ax.transData.transform(np.column_stack([px,py])));checks.append(dict(mesh_clip_on=surface.get_clip_on(),projected_x_min=float(xy[:,0].min()),projected_x_max=float(xy[:,0].max()),projected_y_min=float(xy[:,1].min()),projected_y_max=float(xy[:,1].max()),heat_id=h,native_max=float(D.max()),display_max=float(Z.max()),time_end_s=float(t[-1]),wavelength_start=float(w[0]),wavelength_end=float(w[-1]),rows=len(ii),columns=len(jj)))
 fig.text(.50,.042,'Wavelength (nm)    |    Time (s)    |    Recorded intensity (counts)',ha='center',fontsize=8)
 cax=fig.add_axes([.23,.016,.54,.009]);fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0,65535),cmap='viridis'),cax=cax,orientation='horizontal',ticks=[0,32768,65535]);cax.tick_params(labelsize=7,pad=1)
 pd.DataFrame(checks).to_csv(O/'surface_sampling_audit.csv',index=False);CHECKS['all_surface_maxima_retained']=all(x['native_max']==x['display_max'] for x in checks)
 CHECKS['no_surface_axes_clipping']=all(not x['mesh_clip_on'] for x in checks)
 CHECKS['all_projected_meshes_within_figure']=all(0<=x['projected_x_min']<=x['projected_x_max']<=1 and 0<=x['projected_y_min']<=x['projected_y_max']<=1 for x in checks)
 save(fig,'F07_ten_heat_surfaces')

def main():
 prior.O=O;prior.save=save;prior.F01();prior.F05()
 repair_atlas();surfaces()
 old=pd.read_csv(O/'asset_manifest.csv');new=pd.DataFrame(MAN);pd.concat([old[~old.figure.isin(new.figure)],new],ignore_index=True).to_csv(O/'asset_manifest.csv',index=False);(O/'checks.json').write_text(json.dumps(CHECKS,indent=2)+'\n')
 from simplify import main as simplify_current
 simplify_current()
if __name__=='__main__':main()
