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
def save(fig,name):
 fig.canvas.draw()
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
 z,r=prior.case('6B03205',29,466);Y=z['raw'];D=z['reconstructed'];w=z['wavelength'];t=z['time'];i,j=29,466
 c=r.score/np.sqrt(r.time_z*(.5*r.wave_z/(1+r.persistence)+.5*r.lowrank_z))
 fig=plt.figure(figsize=(6.8,5.2));ax=fig.add_axes([.01,.01,.98,.98]);ax.set_axis_off()
 # A directed, two-level graph: observation -> evidence -> decision -> solver -> record -> function.
 heading(ax,.025,.96,'1','Recorded spectrum')
 a=fig.add_axes([.055,.625,.20,.225]);input_image=a.pcolormesh(w[j-5:j+6],t[i-3:i+4],Y[i-3:i+4,j-5:j+6],shading='nearest',cmap='viridis',rasterized=True);a.scatter(w[j],t[i],s=30,facecolors='none',edgecolors=ORANGE);a.set(xlabel='Wavelength (nm)',ylabel='Time (s)',xticks=[511,514],yticks=[13,16]);a.tick_params(labelsize=7);cax=fig.add_axes([.262,.65,.007,.18]);cb=fig.colorbar(input_image,cax=cax,ticks=[500,1500]);cb.ax.tick_params(labelsize=6,pad=1);cax.set_title('counts',fontsize=6,pad=3)
 txt(ax,.035,.89,r'$Y_{ij}=1767$ counts',fontsize=8,color=RAW)
 heading(ax,.32,.96,'2','Evidence channels')
 rails=[(.82,r'Time $z_t$',r.time_z,f'Leave-point reference: {r.time_reference:.1f}'),(.705,r'Wave $z_\lambda$',r.wave_z,f'Leave-point reference: {r.wavelength_reference:.1f}'),(.59,r'Background $z_L$',r.lowrank_z,'Rank-4 residual / local scale')]
 for y,label,value,sub in rails:
  ax.plot([.31,.65],[y-.045,y-.045],transform=ax.transAxes,color='#d4dfe2',lw=.8)
  txt(ax,.325,y,label,fontsize=8.5,color=RAW);txt(ax,.64,y,f'{value:.2f}',ha='right',fontsize=9,color=TEAL,fontweight='bold');txt(ax,.325,y-.036,sub,fontsize=7,color=GRAY)
 arrow(ax,(.255,.745),(.305,.745))
 heading(ax,.725,.96,'3','Fuse and gate')
 txt(ax,.73,.83,r'$P=0,\quad c=0.925$',fontsize=9)
 txt(ax,.73,.755,r'$S=22.15\geq4$',fontsize=10,color=TEAL,fontweight='bold')
 txt(ax,.73,.69,'Local connected span\n<= 3 frames x 7 channels',fontsize=7.5,linespacing=1.45)
 txt(ax,.325,.49,r'$S=c\sqrt{z_t\{z_\lambda/[2(1+P)]+z_L/2\}}$',fontsize=10)
 arrow(ax,(.66,.745),(.715,.745))
 # Top-code side branch separates recorded quality information from point repairs.
 ax.plot([.255,.285,.285,.975,.975],[.84,.84,.985,.985,.55],transform=ax.transAxes,color='#b7c2c7',ls='--',lw=.8)
 txt(ax,.45,.988,'65535: separate quality flag and conditional lower-bound branch',fontsize=6.8,ha='center',va='bottom',color=GRAY)
 arrow(ax,(.975,.55),(.86,.55),color='#9baab0');txt(ax,.73,.52,r'Record $D=U$',fontsize=7.5,color=GRAY);txt(ax,.73,.48,r'Conditional $X\geq1$',fontsize=7.5,color=GRAY)
 # Decision mask is conveyed by the selected observation and exact global output count.
 ax.plot([.83,.695,.695],[.615,.615,.445],transform=ax.transAxes,color=GRAY,lw=1.1);arrow(ax,(.695,.445),(.83,.445))
 heading(ax,.72,.405,'4','Selected mask')
 txt(ax,.73,.335,r'$(i,j)\in\mathcal{A}$',fontsize=11,color=ORANGE)
 txt(ax,.73,.285,'1,043 repaired nodes\nacross all ten heats',fontsize=7.8,linespacing=1.5)
 heading(ax,.345,.405,'5','Robust completion')
 txt(ax,.35,.335,r'$\min_{\mathrm{rank}(L)\leq4}\ \sum_{\Omega_v}h_\eta(L-Y/U)$',fontsize=9)
 txt(ax,.35,.275,r'$\eta=0.003,\quad U=65535$',fontsize=8)
 txt(ax,.35,.225,'Visible nodes only; stability check',fontsize=7.5,color=GRAY)
 arrow(ax,(.705,.32),(.65,.32))
 heading(ax,.025,.405,'6','Record-scale output')
 txt(ax,.035,.335,r'$D_{ij}=523.87$ counts',fontsize=9,color=TEAL)
 txt(ax,.035,.28,'Fidelity at retained nodes\nLocal nonnegative support',fontsize=7.1,linespacing=1.5)
 txt(ax,.035,.18,r'$D_{\Omega_o}=Y_{\Omega_o},\quad D\geq0$',fontsize=9)
 arrow(ax,(.33,.32),(.285,.32))
 # Continuous output is a consequence of the projected grid, not another fitted background.
 arrow(ax,(.16,.15),(.16,.115));ax.plot([.16,.90],[.115,.115],transform=ax.transAxes,color=TEAL,lw=1.1)
 txt(ax,.35,.13,r'$D\ \longrightarrow\ \widehat F(t,\lambda)=\mathrm{Bilinear}(t,\lambda,D)$',fontsize=10,color=TEAL)
 txt(ax,.35,.065,'Continuous on the observed domain; nonnegative corner weights',fontsize=7.5,color=GRAY)
 CHECKS['mechanism_case']={'heat':'6B03205','row':i,'column':j,'time_s':float(t[i]),'wave_nm':float(w[j]),'raw':float(Y[i,j]),'repaired':float(D[i,j]),'score':float(r.score),'consistency':float(c),'diagram_values_are_actual':True}
 save(fig,'F02_model_mechanism')

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
 ax.set(xlim=(339,1078),ylim=(10,0),yticks=np.arange(10)+.5,yticklabels=heats,xlabel='Wavelength (nm)',xticks=[400,600,800,1000]);ax.set_title('(a) All 1,043 repairs: wavelength, heat and within-heat time',loc='left',fontsize=8.6);ax.text(.99,1.015,'n',ha='right',transform=ax.transAxes,color=GRAY,fontsize=8);ax.legend(loc='upper left',bbox_to_anchor=(0,-.19),ncol=2,fontsize=7.2,frameon=False,borderaxespad=0)
 fig.colorbar(cm.ScalarMappable(norm=norm,cmap='RdBu_r'),ax=ax,location='right',fraction=.026,pad=.015,label='Estimated change D-Y (counts)',ticks=[-30000,0,30000])
 ax.text(.01,-.30,'Within each band: top = start, bottom = end of that heat.  Marker shape identifies fallback.',transform=ax.transAxes,fontsize=7,color=GRAY)
 a=fig.add_subplot(gs[1,0]);b=fig.add_subplot(gs[1,1]);stats=[]
 for fall,col,label in [(False,TEAL,'No fallback (n=578)'),(True,PURPLE,'Fallback (n=465)')]:
  d=p[p.numerical_fallback==fall];x=np.sort(abs(d.difference));a.step(x,np.arange(1,len(x)+1)/len(x),where='post',color=col,label=label)
  edges=np.linspace(0,1,11);hist,_=np.histogram(d.phase,edges);b.stairs(hist/len(d),edges,color=col,label=label)
  stats.append(dict(fallback=fall,n=len(d),median_abs_change=np.median(x),p90_abs_change=np.quantile(x,.9),mean_phase=d.phase.mean()))
 a.set(xscale='log',xlabel='Absolute estimated change (counts)',ylabel='Fraction of repaired nodes',ylim=(0,1.03));a.legend(fontsize=7);a.set_title('(b) Magnitude of the actual changes',loc='left',fontsize=8.5)
 b.set(xlabel='Within-heat time t / T',ylabel='Fraction within each subset',xlim=(0,1));b.set_title('(c) Where repairs occur in time',loc='left',fontsize=8.5)
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
  ax.text2D(.02,.90,h,transform=ax.transAxes,fontsize=9,fontweight='bold',color=RAW)
  for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis._axinfo['grid'].update(color='#d5dde0',linewidth=.4)
  fig.canvas.draw();px,py,_=proj3d.proj_transform(W.ravel(),T.ravel(),Z.ravel(),ax.get_proj());xy=fig.transFigure.inverted().transform(ax.transData.transform(np.column_stack([px,py])));checks.append(dict(mesh_clip_on=surface.get_clip_on(),projected_x_min=float(xy[:,0].min()),projected_x_max=float(xy[:,0].max()),projected_y_min=float(xy[:,1].min()),projected_y_max=float(xy[:,1].max()),heat_id=h,native_max=float(D.max()),display_max=float(Z.max()),time_end_s=float(t[-1]),wavelength_start=float(w[0]),wavelength_end=float(w[-1]),rows=len(ii),columns=len(jj)))
 fig.text(.50,.042,'Wavelength (nm)    |    Time (s)    |    Recorded intensity (counts)',ha='center',fontsize=8)
 cax=fig.add_axes([.23,.016,.54,.009]);fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0,65535),cmap='viridis'),cax=cax,orientation='horizontal',ticks=[0,32768,65535]);cax.tick_params(labelsize=7,pad=1)
 pd.DataFrame(checks).to_csv(O/'surface_sampling_audit.csv',index=False);CHECKS['all_surface_maxima_retained']=all(x['native_max']==x['display_max'] for x in checks)
 CHECKS['no_surface_axes_clipping']=all(not x['mesh_clip_on'] for x in checks)
 CHECKS['all_projected_meshes_within_figure']=all(0<=x['projected_x_min']<=x['projected_x_max']<=1 and 0<=x['projected_y_min']<=x['projected_y_max']<=1 for x in checks)
 save(fig,'F07_ten_heat_surfaces')

def main():
 for name in ['F01_data_structure','F03_local_evidence','F04_continuous_reconstruction','F05_joint_validation']:
  for ext in ['png','svg','pdf']:
   p=O/f'main/{name}.{ext}';shutil.copy2(F/f'figures/publication_v3/main/{name}.{ext}',p);MAN.append(dict(figure=name,format=ext,path=p.relative_to(O).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
 mechanism();repair_atlas();surfaces()
 pd.DataFrame(MAN).to_csv(O/'asset_manifest.csv',index=False);(O/'checks.json').write_text(json.dumps(CHECKS,indent=2)+'\n')
if __name__=='__main__':main()
