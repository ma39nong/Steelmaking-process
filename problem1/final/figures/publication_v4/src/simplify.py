"""Simplify the three core figures in the sole active publication set."""
from pathlib import Path
import sys,json,hashlib,re,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.text import Text
O=Path(__file__).resolve().parents[1];F=O.parents[1];P=F.parent
sys.path.insert(0,str(F/'supplement_v1/src'));import publication_figures as prior
plt.rcParams.update({'font.size':8.2,'axes.labelsize':8.2,'axes.titlesize':8.5,'legend.fontsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'svg.hashsalt':'problem1-active-simple'})
RAW=prior.RAW;REPAIR=prior.REPAIR;REF='#7b858c';MAN=[];PANELS=[]
(O/'appendix').mkdir(exist_ok=True)
import typography

def save(fig,name,section='main'):
 typography.apply(fig,name);fig.canvas.draw();assert not any(re.search(r'[\u3400-\u9fff]',t.get_text()) for t in fig.findobj(match=Text))
 for ext in ['png','svg','pdf']:
  p=O/section/f'{name}.{ext}';fig.savefig(p,dpi=600,bbox_inches='tight',pad_inches=.055,metadata={'Creator':'Steelmaking-process current scientific figures'})
  if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
  MAN.append(dict(figure=name,format=ext,path=p.relative_to(O).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
 plt.close(fig);print('rendered',name,flush=True)

def mechanism():
 fig=plt.figure(figsize=(6.8,5.4));ax=fig.add_axes([.015,.015,.97,.97]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off');box_artists=[]
 def box(x,y,w,h,label,fill='white'):
  rect=Rectangle((x,y),w,h,fc=fill,ec='#4d626e',lw=.8);ax.add_patch(rect);text=ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=8.2,linespacing=1.1,color=RAW);box_artists.append((rect,text))
 def wire(points,head=True):
  if len(points)>2:ax.plot([x for x,y in points[:-1]],[y for x,y in points[:-1]],color=REF,lw=.9)
  ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='-|>' if head else '-',mutation_scale=9,color=REF,lw=.9,shrinkA=0,shrinkB=1.5))
 box(.025,.855,.25,.075,'Raw flame spectra', '#f0f5f6');box(.375,.855,.25,.075,'Quality check\nand classification')
 wire([(.275,.8925),(.375,.8925)])
 box(.755,.855,.22,.075,'Suspected\nupper clipping');wire([(.625,.8925),(.755,.8925)])
 # Three parallel ordinary-observation evidence branches.
 for x,label in [(.025,'Temporal\nevidence'),(.255,'Wavelength\nevidence'),(.485,'Low-rank\nevidence')]:
  box(x,.68,.18,.075,label);wire([(.5,.855),(.5,.80),(x+.09,.80),(x+.09,.755)])
 box(.20,.555,.42,.075,'Adaptive evidence fusion\nand anomaly decision')
 for x in [.115,.345,.575]:wire([(x,.68),(x,.655),(.41,.655),(.41,.63)])
 for y,label in [(.445,'Anomaly mask'),(.335,'Rank-4 Huber robust\nlow-rank completion'),(.225,'Stability check\nand fallback if needed')]:box(.20,y,.42,.075,label)
 wire([(.41,.555),(.41,.52)]);wire([(.41,.445),(.41,.41)]);wire([(.41,.335),(.41,.30)])
 # Separate recorded-range retention from conditional latent estimation.
 box(.755,.70,.22,.075,'Retain recorded\nrange information');box(.755,.565,.22,.075,'Conditional latent\nintensity estimate')
 wire([(.865,.855),(.865,.775)]);wire([(.865,.83),(.99,.83),(.99,.6025),(.975,.6025)])
 box(.025,.065,.22,.085,'Observation fidelity\nand local\nnonnegative support')
 box(.28,.065,.19,.085,'Formal repaired\nmatrix', '#f0f5f6');box(.505,.065,.22,.085,'Bilinear continuous\nreconstruction');box(.76,.065,.22,.085,'Continuous\nspectral surface','#f0f5f6')
 wire([(.41,.225),(.41,.195),(.135,.195),(.135,.15)])
 wire([(.245,.1075),(.28,.1075)]);wire([(.47,.1075),(.505,.1075)]);wire([(.725,.1075),(.76,.1075)])
 wire([(.755,.7375),(.69,.7375),(.69,.172),(.375,.172),(.375,.15)])
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 for rect,text in box_artists:
  a=rect.get_window_extent(renderer);b=text.get_window_extent(renderer);assert a.x0+1<=b.x0 and b.x1<=a.x1-1 and a.y0<=b.y0 and b.y1<=a.y1,text.get_text()
 (O/'flowchart_validation.json').write_text(json.dumps({'all_text_inside_boxes':True,'box_count':len(box_artists),'font_size_pt':8.2,'labels':[text.get_text() for rect,text in box_artists]},indent=2)+'\n')
 save(fig,'F02_model_mechanism')

def legend():return [Line2D([],[],color=RAW,lw=1.5,marker='o',ms=2.5,label='Raw'),Line2D([],[],color=REPAIR,lw=1.1,ls='--',label='Repaired'),Line2D([],[],color=REF,lw=0,marker='x',ms=5,label='Neighborhood reference')]

def evidence(ax,h,i,j,direction,panel):
 z,r=prior.case(h,i,j);Y=z['raw'];D=z['reconstructed'];t=z['time'];w=z['wavelength'];ix=slice(i-4,i+5) if direction=='time' else slice(j-7,j+8)
 x=t[ix] if direction=='time' else w[ix];raw=Y[ix,j] if direction=='time' else Y[i,ix];rep=D[ix,j] if direction=='time' else D[i,ix]
 target=t[i] if direction=='time' else w[j];reference=r.time_reference if direction=='time' else r.wavelength_reference
 ax.plot(x,raw,color=RAW,lw=1.6,marker='o',ms=2.8,zorder=2);ax.plot(x,rep,color=REPAIR,lw=1.1,ls='--',zorder=3);ax.scatter(target,reference,marker='x',s=35,color=REF,zorder=5);ax.axvline(target,color='#c4cbce',ls=':',lw=.65,zorder=1)
 ax.set(xlabel='Time (s)' if direction=='time' else 'Wavelength (nm)',ylabel='Intensity (counts)');ax.set_title(f'({panel}) '+('Temporal evidence' if direction=='time' else 'Wavelength evidence'),loc='center')
 for xx,rr,dd in zip(x,raw,rep):PANELS.append(dict(panel=panel,heat_id=h,direction=direction,coordinate=xx,raw=rr,repaired=dd,target_coordinate=target,neighborhood_reference=reference))
 return z,r

def local_evidence():
 fig,axs=plt.subplots(2,2,figsize=(6.8,4.9));fig.subplots_adjust(left=.10,right=.985,bottom=.17,top=.89,hspace=.65,wspace=.30)
 evidence(axs[0,0],'6B03205',29,466,'time','a');evidence(axs[0,1],'6B03205',29,466,'wave','b')
 evidence(axs[1,0],'6C02108',608,1211,'time','c');evidence(axs[1,1],'6C02108',608,1211,'wave','d')
 fig.text(.54,.965,'Isolated anomaly | 6B03205 | 14.5 s, 512.659 nm',ha='center',fontsize=8.5,color=RAW)
 fig.text(.54,.515,'Persistent feature | 6C02108 | 304.0 s, 766.413 nm',ha='center',fontsize=8.5,color=RAW)
 fig.legend(handles=legend(),loc='lower center',bbox_to_anchor=(.54,.015),ncol=3,frameon=False)
 save(fig,'F03_local_evidence');pd.DataFrame(PANELS).to_csv(O/'F03_panel_data.csv',index=False)

def interpolation():
 d=pd.read_csv(P/'v2/results/interpolation_independent_macro.csv',header=[0,1],index_col=[0,1]);fig,axs=plt.subplots(1,2,figsize=(6.8,2.6));fig.subplots_adjust(left=.10,right=.985,top=.88,bottom=.30,wspace=.30);rows=[]
 for ax,title,xx,scenes in [(axs[0],'(a) Time downsampling',[1,1.5,2],['time_stride_2','time_stride_3','time_stride_4']),(axs[1],'(b) Wavelength downsampling',[2,3],['wave_stride_2','wave_stride_3'])]:
  for model,col,label in [('linear',REPAIR,'Bilinear'),('pchip',RAW,'PCHIP'),('cubic',prior.ANOM,'Cubic / tensor')]:
   yy=[d.loc[(model,s),('rmse','mean')] for s in scenes];ax.plot(xx,yy,'o-',color=col,ms=4,label=label)
   rows.extend(dict(model=model,scenario=s,spacing=x,rmse=y) for s,x,y in zip(scenes,xx,yy))
  ax.set_title(title,loc='center');ax.set(xlabel='Retained interval (s)' if ax==axs[0] else 'Channel stride',ylabel='Held-out RMSE (counts)',xticks=xx)
 fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,.005),ncol=3,frameon=False)
 save(fig,'F04_continuous_reconstruction');pd.DataFrame(rows).to_csv(O/'F04_plot_data.csv',index=False)

def appendices():
 # Original global locator and the three other reviewed spikes remain research evidence.
 fig=plt.figure(figsize=(6.8,4.9),layout='constrained');gs=fig.add_gridspec(2,3);ax=fig.add_subplot(gs[0,:]);z=prior.load('6B03205');i,j=29,466;ii=np.unique(np.r_[np.arange(0,len(z['time']),4),i]);jj=np.unique(np.r_[np.arange(0,len(z['wavelength']),4),j]);im=ax.pcolormesh(z['wavelength'][jj],z['time'][ii],z['raw'][np.ix_(ii,jj)],cmap='viridis',vmin=0,vmax=65535,shading='nearest',rasterized=True);ax.scatter(z['wavelength'][j],z['time'][i],s=30,facecolors='none',edgecolors=prior.ANOM);ax.set(xlabel='Wavelength (nm)',ylabel='Time (s)',title='(a) 6B03205: representative case location');fig.colorbar(im,ax=ax,label='Recorded intensity (counts)',fraction=.018,pad=.02)
 review=pd.read_csv(F/'results/four_reviewed_actual_repairs.csv').iloc[1:]
 for k,r in enumerate(review.itertuples()):
  a=fig.add_subplot(gs[1,k]);z,e=prior.case(r.heat_id,r.row,r.column);sl=slice(r.column-6,r.column+7);a.plot(z['wavelength'][sl],z['raw'][r.row,sl],color=RAW,lw=1.4);a.plot(z['wavelength'][sl],z['reconstructed'][r.row,sl],color=REPAIR,ls='--');a.scatter(z['wavelength'][r.column],e.wavelength_reference,color=REF,marker='x',s=25);a.set(xlabel='Wavelength (nm)',ylabel='Intensity (counts)');a.set_title(f'({chr(98+k)}) {r.heat_id}',loc='center')
 fig.legend(handles=legend(),loc='outside lower center',ncol=3,frameon=False)
 save(fig,'A11_case_locator_and_spikes','appendix')
 z=prior.load('6B03205');i,j=29,466;a,b=.35,.6;corn=np.array([z['reconstructed'][i,j],z['reconstructed'][i+1,j],z['reconstructed'][i,j+1],z['reconstructed'][i+1,j+1]]);weights=np.array([(1-a)*(1-b),a*(1-b),(1-a)*b,a*b])
 fig,ax=plt.subplots(figsize=(3.5,2.9));ax.add_patch(Rectangle((0,0),1,1,fill=False,ec=REF,lw=.8))
 for (x,y),v,wt in zip([(0,0),(1,0),(0,1),(1,1)],corn,weights):
  ax.scatter(x,y,s=25,color=RAW);ax.annotate(f'{v:.1f} counts\nw={wt:.2f}',(x,y),xytext=(5,5 if y==0 else -28),textcoords='offset points',fontsize=8)
 ax.scatter(a,b,color=REPAIR,s=35);ax.annotate(f'Query: {weights@corn:.2f} counts',(a,b),xytext=(4,8),textcoords='offset points',fontsize=8,color=REPAIR)
 ax.set(xlim=(-.15,1.6),ylim=(-.22,1.22),xticks=[0,1],xticklabels=[f'{z["time"][i]:.1f}',f'{z["time"][i+1]:.1f}'],yticks=[0,1],yticklabels=[f'{z["wavelength"][j]:.3f}',f'{z["wavelength"][j+1]:.3f}'],xlabel='Time (s)',ylabel='Wavelength (nm)');ax.spines[['top','right']].set_visible(False)
 save(fig,'A12_bilinear_cell','appendix');pd.DataFrame(dict(value=corn,weight=weights)).to_csv(O/'A12_cell_values.csv',index=False)

def main():
 prior.O=O;prior.save=save
 for k,h in enumerate(prior.HEATS):
  fig=plt.figure(figsize=(6.8,4.5));ax=fig.add_subplot(111,projection='3d');prior.draw_surface(ax,prior.load(h),h);save(fig,f'A{k+1:02d}_{h}_surface','appendix')
 prior.appendices()
 mechanism();local_evidence();interpolation();appendices()
 previous=pd.read_csv(O/'asset_manifest.csv');new=pd.DataFrame(MAN);pd.concat([previous[~previous.figure.isin(new.figure)],new],ignore_index=True).to_csv(O/'asset_manifest.csv',index=False)
 (O/'figure_layout.json').write_text(json.dumps({'F02':'flowchart; no case values or formulas','F03':{'rows':2,'columns':2,'panel_roles':['isolated/time','isolated/wavelength','persistent/time','persistent/wavelength']},'F04':{'rows':1,'columns':2,'shared_legend':True},'appendix':['A11_case_locator_and_spikes','A12_bilinear_cell'],'figure_width_inches':6.8,'font_size_pt':8.2,'font_family':'Times New Roman','panel_title_alignment':'center'},indent=2)+'\n')
if __name__=='__main__':
 main()
 p=O/'checks.json';d=json.loads(p.read_text());d.pop('mechanism_case',None);d['current_core_layout']=json.loads((O/'figure_layout.json').read_text());d['current_core_case_evidence']='F03_panel_data.csv';p.write_text(json.dumps(d,indent=2)+'\n')
