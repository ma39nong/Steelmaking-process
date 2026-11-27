"""Horizontal scientific architecture using frozen measured and repaired spectra."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyArrowPatch
from matplotlib.colors import ListedColormap,Normalize

def draw(O,prior,save):
 F=O.parents[1];sys.path.insert(0,str(F.parent/'v3/src'));from interface import Spectrum
 h='6C02058';source=F/f'data/{h}_interface.npz';z=prior.load(h)
 # Fixed, real 24 x 40 window; it illustrates the frozen repair mask, not a synthetic example.
 i,j,n,m=656,708,24,40;t=z['time'][i:i+n];w=z['wavelength'][j:j+m]
 Y=z['raw'][i:i+n,j:j+m];D=z['reconstructed'][i:i+n,j:j+m];A=z['statistical_repair'][i:i+n,j:j+m]
 with np.load(source) as full:V=full['valid'][i:i+n,j:j+m]
 tt=np.linspace(t[0],t[-1],2*n-1);ww=np.linspace(w[0],w[-1],2*m-1);Z=Spectrum(t,w,D).evaluate(tt,ww)
 fig=plt.figure(figsize=(6.8,3.35));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off');boxes=[];edges=[]
 def box(x,y,b,h,label):
  rect=Rectangle((x,y),b,h,facecolor='white',edgecolor='black',linewidth=.55);ax.add_patch(rect)
  txt=ax.text(x+b/2,y+h/2,label,ha='center',va='center',fontsize=8.2,linespacing=1.05,color='black');boxes.append((rect,txt))
 def wire(points,source,target):
  for a,b in zip(points[:-2],points[1:-1]):ax.plot([a[0],b[0]],[a[1],b[1]],color='black',lw=.55)
  ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='-|>',mutation_scale=6.5,lw=.55,color='black',shrinkA=0,shrinkB=1))
  edges.append({'source':source,'target':target,'points':points})
 # Parallel evidence, a single fusion node, and a compact reconstruction module.
 for y,label,key in [(.695,'Temporal\n'+r'$z_t$','time'),(.535,'Wavelength\n'+r'$z_\lambda$','wavelength'),(.375,'Low-rank\n'+r'$z_L$','lowrank')]:
  box(.19,y,.115,.09,label);wire([(.145,.58),(.17,.58),(.17,y+.045),(.19,y+.045)],'Y',key)
  wire([(.305,y+.045),(.32,y+.045),(.32,.58),(.335,.58)],key,'score')
 box(.335,.505,.075,.15,'Joint score\nand local\nscreening')
 wire([(.41,.58),(.445,.58)],'score','A')
 for y,label,key in [(.695,'Rank-4 Huber\ncompletion','fit'),(.535,'Stability check\nand fallback','stable'),(.375,'Fidelity and local\nnonnegative support','fidelity')]:box(.54,y,.195,.09,label)
 wire([(.50,.58),(.52,.58),(.52,.74),(.54,.74)],'A','fit')
 wire([(.6375,.695),(.6375,.625)],'fit','stable');wire([(.6375,.535),(.6375,.465)],'stable','fidelity')
 wire([(.145,.44),(.158,.44),(.158,.27),(.525,.27),(.525,.42),(.54,.42)],'retained_Y','fidelity')
 ax.text(.338,.30,'Retained records, including top codes',ha='center',va='center',fontsize=7.3)
 wire([(.735,.42),(.755,.42),(.755,.58),(.77,.58)],'fidelity','D')
 wire([(.825,.58),(.855,.58)],'D','bilinear_surface')
 # The latent estimate terminates separately; it never feeds the recorded-scale surface.
 box(.20,.085,.17,.105,'Suspected\ntop codes');box(.415,.085,.205,.105,'Conditional lower-bound\nestimation');box(.665,.085,.195,.105,'Separate latent\nintensity estimate')
 wire([(.158,.27),(.158,.1375),(.20,.1375)],'Y','sat')
 wire([(.37,.1375),(.415,.1375)],'sat','conditional');wire([(.62,.1375),(.665,.1375)],'conditional','latent_output')
 # Scientific thumbnails use identical true coordinates and intensity mapping.
 norm=Normalize(0,65535)
 def matrix(bounds,array,title,mask=False,ylabels=False):
  a=fig.add_axes(bounds);a.pcolormesh(w,t,array,shading='nearest',cmap=ListedColormap(['white','black']) if mask else 'viridis',norm=Normalize(0,1) if mask else norm,rasterized=not mask)
  a.tick_params(labelsize=6.8,length=2,pad=1,width=.4);a.set_xticks([w[0],w[-1]],labels=[f'{w[0]:.0f}',f'{w[-1]:.0f}']);a.set_xlabel(r'$\lambda$ (nm)',fontsize=7.2,labelpad=1)
  if ylabels:a.set_yticks([t[0],t[-1]],labels=[f'{t[0]:g}',f'{t[-1]:g}']);a.set_ylabel('Time (s)',fontsize=7.2,labelpad=1)
  else:a.set_yticks([])
  for s in a.spines.values():s.set_visible(True);s.set_linewidth(.4);s.set_color('black')
  fig.text(bounds[0]+bounds[2]/2,.86,title,ha='center',va='center',fontsize=8.2)
  return a
 matrix([.025,.40,.12,.38],Y,r'Raw $Y$',ylabels=True)
 matrix([.445,.40,.055,.38],np.ma.masked_where(~V,Y),r'$Y_{\Omega_v}$')
 matrix([.77,.40,.055,.38],D,r'Repaired $D$')
 fig.text(.30,.86,'Adaptive detection',ha='center',va='center',fontsize=8.2)
 fig.text(.6375,.86,'Robust reconstruction',ha='center',va='center',fontsize=8.2)
 surface=fig.add_axes([.855,.385,.145,.40],projection='3d');W,T=np.meshgrid(ww,tt);surface.plot_surface(W,T,Z,cmap='viridis',norm=norm,rcount=len(tt),ccount=len(ww),linewidth=0,antialiased=False,rasterized=True)
 surface.view_init(25,-125);surface.set_proj_type('ortho');surface.set_box_aspect((1.25,1,.8),zoom=1.05);surface.set_xticks([]);surface.set_yticks([]);surface.set_zticks([]);surface.set_xlabel(r'$\lambda$',labelpad=-9,fontsize=7.2);surface.set_ylabel(r'$t$',labelpad=-9,fontsize=7.2)
 for a in [surface.xaxis,surface.yaxis,surface.zaxis]:a.pane.fill=False;a.line.set_linewidth(.4);a._axinfo['grid']['linewidth']=0
 fig.text(.9275,.86,r'$\widehat F(t,\lambda)$',ha='center',va='center',fontsize=8.2)
 fig.text(.9275,.335,'Bilinear\nreconstruction',ha='center',va='center',fontsize=7.6)
 # Low/high color swatch is solely the physical measurement scale, shared by Y, D and surface.
 cax=fig.add_axes([.04,.215,.09,.013]);c=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax,orientation='horizontal',ticks=[0,65535]);cax.tick_params(labelsize=6.6,length=1.5,pad=1,width=.4);c.outline.set_linewidth(.35);cax.set_xlabel('Intensity (counts)',fontsize=7.1,labelpad=1)
 cax.get_xticklabels()[0].set_ha('left');cax.get_xticklabels()[-1].set_ha('right')
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 branch_x=ax.transData.transform((.158,.18))[0]
 colorbar_text=[*cax.get_xticklabels(),cax.xaxis.label]
 assert all(text.get_window_extent(renderer).x1 < branch_x-3 for text in colorbar_text),'Colorbar text touches the saturation branch'
 colorbar_clearance_pt=min(branch_x-text.get_window_extent(renderer).x1 for text in colorbar_text)*72/fig.dpi
 for rect,txt in boxes:
  a=rect.get_window_extent(renderer);b=txt.get_window_extent(renderer)
  assert a.x0<=b.x0 and b.x1<=a.x1 and a.y0<=b.y0 and b.y1<=a.y1,txt.get_text()
 graph={};
 for e in edges:graph.setdefault(e['source'],[]).append(e['target'])
 def reaches(start,target):
  todo=[start];seen=set()
  while todo:
   q=todo.pop()
   if q==target:return True
   if q not in seen:seen.add(q);todo.extend(graph.get(q,[]))
  return False
 assert reaches('retained_Y','D') and reaches('D','bilinear_surface') and not reaches('conditional','D') and not reaches('latent_output','bilinear_surface')
 def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
 segments=[(a,b) for e in edges for a,b in zip(e['points'][:-1],e['points'][1:])]
 proper_crossings=[]
 for k,(a,b) in enumerate(segments):
  for c,d in segments[k+1:]:
   if orient(a,b,c)*orient(a,b,d)<-1e-12 and orient(c,d,a)*orient(c,d,b)<-1e-12:proper_crossings.append([a,b,c,d])
 assert not proper_crossings,proper_crossings

 def digest(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
 pd.DataFrame({'time_s':np.repeat(t,m),'wavelength_nm':np.tile(w,n),'raw':Y.ravel(),'repair_mask':A.astype(int).ravel(),'visible':V.astype(int).ravel(),'repaired':D.ravel()}).to_csv(O/'F02_plot_nodes.csv',index=False)
 provenance={'heat_id':h,'source':str(source.relative_to(F.parents[1])),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':[i,i+n],'columns':[j,j+m],'time_range_s':[float(t[0]),float(t[-1])],'wavelength_range_nm':[float(w[0]),float(w[-1])],'raw_sha256':digest(Y),'repair_mask_sha256':digest(A),'visible_sha256':digest(V),'displayed_middle_matrix':'measured Y on frozen visible set; excluded nodes shown white','repaired_sha256':digest(D),'surface_sha256':digest(Z),'surface_shape':list(Z.shape),'surface_source':'Spectrum(t,w,D).evaluate on twice refined true-coordinate grid','recorded_scale_only':True,'conditional_estimate_enters_official_surface':False,'coordinates':'original observed time and nonuniform wavelength; pcolormesh cell centers'}
 (O/'F02_data_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 (O/'flowchart_validation.json').write_text(json.dumps({'all_text_inside_boxes':True,'box_count':len(boxes),'font_size_pt':8.2,'layout':'left-to-right, real matrices and surface','all_box_fills':'white','all_box_edges':'black','middle_panel':'actual visible-node matrix Y_Omega_v','labels':[text.get_text() for rect,text in boxes],'edges':edges,'retained_observations_reach_D':True,'conditional_branch_separate':True,'colorbar_branch_clearance_pt':colorbar_clearance_pt,'proper_arrow_crossings':proper_crossings,'shared_bus_junctions_are_intentional':True},indent=2)+'\n')
 save(fig,'F02_model_mechanism')
