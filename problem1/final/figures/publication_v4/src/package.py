"""Verify the sole current manuscript and figure set; never regenerate historical prose."""
from pathlib import Path
import re,json,hashlib,subprocess
import pandas as pd,numpy as np
from PIL import Image
O=Path(__file__).resolve().parents[1];F=O.parents[1];R=F/'reports';ROOT=F.parents[1];CHECK={}
def verify(name,value):
 assert value,name
 CHECK[name]=True

def main():
 for r in pd.read_csv(O/'source_manifest.csv').itertuples():verify('source:'+r.path,hashlib.sha256((ROOT/r.path).read_bytes()).hexdigest()==r.sha256)
 doc=R/'问题一_中文论文_V4_机制与结果图优化版.md';s=doc.read_text();old=subprocess.run(['git','show','bd7cfa8:'+str(doc.relative_to(ROOT))],capture_output=True,check=True).stdout.decode()
 verify('full_math_blocks_preserved',re.findall(r'\$\$(.*?)\$\$',s,re.S)==re.findall(r'\$\$(.*?)\$\$',old,re.S))
 verify('equations_preserved',re.findall(r'\\tag\{[^}]+\}',s)==re.findall(r'\\tag\{[^}]+\}',old));verify('seven_sections',len(re.findall(r'^## 5\.',s,re.M))==7);verify('seven_body_figures',len(re.findall(r'^!\[',s,re.M))==7)
 verify('removed_F04_cd_body_calls',not re.search(r'图5-4\([cd]\)|图5-4（[cd]）',s));verify('new_appendix_calls','附图A11' in s and '附图A12' in s)
 for p in [doc,R/'问题一_V3_图表与研究附录.md',ROOT/'README.md',F/'README.md']:
  for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
   if not link.startswith(('http','mailto:')):verify('link:'+str(p.relative_to(ROOT))+':'+link,(p.parent/link.split('#')[0]).exists())
 for p in [ROOT/'README.md',F/'README.md']:
  verify('one_current_entry:'+str(p), '问题一_中文论文_V4_机制与结果图优化版.md' in p.read_text() and '问题一_中文论文_V3_研究深化与图表重构版.md' not in p.read_text())
 # Verify the four panels directly against final matrices and actual leave-point references.
 d=pd.read_csv(O/'F03_panel_data.csv');targets={'a':('6B03205',29,466,'time'),'b':('6B03205',29,466,'wave'),'c':('6C02108',608,1211,'time'),'d':('6C02108',608,1211,'wave')}
 verify('F03_four_panels',set(d.panel)==set(targets))
 for panel,(h,i,j,direction) in targets.items():
  g=d[d.panel==panel]
  with np.load(F/f'data/{h}_interface.npz') as z:
   y=z['raw'];D=z['reconstructed'];t=z['time'];w=z['wavelength'];coord=t if direction=='time' else w;ix=np.searchsorted(coord,g.coordinate);raw=y[ix,j] if direction=='time' else y[i,ix];rep=D[ix,j] if direction=='time' else D[i,ix]
   verify('F03_raw:'+panel,np.array_equal(g.raw,raw));verify('F03_repaired:'+panel,np.allclose(g.repaired,rep,rtol=0,atol=1e-8))
   alpha=(t[i]-t[i-1])/(t[i+1]-t[i-1]);beta=(w[j]-w[j-1])/(w[j+1]-w[j-1]);reference=(1-alpha)*y[i-1,j]+alpha*y[i+1,j] if direction=='time' else (1-beta)*y[i,j-1]+beta*y[i,j+1]
   verify('F03_actual_reference:'+panel,np.allclose(g.neighborhood_reference,reference,rtol=0,atol=1e-7));verify('F03_target_class:'+panel,bool(z['statistical_repair'][i,j])==(panel in ['a','b']))
 # F02 thumbnails and visible nodes are tied directly to the unchanged real matrices.
 meta=json.loads((O/'F02_data_provenance.json').read_text());source=ROOT/meta['source'];verify('F02_real_source',hashlib.sha256(source.read_bytes()).hexdigest()==meta['source_sha256']);i,end=meta['rows'];j,stop=meta['columns'];nodes=pd.read_csv(O/'F02_plot_nodes.csv')
 with np.load(source) as z:
  shape=(end-i,stop-j);t=z['time'][i:end];w=z['wavelength'][j:stop]
  for col,key in [('raw','raw'),('repaired','reconstructed'),('repair_mask','statistical_repair'),('visible','valid')]:verify('F02_actual:'+col,np.allclose(nodes[col].to_numpy().reshape(shape),z[key][i:end,j:stop],rtol=0,atol=0))
  verify('F02_actual_time',np.array_equal(nodes.time_s.to_numpy().reshape(shape)[:,0],t));verify('F02_actual_wavelength',np.allclose(nodes.wavelength_nm.to_numpy().reshape(shape)[0],w,rtol=0,atol=1e-10))
  import importlib.util
  spec=importlib.util.spec_from_file_location('f02_interface',F.parent/'v3/src/interface.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  expected=module.Spectrum(t,w,z['reconstructed'][i:end,j:stop]).evaluate(np.linspace(t[0],t[-1],2*len(t)-1),np.linspace(w[0],w[-1],2*len(w)-1))
  verify('F02_actual_bilinear_surface',hashlib.sha256(np.ascontiguousarray(expected).tobytes()).hexdigest()==meta['surface_sha256'])
 flow=json.loads((O/'flowchart_validation.json').read_text());verify('F02_no_crossing',not flow['proper_arrow_crossings']);verify('F02_colorbar_clearance',flow['colorbar_branch_clearance_pt']>3);verify('F02_retained_path',flow['retained_observations_reach_D']);verify('F02_latent_separate',flow['conditional_branch_separate'] and not meta['conditional_estimate_enters_official_surface'])
 # All 15 interpolation points must match the existing independent experiment.
 e=pd.read_csv(F.parent/'v2/results/interpolation_independent_macro.csv',header=[0,1],index_col=[0,1]);points=pd.read_csv(O/'F04_plot_data.csv');verify('F04_15_points',len(points)==15)
 for r in points.itertuples():verify('F04_exact:'+r.model+r.scenario,abs(r.rmse-e.loc[(r.model,r.scenario),('rmse','mean')])<1e-8)
 cell=pd.read_csv(O/'A12_cell_values.csv');verify('cell_weights',np.isclose(cell.weight.sum(),1) and (cell.weight>=0).all());verify('cell_original_values',np.allclose(cell.value,pd.read_csv(F/'archive/figures/publication_v3/bilinear_cell_values.csv').value))
 # Style-only re-render: verify plotted data before/after typography, plus frozen numerical sources above.
 ty=json.loads((O/'typography_validation.json').read_text());verify('actual_Times_New_Roman','Times.TTF' in ty['font_file']);verify('all_26_restyled',len(ty['figures'])==26)
 for name,r in ty['figures'].items():
  verify('plot_data_unchanged:'+name,r['plot_data_before']==r['plot_data_after']);verify('centered_titles:'+name,all(t['x']==.5 and t['ha']=='center' for t in r['centered_titles']))
 manifest=pd.read_csv(O/'asset_manifest.csv');verify('one_complete_figure_set',manifest.figure.nunique()==26 and len(manifest)==78)
 verify('flowchart_text_fits',json.loads((O/'flowchart_validation.json').read_text())['all_text_inside_boxes'])
 for r in manifest.itertuples():
  p=O/r.path;verify('asset:'+r.path,hashlib.sha256(p.read_bytes()).hexdigest()==r.sha256)
  if r.format=='pdf':
   fonts=subprocess.run(['pdffonts',str(p)],capture_output=True,text=True,check=True).stdout;verify('Times_embedded:'+r.figure,'TimesNewRoman' in fonts and 'DejaVu' not in fonts)
  if r.format=='svg':verify('English:'+r.figure,not re.search(r'[\u3400-\u9fff]',p.read_text()))
  if r.format=='png':
   im=Image.open(p);verify('600dpi:'+r.figure,abs(im.info['dpi'][0]-600)<1)
 verify('seven_current_previews',len(list((O/'review').glob('*.png')))==7)
 verify('no_stale_previews',all(p.name.endswith('_caption_page.png') for p in (O/'review').glob('*.png')))
 review=json.loads((O/'visual_review.json').read_text());verify('visual_review_completed',review['completed'])
 for path,digest in review['reviewed_sha256'].items():verify('review_current:'+path,hashlib.sha256((O/path).read_bytes()).hexdigest()==digest)
 for rel in ['problem1/final/src/detectors.py','problem1/final/config/protocol.json','problem1/v3/src/model3.py','problem1/v3/src/stable_completion.py','problem1/v3/src/interface.py','problem1/v2/src/convex_censor.py']:
  base=subprocess.run(['git','show','bd7cfa8:'+rel],capture_output=True,check=True).stdout;verify('frozen_model:'+rel,(ROOT/rel).read_bytes()==base)
 captions=re.findall(r'^图5-[1-7]　[^\n]+',s,re.M);verify('seven_captions',len(captions)==7)
 appendix=(R/'问题一_V3_图表与研究附录.md').read_text();extra=re.findall(r'^附图A1[12]　[^\n]+',appendix,re.M)
 (O/'figure_captions.md').write_text('# 当前正式图注\n\n'+'\n\n'.join(captions+extra)+'\n')
 (O/'validation.json').write_text(json.dumps({'passed':len(CHECK),'failed':0,'model_or_data_changed':False,'one_current_manuscript':str(doc.relative_to(F)),'active_figure_set':str(O.relative_to(F)),'changed_figures':['F02_model_mechanism'],'change_scope':'Horizontal scientific architecture from frozen data','appendix_added':[],'formal_width_inches':6.8,'rendering_review_completed':True},indent=2)+'\n')
 print('CURRENT_PACKAGE_VERIFIED',len(CHECK))
if __name__=='__main__':main()
