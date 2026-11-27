"""Independent provenance/integrity checks for supplementary publication artifacts."""
from pathlib import Path
import re,json,hashlib,subprocess
import numpy as np,pandas as pd
from PIL import Image,ImageOps,ImageDraw
S=Path(__file__).resolve().parents[1];F=S.parent;ROOT=F.parents[1];O=F/'figures/publication_v3';results={}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1024**2),b''):h.update(x)
 return h.hexdigest()
def check(name,condition):
 assert condition,name
 results[name]=True

def main():
 hashes=pd.read_csv(S/'results/input_hashes.csv')
 for r in hashes.itertuples():check('sha:'+r.path,sha(ROOT/r.path)==r.sha256)
 # All historical scientific locks still enforced. Navigation exceptions were already
 # present at baseline; compare exact git objects rather than rewriting any lock.
 history=[]
 for version in ['final','v3','v3_1']:
  for name,digest in json.loads((ROOT/f'problem1/{version}/results/history_immutable.json').read_text()).items():
   actual=sha(ROOT/name)
   if actual!=digest:
    base=subprocess.run(['git','show',f'4d28373:{name}'],capture_output=True,check=True).stdout
    check('existing_navigation_exception:'+version+name,name in ['README.md','problem1/README.md'] and hashlib.sha256(base).hexdigest()!=digest)
    history.append(dict(version=version,path=name,locked_sha256=digest,baseline_sha256=hashlib.sha256(base).hexdigest(),current_sha256=actual,preexisting_mismatch=True))
   else:check('history:'+version+name,True)
 pd.DataFrame(history).to_csv(S/'results/historical_navigation_lock_exceptions.csv',index=False)
 for p in [*F.glob('src/*.py'),*F.glob('config/*.json')]:
  base=subprocess.run(['git','show',f'4d28373:{p.relative_to(ROOT)}'],capture_output=True,check=True).stdout
  check('frozen_source:'+p.name,p.read_bytes()==base)
 rel=pd.read_csv(S/'results/repair_reliability_points.csv');native=rel[rel.statistical_repair];check('actual_partition',len(native)==1043 and native.fallback.sum()==465 and (~native.fallback).sum()==578)
 sums=[]
 for h in native.heat_id.unique():
  d=np.load(F/f'data/{h}_interface.npz');actual=native[native.heat_id==h];check('repair_ids:'+h,set(zip(actual.row,actual.column))==set(zip(*np.where(d['statistical_repair']))));check('unchanged_observations:'+h,np.array_equal(d['raw'][~d['statistical_repair']],d['reconstructed'][~d['statistical_repair']]))
  check('nonnegative:'+h,np.isfinite(d['reconstructed']).all() and (d['reconstructed']>=0).all());sums.append(d['raw'].size)
 check('observation_count',sum(sums)==27816683)
 peaks=pd.read_csv(S/'results/introduced_peaks_17_locations.csv');replay=pd.read_csv(S/'results/introduced_peak_replay_cases.csv');check('historical_peak_replay',len(peaks)==17 and (replay.expected==replay.replayed).all() and len(replay)==13)
 sc=pd.read_csv(O/'surface_sampling_audit.csv');check('all_surface_peaks_retained',len(sc)==10 and np.array_equal(sc.native_max,sc.new_sample_max) and sc.all_reviewed_points_retained.all())
 # New region pairs do not overlap each other or historical 96-frame local slots.
 idx=pd.read_csv(S/'results/new_panel_index.csv')
 for h,g in idx.groupby('heat_id'):
  t=np.load(F/f'data/{h}_interface.npz')['time'];starts=np.linspace(0,len(t)-96,6).astype(int)
  for r in g.itertuples():check(f'region_gap:{h}:{r.split}:{r.column_start}',all(r.row_start+32<=a or r.row_start>=a+96 for a in starts))
  a=g[g.split=='development'].row_start.iloc[0];b=g[g.split=='confirmation'].row_start.iloc[0];check('split_disjoint:'+h,abs(a-b)>=32)
 proof=O/'review';proof.mkdir(exist_ok=True)
 manifest=pd.read_csv(O/'asset_manifest.csv');check('24_figures_all_formats',len(manifest)==72 and manifest.figure.nunique()==24)
 render=[]
 for r in manifest.itertuples():
  p=F/r.path;check('asset:'+r.path,sha(p)==r.sha256)
  if r.format=='svg':check('english_labels:'+r.figure,not re.search(r'[\u3400-\u9fff]',p.read_text()))
  if r.format=='png':
   im=Image.open(p);check('600dpi:'+r.figure,abs(im.info['dpi'][0]-600)<1);render.append(dict(figure=r.figure,width=im.width,height=im.height,dpi=im.info['dpi'][0]))
   if r.section=='main':
    # A4-width proof with a 6.8-inch content area at 120dpi, to inspect true scale.
    im=im.convert('RGB');im.thumbnail((816,1210));page=Image.new('RGB',(992,1403),'white');page.paste(im,((992-im.width)//2,60));ImageDraw.Draw(page).text((88,24),f'{r.figure} | A4 proof, 6.8 inch maximum width',fill='black');page.save(proof/f'{r.figure}_page.png')
 pd.DataFrame(render).to_csv(O/'render_dimensions.csv',index=False)
 paper=F/'reports/问题一_中文论文_V3_研究深化与图表重构版.md';s=paper.read_text();old=(F/'reports/问题一_中文论文_V2_英文配图版.md').read_text()
 check('seven_sections',len(re.findall(r'^## 5\.[1-7] ',s,re.M))==7);check('seven_main_images',len(re.findall(r'^!\[',s,re.M))==7);check('equations_preserved',re.findall(r'\\tag\{[^}]+\}',s)==re.findall(r'\\tag\{[^}]+\}',old))
 for p in [paper,F/'reports/问题一_V3_图表与研究附录.md']:
  for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
   if not link.startswith(('http','mailto:')):check('link:'+str(p.relative_to(F))+':'+link,(p.parent/link.split('#')[0]).exists())
 # Fair best-value annotations and absence of undefined P0 precision.
 t=F/'tables_v3';check('undefined_precision_documented','—*' in (t/'T3a.md').read_text());check('pchip_1s_best','**3299.81**' in (t/'A3.md').read_text());check('squared_model_best_visible','**710.22 ± 270.12**' in (t/'T4.md').read_text());check('weighted_ppm',round(1043/sum(sums)*1e6,2)==37.50)
 for p in (t/'latex').glob('*.tex'):check('three_line:'+p.stem,all(x in p.read_text() for x in [r'\toprule',r'\midrule',r'\bottomrule']) and not re.search(r'\{[lr|]*\|[lr|]*\}',p.read_text()))
 out={'passed':len(results),'failed':0,'navigation_lock_exceptions':len(history),'frozen_science_and_data_unchanged':True,'figure_count':24,'main_figures':7,'appendix_figures':17,'actual_repairs':1043,'optional_spline':'not undertaken'}
 (S/'results/validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
