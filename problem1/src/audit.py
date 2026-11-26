"""Read ten independent physical wavelength grids, enforce replacement hash, audit before fitting."""
from pathlib import Path
import hashlib,json,logging,sys
import numpy as np,pandas as pd
from pyxlsb import open_workbook
ROOT=Path(__file__).resolve().parents[2];P=ROOT/'problem1'
HEATS=['6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108']
NEW_HASH='8c73aaefadb9c8d8c3f83b193d9b4583dcfdfb239446c69b7eacd4120820783c'
def source_files():
 files=[]
 for h in HEATS:
  choices=list((ROOT/'Data').glob(f'*{h}*.xlsb'))
  if h==HEATS[0]:choices=[p for p in choices if hashlib.sha256(p.read_bytes()).hexdigest()==NEW_HASH]
  if len(choices)!=1:raise ValueError(f'{h}: expected exactly one valid source, got {choices}')
  files.append(choices[0])
 return files

def read_heat(path):
 with open_workbook(str(path)) as wb:
  if len(wb.sheets)!=1:raise ValueError(f'Unexpected sheets {wb.sheets}')
  sheet=wb.sheets[0]
  with wb.get_sheet(1) as sh:
   rows=sh.rows();header=[c.v for c in next(rows)]
   indices=[j for j,v in enumerate(header) if j>0 and isinstance(v,(int,float))]
   wl=np.asarray([header[j] for j in indices],float)
   if len(wl)!=2047 or not np.all(np.diff(wl)>0):raise ValueError('Invalid physical wavelength header')
   values=np.asarray([[c.v for c in row] for row in rows],float)
 seq=values[:,0];y=values[:,indices];t=(seq-seq[0])*.5
 if not np.all(np.diff(t)>0):raise ValueError('Non-increasing Sequence; resolve duplicates before fit')
 return sheet,header,seq,t,wl,y,values[:,max(indices)+1:]

def run():
 for rel in ['data/raw','data/processed','data/masks','data/reconstructed','reports','tests','results/audit','results/eda','results/anomaly_detection','results/model_comparison','results/validation','results/final','figures/raw','figures/cleaned','figures/reconstructed','figures/comparison']:
  (P/rel).mkdir(parents=True,exist_ok=True)
 logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s',handlers=[logging.StreamHandler(),logging.FileHandler(P/'results/audit/audit.log')])
 manifest=[];quality=[];wave=[];initial=[]
 for path,h in zip(source_files(),HEATS):
  digest=hashlib.sha256(path.read_bytes()).hexdigest();logging.info('Reading %s',h)
  sheet,head,seq,t,wl,y,extra=read_heat(path)
  target=P/'data/raw'/path.name
  if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()!=digest:raise RuntimeError(f'Read-only raw copy has different hash: {target}; preserve and resolve version explicitly')
  if not target.exists():
   import shutil
   shutil.copy2(path,target);target.chmod(0o444)
  np.savez_compressed(P/'data/processed'/f'{h}.npz',sequence=seq,time=t,wavelength=wl,Y=y,process_audit_only=extra)
  finite=np.isfinite(y);sat=y==65535
  from scipy.ndimage import label
  _,nblocks=label(sat)
  common=dict(heat_id=h,records=len(t),total_columns=len(head),wavelength_channels=len(wl),wavelength_start_nm=wl[0],wavelength_end_nm=wl[-1],time_start_s=t[0],time_end_s=t[-1])
  manifest.append(dict(common,source_file=str(path.relative_to(ROOT)),sheet_name=sheet,sha256=digest,data_validity='valid' if finite.all() else 'nonfinite_requires_mask',version='replacement_SHA256_verified' if h==HEATS[0] else 'provided_attachment'))
  row=dict(common,sequence_gaps=int(np.sum(np.diff(seq)!=1)),duplicate_sequences=int(len(seq)-len(np.unique(seq))),duplicate_spectra=int(len(y)-len(np.unique(y,axis=0))),nan_count=int(np.isnan(y).sum()),inf_count=int(np.isinf(y).sum()),minimum=float(np.nanmin(y)),maximum=float(np.nanmax(y)),zeros=int((y==0).sum()),negative=int((y<0).sum()),sat65535=int(sat.sum()),sat_fraction=float(sat.mean()),saturation_components=nblocks,constant_channels=int(np.sum(np.ptp(y,axis=0)==0)),constant_rows=int(np.sum(np.ptp(y,axis=1)==0)),unchanged_adjacent_cells=int(np.sum(np.diff(y,axis=0)==0)))
  quality.append(row);d=np.diff(wl)
  wave.append(dict(heat_id=h,start_nm=wl[0],end_nm=wl[-1],delta_min_nm=d.min(),delta_max_nm=d.max(),delta_mean_nm=d.mean(),uniform=bool(np.allclose(d,d.mean())),statement_endpoint_nm=1027.678,endpoint_difference_nm=wl[-1]-1027.678))
  pd.DataFrame(dict(wavelength_nm=wl,mean=np.nanmean(y,axis=0),std=np.nanstd(y,axis=0),minimum=np.nanmin(y,axis=0),maximum=np.nanmax(y,axis=0),zero_count=(y==0).sum(0),sat_count=sat.sum(0))).to_csv(P/f'results/audit/{h}_channels.csv',index=False)
  pd.DataFrame(dict(sequence=seq,time_s=t,mean=np.nanmean(y,axis=1),std=np.nanstd(y,axis=1),maximum=np.nanmax(y,axis=1),sat_count=sat.sum(1),zero_count=(y==0).sum(1))).to_csv(P/f'results/audit/{h}_frames.csv',index=False)
  pd.DataFrame(dict(sequence=seq,**{f'process_{j}_audit_only':extra[:,j] for j in range(extra.shape[1])})).to_csv(P/f'results/audit/{h}_process_audit.csv',index=False)
  initial.append(dict(heat_id=h,invalid_count=int((~finite).sum()),sat_candidate_count=int(sat.sum()),all_zero_rows=int(np.sum(np.all(y==0,axis=1))),all_zero_columns=int(np.sum(np.all(y==0,axis=0)))))
  logging.info('%s shape=%s max=%g saturation=%d',h,y.shape,row['maximum'],row['sat65535'])
 pd.DataFrame(manifest).to_csv(ROOT/'data_manifest.csv',index=False)
 pd.DataFrame(quality).to_csv(P/'results/audit/data_quality_summary.csv',index=False)
 pd.DataFrame(wave).to_csv(P/'results/audit/wavelength_summary.csv',index=False)
 pd.DataFrame(initial).to_csv(P/'results/audit/anomaly_initial_statistics.csv',index=False)
 (P/'results/audit/COMPLETE.json').write_text(json.dumps({'heats':HEATS,'manifest_sha256':hashlib.sha256((ROOT/'data_manifest.csv').read_bytes()).hexdigest()},indent=2))
if __name__=='__main__':run()
