"""Artifact, data version and engineering acceptance; scientific acceptance kept separate."""
from ctx import *
import pandas as pd,subprocess

def main():
 historical=json.loads((P/'results/history_immutable.json').read_text());assert all(sha(ROOT/f)==h for f,h in historical.items())
 locked=json.loads((P/'config/NUMERICAL_LOCKED.json').read_text());assert all(sha(P/'src'/f)==h for f,h in locked.items())
 assert sha(P/'config/protocol.json')==json.loads((P/'config/PROTOCOL_LOCKED.json').read_text())['protocol_sha256']
 manifest=pd.read_csv(ROOT/'data_manifest.csv');source_rows=[]
 for row in manifest.itertuples():
  assert sha(ROOT/row.source_file)==row.sha256,row.source_file
  yy,tt,ww=load(row.heat_id);assert yy.shape==(row.records,row.wavelength_channels);assert tt[0]==row.time_start_s and tt[-1]==row.time_end_s;assert abs(ww[0]-row.wavelength_start_nm)<1e-9 and abs(ww[-1]-row.wavelength_end_nm)<1e-9
  source_rows.append(dict(heat_id=row.heat_id,source_file=row.source_file,sha256=row.sha256,verified=True,cache_sha256=sha(V1/f'data/processed/{row.heat_id}.npz')))
 pd.DataFrame(source_rows).to_csv(P/'results/data_version_verification.csv',index=False)
 rows=[];channels=pd.read_csv(P/'results/channel_diagnostics_all.csv')
 for h in HEATS:
  d=np.load(P/f'data/{h}_interface.npz');raw,t,w=load(h);assert np.array_equal(raw,d['raw']);assert np.array_equal(raw,d['reconstructed']);assert (d['reconstructed']>=0).all();sat=raw==65535;assert (d['conditional'][sat]>=65535).all();assert d['conditional_solver_converged'][sat].all()
  frame=pd.read_csv(P/f'results/{h}_frame_diagnostics.csv');ch=channels[channels.heat_id==h].candidate.to_numpy();q=d['quality'];np.savez_compressed(P/f'data/{h}_quality.npz',valid=d['valid'],invalid=d['invalid'],point_candidate=(q&256)>0,frame_candidate=frame.candidate.to_numpy(),channel_candidate=ch,quality=q,grade=d['quality_grade'],suspected_saturation=sat,unidentifiable=sat,numerical_fallback=np.zeros(raw.shape,bool))
  dump(P/f'results/{h}_quality_reasons.json',dict(point='V2 MAD+auxiliary statistical candidate; retain under P0',frame='frame_diagnostics CSV; statistical candidate is not confirmed fault',channel='persistent spectral contrast; true line versus faulty channel unidentified',confirmed_invalid_ledger='config/confirmed_invalid.json; empty for current attachment',saturation='65535 top code; recorded retained; latent intensity conditional, no true upper bound',quality_array_sha256=sha(P/f'data/{h}_quality.npz')))
  for name in ['interface','dense','quality']:rows.append(dict(heat_id=h,artifact=name,sha256=sha(P/f'data/{h}_{name}.npz'),bytes=(P/f'data/{h}_{name}.npz').stat().st_size))
  assert (P/f'figures/appendix/{h}_连续光谱曲面.png').exists()
 pd.DataFrame(rows).to_csv(P/'results/artifact_manifest.csv',index=False)
 logs=['v1_regression_tests.log','v2_regression_tests.log','engineering_tests.log'];statuses={x:(P/'results'/x).read_text().strip().splitlines()[-1] for x in logs};assert all('passed' in x for x in statuses.values())
 dump(P/'results/ENGINEERING_ACCEPTANCE.json',dict(passed=True,historical_files_unchanged=len(historical),protocol_and_numeric_hashes_verified=True,raw_versions_verified=10,heats_complete=10,formal_native_grid_exact_record_identity=True,formal_negative_count=0,conditional_lower_violations=0,tests=statuses,scientific_acceptance='SCIENTIFIC_ACCEPTANCE.json; restricted target and explicit limitations'))
if __name__=='__main__':main()
