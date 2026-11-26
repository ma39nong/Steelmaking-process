"""Refuse stale checkpoint reuse after data, model, detector, or experiment protocol changes."""
import json,hashlib
from pathlib import Path
from audit import ROOT,P,source_files

def signature():
 manifest=__import__('pandas').read_csv(ROOT/'data_manifest.csv');expected=dict(zip(manifest.heat_id,manifest.sha256))
 for path,h in zip(source_files(),manifest.heat_id):
  if hashlib.sha256(path.read_bytes()).hexdigest()!=expected[h]:raise RuntimeError(f'{h}: source changed; rerun audit and invalidate experiment checkpoints')
 inputs=[ROOT/'data_manifest.csv',P/'config/experiment.json']+[P/f'src/{f}' for f in ['audit.py','models.py','anomalies.py','experiment.py']]
 return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}

def check_or_stamp():
 now=signature();path=P/'results/validation/run_signature.json'
 if path.exists() and json.loads(path.read_text())!=now:raise RuntimeError('Data/code/config changed: refuse stale development/test checkpoints. Archive checkpoints, remove COMPLETE/selection, then rerun.')
 path.write_text(json.dumps(now,indent=2))
