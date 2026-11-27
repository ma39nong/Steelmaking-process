"""Verify archive relocations, current links and preservation of numerical artifacts."""
from pathlib import Path
import csv,hashlib,json,re,subprocess,os,sys
ROOT=Path(__file__).resolve().parents[1];F=ROOT/'problem1/final';A=F/'archive'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 mapping=json.loads((A/'relocation_map.json').read_text());rows=list(csv.DictReader((A/'relocation_manifest.csv').open()));assert len(rows)==mapping['files']
 for r in rows:
  p=ROOT/r['archive_path'];assert p.is_file(),p;assert sha(p)==r['archived_sha256'],p;assert not (ROOT/r['original_path']).exists(),r['original_path']
  if p.suffix!='.md':assert r['original_sha256']==r['archived_sha256'],p
 broken=[];links=0
 for p in [ROOT/'README.md',*(F.rglob('*.md'))]:
  if p.is_relative_to(F/'data'):continue
  for dest in re.findall(r'\]\(([^)\n]+)\)',p.read_text()):
   if dest.startswith(('http','mailto:','#','/')):continue
   links+=1
   if not (p.parent/dest.split('#')[0].strip('<>')).exists():broken.append((str(p.relative_to(ROOT)),dest))
 assert not broken,broken
 assert {p.name for p in (F/'figures').iterdir() if p.is_dir()}=={'publication_v4'}
 # Every tracked numerical result and every data/config/lock artifact retains baseline bytes.
 base=mapping['source_commit'];files=subprocess.check_output(['git','ls-tree','-r','--name-only','-z',base],cwd=ROOT).decode().rstrip('\0').split('\0');count=0
 for name in files:
  p=Path(name)
  if p.suffix in {'.json','.csv','.npz','.bib'} and any(x in p.parts for x in ['results','config','data','delivery','references']) and '/figures/' not in name:
   old=subprocess.check_output(['git','show',base+':'+name],cwd=ROOT)
   current=(ROOT/name).read_bytes()
   if p.suffix=='.csv':old=old.replace(b'\r\n',b'\n');current=current.replace(b'\r\n',b'\n')
   assert old==current,name;count+=1
 figure=json.loads((F/'figures/publication_v4/validation.json').read_text());assert figure['failed']==0
 test=subprocess.run([sys.executable,'-m','pytest','problem1/final/tests','-q'],cwd=ROOT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'},capture_output=True,text=True);assert test.returncode==0,test.stdout+test.stderr
 test_summary=test.stdout.strip().splitlines()[-1]
 out={'archived_files':len(rows),'archived_original_bytes':sum(int(x['bytes']) for x in rows),'archived_binary_and_numeric_bytes_preserved':True,'relative_links_checked':links,'broken_links':broken,'frozen_numeric_and_config_files_unchanged':count,'current_figure_directories':['publication_v4'],'figure_package_checks':figure['passed'],'tests':test_summary,'baseline_test_issue':'Before cleanup: 1 failed, 7 passed; old history check compared editable AGENTS.md against original snapshot. Original manifests retained; historical Git snapshot and current scientific artifacts now checked separately, with CSV checkout line-ending normalization.'}
 (A/'cleanup_validation.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n');print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
