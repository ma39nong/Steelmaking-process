"""Run disjoint heats in an isolated results directory, then atomically publish completed caches.
Same numerical code/configuration; no agent delegation and no model changes.
"""
import sys,shutil,os
from common import P
import full_confirmation as run
actual=P
heats=sys.argv[1:]
for h in heats:
 if (actual/f'results/{h}_full_confirmation.json').exists():continue
 stage=actual/f'data/workers/{h}';(stage/'results').mkdir(parents=True,exist_ok=True);shutil.copytree(actual/'config',stage/'config',dirs_exist_ok=True)
 run.P=stage;run.HEATS=[h];run.main();source=stage/f'results/{h}_full_confirmation.json';target=actual/f'results/{h}_full_confirmation.json';temp=target.with_suffix('.publishing');shutil.copyfile(source,temp);os.replace(temp,target)
 print('published full native confirmation',h,flush=True)
