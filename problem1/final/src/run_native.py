from concurrent.futures import ProcessPoolExecutor
from core import *
from fit_native import one
import time

def run(h):
 while not (P/f'data/{h}_masks.npz').exists():time.sleep(1)
 one(h)
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=3) as ex:list(ex.map(run,HEATS))
