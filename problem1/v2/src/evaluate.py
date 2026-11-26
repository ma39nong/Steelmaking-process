"""Usage: .venv/bin/python problem1/v2/src/evaluate.py HEAT --time 100 101 --wave 400 589"""
import argparse
from common import P
from interpolation import load_surface
p=argparse.ArgumentParser();p.add_argument('heat');p.add_argument('--time',type=float,nargs='+',required=True);p.add_argument('--wave',type=float,nargs='+',required=True);a=p.parse_args();m=load_surface(P/f'data/{a.heat}_model.npz');print(m.evaluate(a.time,a.wave))
