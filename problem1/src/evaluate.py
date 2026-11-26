"""Evaluate saved continuous models at arbitrary coordinates inside a heat's observed domain."""
import argparse
import numpy as np
from models import load_surface
p=argparse.ArgumentParser();p.add_argument('model');p.add_argument('--times',nargs='+',type=float,required=True);p.add_argument('--wavelengths',nargs='+',type=float,required=True);p.add_argument('--output',required=True)
if __name__=='__main__':
 a=p.parse_args();m=load_surface(a.model);np.savez_compressed(a.output,time=a.times,wavelength=a.wavelengths,intensity=m.evaluate(a.times,a.wavelengths))
