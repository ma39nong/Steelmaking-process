"""Read the frozen Problem 1 data and query its existing continuous interface."""
import argparse
from pathlib import Path
import sys
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--heat', default='6C02058')
    parser.add_argument('--time', type=float, default=445.5)
    parser.add_argument('--wavelength', type=float, default=767.07)
    args = parser.parse_args()
    path = args.root / 'problem1/final/data' / f'{args.heat}_interface.npz'
    with np.load(path, allow_pickle=False) as data:
        t, w = data['time'], data['wavelength']
        print('heat:', args.heat)
        print('shape:', (len(t), len(w)))
        print('time domain (s):', (float(t[0]), float(t[-1])))
        print('wavelength domain (nm):', (float(w[0]), float(w[-1])))
        print('actual repaired cells:', int(data['statistical_repair'].sum()))
        print('suspected saturated cells:', int(data['suspected_saturation'].sum()))
    sys.path.insert(0, str(args.root / 'problem1/v3/src'))
    from interface import Spectrum
    model = Spectrum.load(path)
    value = model.evaluate([args.time], [args.wavelength])
    print('recorded-scale query (counts):', float(value[0, 0]))
    print('offline reconstruction: use raw and causal preprocessing for online prediction')


if __name__ == '__main__':
    main()
