#!/usr/bin/env python
"""
Generate some simulated normalized data to test the fitting functions
STATUS: in progress
"""
import os
import numpy as np
from numpy.core.defchararray import index
from numpy.ma.core import default_fill_value
import pandas as pd
import argparse
import sys
import itertools
import scipy.stats as st
from joblib import Parallel, delayed
import lmfit
import logging
from fittinglibs import (plotting, fitting, fileio, seqfun, distribution, objfunctions, initfits, processing)

parser = argparse.ArgumentParser()
parser.add_argument('--out_params', type=str, help='Output file name')
parser.add_argument('--out_series', type=str, help='Output file name')
parser.add_argument('--xdata', type=str, help='xdata file')

def simulate_series(xdata, params_mean, params_std, n_cluster=100, func='melt_curve'):
    function = getattr(objfunctions, func)
    param_names = function(return_param_names=True)
    assert set(param_names) == set(params_mean.keys())

    simulated_params = pd.DataFrame(index=np.arange(n_cluster), columns=param_names)
    simulated_series = pd.DataFrame(index=np.arange(n_cluster), columns=['.2f'%x for x in xdata])
    for param in param_names:
        simulated_params[param] = np.random.normal(params_mean[param], params_std[param], size=(n_cluster, 1))

    simulated_series = simulated_params.apply(lambda p: function(p, xdata))

    return simulated_params, simulated_series

if __name__ == "__main__":
    args = parser.parse_args()

    xdata = pd.read_csv(args.xdata, header=None)
    params_mean = {'fmin': 0, 'dH': -37, 'Tm': 40, 'fmax':1}
    params_std = {'fmin': 0.1, 'dH': 20, 'Tm': 20, 'fmax':0.1}
    simulated_params, simulated_series = simulate_series(xdata, params_mean, params_std, n_cluster=100, func='melt_curve')

    simulated_params.to_csv(args.out_params)
    simulated_series.to_csv(args.out_series)