#!/usr/bin/env python
""" Fit all single clusters with minimal constraints.

Fits all single clusters.

Input:
CPsignal file
xvalues file

Output:
normalized binding series file
fit results

Sarah Denny

"""

import os
import numpy as np
import pandas as pd
import argparse
import sys
import itertools
import scipy.stats as st
from joblib import Parallel, delayed
import lmfit
import logging
from fittinglibs import (plotting, fitting, fileio, seqfun, distribution, objfunctions, initfits, processing)

### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description='fit single clusters to binding curve')
processing.add_common_args(parser.add_argument_group('common arguments'), required_x=True)

group = parser.add_argument_group('optional arguments for single cluster fitting')
group.add_argument('--subset',action="store_true", default=False,
                    help='if flagged, will only do a subset of the data for test purposes')
group.add_argument('--subset_num', default=5000, type=int,
                    help='do at most this many single clusters when the subset flag is true. default=5000')

group = parser.add_argument_group('arguments about fitting function')
group.add_argument('--func', default = 'binding_curve',
                   help='fitting function. default is "binding_curve", referring to module names in fittinglibs.objfunctions.')
group.add_argument('--params_name', nargs='+', help='name of param(s) to edit.')
group.add_argument('--params_init', nargs='+', type=float, help='new initial val(s) of param(s) to edit.')
group.add_argument('--params_vary', nargs='+', type=int, help='whether to vary val(s) of param(s) to edit.')
group.add_argument('--params_lb', nargs='+', type=float, help='new lowerbound val(s) of param(s) to edit.')
group.add_argument('--params_ub', nargs='+', type=float, help='new upperbound val(s) of param(s) to edit.')

#group.add_argument('--ft_only',action="store_true", default=False,
#                    help='if flagged, do not fit, but save the fit parameters')


def splitAndFit(fitParams, bindingSeries, numCores, index=None):
    """ Given a table of binding curves, parallelize fit. """
    if index is None:
        index = bindingSeries.index    
    logging.info('Fitting binding curves:')
    fits = (Parallel(n_jobs=numCores, verbose=10)
            (delayed(fitting.perCluster)(fitParams, bindingSeries.loc[idx])
             for idx in index))
    
    return pd.concat({idx:val for idx, val in zip(index, fits)}).unstack()
# define functions

    
def checkFitResults(fitResults):
    # did any of the stde work?
    param_names = ['fmax', 'dG', 'fmin']
    numClusters = fitResults.dropna(subset=param_names).shape[0]
    logging.info('%4.2f%% clusters have rsq>50%%'
           %(100*(fitResults.rsq > 0.5).sum()/float(numClusters)))
    logging.info('%4.2f%% clusters have stde in dG < 1'
           %(100*(fitResults.dG_stde < 1).sum()/float(numClusters)))
    logging.info('%4.2f%% clusters have stde in fmax < fmax'
           %(100*(fitResults.fmax_stde < fitResults.fmax).sum()/float(numClusters)))
    logging.info('%4.2f%% clusters have stde != 0 for at least one fit parameters'
           %(100 -100*(fitResults.loc[:, ['%s_stde'%param for param in param_names]]==0).all(axis=1).sum()/float(numClusters)))


if __name__=="__main__":    
    args = parser.parse_args()
     
    # define out file
    if args.out_file is None:
        basename = fileio.stripExtension(args.binding_series)
        if args.subset:
            args.out_file = basename + '_subset.CPfitted.gz'
        else:
            args.out_file = basename + '.CPfitted.gz'

    # load files
    logging.info("Loading binding series...")
    xvalues = np.loadtxt(args.xvalues)
    bindingSeries = fileio.loadFile(args.binding_series)  
    min_xval_col = pd.Series(xvalues, index=bindingSeries.columns).idxmin()
    
    # Initialize the fit parameters class.
    # This includes defining the fitting function, defining the xvalues,
    # and setting a function to change the fmax initial values per cluster.
    # By default, the fmax of each fit should be the max fluorescence of that cluster (per cluster fmax)
    fitParams = initfits.FitParams(args.func, xvalues, before_fit_ops=[('fmax', 'initial', np.max)])
    
    # process input args
    args = initfits.process_new_params(args)
    for param_name, param_init, param_lb, param_ub, param_vary in zip(args.params_name, args.params_init, args.params_lb, args.params_ub, args.params_vary):
        # update the initial values, upper and lower bounds for each fit param, or add additional fit params
        if param_name:
            fitParams.update_init_params(**{param_name:{'initial':param_init, 'lowerbound':param_lb, 'upperbound':param_ub, 'vary':bool(param_vary)}})
        
    # only table with at least 4 entries (with three free params)
    index_all = bindingSeries.dropna(axis=0, thresh=4).index.tolist()
    if args.subset:
        # take a subset of indices
        if len(index_all) > args.subset_num:
            index_all = index_all[:args.subset_num]

    # fit
    fitResults = splitAndFit(fitParams, bindingSeries, args.numCores, index=index_all)
    # save
    fitResults.to_csv(args.out_file, sep='\t', compression='gzip')

    fileio.saveFile(fileio.stripExtension(args.out_file)+'.fitParameters.p', fitParams)
    
    checkFitResults(fitResults)