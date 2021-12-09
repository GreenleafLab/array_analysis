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
Yuxi Ke
"""

import os
import numpy as np
from numpy.ma.core import default_fill_value
import pandas as pd
import argparse
import sys
import itertools
import scipy.stats as st
from joblib import Parallel, delayed
import lmfit
import logging
from fittinglibs.models import MeltCurveModel
from fittinglibs import (plotting, fitting, fileio, seqfun, distribution, objfunctions, initfits, processing)
### MAIN ###

################ Parse input parameters ################
logging.basicConfig(level=logging.DEBUG)
#set up command line argument parser
parser = argparse.ArgumentParser(description='fit single clusters to melt curve')
processing.add_common_args(parser.add_argument_group('common arguments'), required_x=True)
parser.add_argument('--mapfile', help='a csv file indicating which columns are the green signal')
parser.add_argument('--parallel', action='store_true', help='if flagged, use pandarallel to fit in parallel')
group = parser.add_argument_group('optional arguments for single cluster fitting')
group.add_argument('--subset',action="store_true", default=False,
                    help='if flagged, will only do a subset of the data for test purposes')
group.add_argument('--subset_num', default=5000, type=int,
                    help='do at most this many single clusters when the subset flag is true. default=5000')

# group = parser.add_argument_group('arguments about fitting function')
# group.add_argument('--func', default = 'melt_curve',
#                    help='fitting function. default is "binding_curve", referring to module names in fittinglibs.objfunctions.')
# group.add_argument('--params_name', nargs='+', help='name of param(s) to edit.')
# group.add_argument('--params_init', nargs='+', type=float, help='new initial val(s) of param(s) to edit.')
# group.add_argument('--params_vary', nargs='+', type=int, help='whether to vary val(s) of param(s) to edit.')
# group.add_argument('--params_lb', nargs='+', type=float, help='new lowerbound val(s) of param(s) to edit.')
# group.add_argument('--params_ub', nargs='+', type=float, help='new upperbound val(s) of param(s) to edit.')

#group.add_argument('--ft_only',action="store_true", default=False,
#                    help='if flagged, do not fit, but save the fit parameters')


def splitAndFit(model, meltSeries, conditions, numCores, index=None):
    """ Given a table of binding curves, parallelize fit. """
    if index is None:
        index = meltSeries.index
    logging.info('Fitting curves:')
    fits = (Parallel(n_jobs=numCores, verbose=10)
            (delayed(fitting.fit_single_clusters_in_df)(model, meltSeries.loc[idx], conditions)
             for idx in index))
    
    return pd.concat({idx:val for idx, val in zip(index, fits)}).unstack()
# define functions

    
def checkFitResults(fitResults):
    # did any of the stde work?
    param_names = ['fmax', 'dH', 'Tm', 'fmin']
    numClusters = fitResults.dropna(subset=param_names).shape[0]
    logging.info('%4.2f%% clusters have rsq>50%%'
           %(100*(fitResults.rsq > 0.5).sum()/float(numClusters)))
    logging.info('%4.2f%% clusters have stde in dH < 1'
           %(100*(fitResults.dH_stde < 1).sum()/float(numClusters)))
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
    logging.info("Loading series...")
    series_df = fileio.loadFile(args.binding_series)  
    xvalues = np.loadtxt(args.xvalues)
    green_conditions = fileio.get_conditions_from_mapfile(args.mapfile, 'green')
    norm_conditions = [s + '_norm' for s in green_conditions]

    # Initialize the model class.
    # This includes defining the fitting function, defining the xvalues,
    # and setting a function to change the fmax initial values per cluster.
    model = MeltCurveModel(T=xvalues, f_margin=np.inf)

        
    # only table with at least 5 entries (with 4 free params)
    index_all = series_df.dropna(axis=0, thresh=5).index.tolist()
    if args.subset:
        # take a subset of indices
        if len(index_all) > args.subset_num:
            index_all = index_all[:args.subset_num]

    # fit
    logging.info("Fitting curves...")
    #fitResults = splitAndFit(model, series_df, conditions=norm_conditions, numCores=args.numCores, index=index_all)
    if args.parallel:
        pandarallel.initialize()
    fitResults = fitting.fit_single_clusters_in_df(model, series_df.loc[index_all,:], conditions=norm_conditions, parallel=args.parallel)

    # save
    fitResults.to_csv(args.out_file, sep='\t', compression='gzip')
    logging.info(f"Saved results to {args.out_file}")
    #fileio.saveFile(fileio.stripExtension(args.out_file)+'.fitParameters.p', fitParams)
    
    #checkFitResults(fitResults)