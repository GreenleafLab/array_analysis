#!/usr/bin/env python
""" Fit all single clusters with minimal constraints.

Fits all varaints and refines.

Input:
    cluster = normalizedSeries,
    variant = datadir + "fitted_single_cluster/" + config["imagingExperiment"] + ".CPvariant",
    xdata = datadir + "series_normalized/" + config["imagingExperiment"] + "_xdata.txt",
    mapfile = config["mapfile"],
    fm = datadir + "fitted_fmax_fmin/%s-fmax_fmin.json" % config["imagingExperiment"]

Output:
    fitted CPvariant file for downstream analysis
    plots

Yuxi Ke
Jan 2022
"""

import os
import json
import numpy as np
from numpy.ma.core import default_fill_value
import pandas as pd
import argparse
import sys
import itertools
import scipy.stats as st
import lmfit
from pandarallel import pandarallel
import logging
from fittinglibs.models import MeltCurveModel
from fittinglibs import (plotting, fitting, fileio, seqfun, distribution, objfunctions, initfits, processing)
### MAIN ###

################ Parse input parameters ################
logging.basicConfig(level=logging.DEBUG)
#set up command line argument parser
parser = argparse.ArgumentParser(description='refine fit variants to melt curve')
processing.add_common_args(parser.add_argument_group('common arguments'), required_x=True)
parser.add_argument('--cluster_annot', help='CPannot')
parser.add_argument('-vf', help='CPvariant single cluster fit result, filtered and bootstrapped to the variant level')
parser.add_argument('-o', '--output', help='fitted CPvariant file')
group = parser.add_argument_group('arguments for fitting')
group.add_argument('--parallel', action='store_true', help='if flagged, use pandarallel to fit in parallel')
group.add_argument('--param', nargs='+', default='dH Tm', metavar="dH Tm",
                   help='param to bootstrap, excluding fmax and fmin. Default is "dH Tm"')
group.add_argument('--fmax_fmin', help='json file with parameters for Fmax and Fmin distributions')
group.add_argument('--figdir', help='Directory for plots')
group.add_argument('--n_bootstraps', action='store', type=int)
group.add_argument('--good_clusters', help='txt file with indices of good clusters')

group = parser.add_argument_group('arguments about filtering functions')
group.add_argument('--cluster_filter', help='query string for good fit quaulity')                
group.add_argument('--variant_filter', help='query string for good fit quaulity')                

group = parser.add_argument_group('additional option arguments')
group.add_argument('-vc', '--variant_col', default='SEQID', help='Name of the column to aggregate the clusters on')
group.add_argument('--subset',action="store_true", default=False,
                    help='if flagged, will only do a subset of the data for test purposes')
group.add_argument('--subset_num', default=5000, type=int,
                    help='do at most this many single clusters when the subset flag is true. default=5000')

    
def checkFitResults(fitResults):
    param_names = ['fmax', 'dH', 'Tm', 'fmin']
    numVariants = fitResults.dropna(subset=param_names).shape[0]
    logging.info('%4.2f%% clusters have ci length in dH < 50'
           %(100*(fitResults.dH_ub_final - fitResults.dH_lb_final < 50).sum()/float(numVariants)))


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
    logging.info("Loading files...")
    cluster_table = fileio.loadFile(args.binding_series)
    good_clusters = pd.read_csv(args.good_clusters, header=None).values.flatten()
    variant_table = pd.read_csv(args.vf, sep='\t').set_index(args.variant_col)
    xvalues = np.loadtxt(args.xvalues)
    annotated_clusters = pd.read_table(args.cluster_annot)[['clusterID', args.variant_col]]
    with open(args.fmax_fmin, 'r') as fh:
        fmax_params_dict = json.load(fh)

    logging.info("Files loaded")

    # filter the tables
    cluster_table = cluster_table.loc[cluster_table.index.isin(good_clusters), :]
    variant_table = variant_table.query(args.variant_filter)

    # fit
    logging.info("Fitting curves...")

    fitResults = fitting.fit_variant_bootstrap_df(cluster_table, variant_table, xvalues, annotated_clusters, fmax_params_dict, n_samples=args.n_bootstraps)
    # save
    fitResults.to_csv(args.output, sep='\t', compression='gzip')
    logging.info(f"Saved results to {args.output}")
    
    # checkFitResults(fitResults)