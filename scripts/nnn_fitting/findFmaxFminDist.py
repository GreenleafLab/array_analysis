#!/usr/bin/env python
""" Find the fmax distribution.

Returns a class describing the fmax.

Parameters:
-----------
variant_table : per-variant DataFrame including columns fmax, dG, and fmin, pvalue
initial_points : per-cluster DataFrame including variant_number, fmax, dG, fmin
affinity_cutoff : maximum dG (kcal/mol) to be included in fit of fmax
use_simulated : (bool) whether to use distribution of median fmax or
    subsampled fmaxes

Returns:
--------
fmaxDistObject :
"""
import argparse
import sys
import os
import numpy as np
import scipy.stats as st
from lmfit import Parameters, minimize
import pandas as pd
import logging
import seaborn as sns
import pickle

from fittinglibs import (fitting, plotting, fileio, processing, distribution, initfits, filterfunctions)



parser = argparse.ArgumentParser(description='Find Fmax Fmin distribution from variant fits')
processing.add_common_args(parser.add_argument_group('common arguments'),
                           required_f=True, required_a=True, required_x=True)

parser.add_argument('-vf', '--variant_file', help='CPvariant single cluster fit result, filtered and bootstrapped to the variant level')
parser.add_argument('-o', '--output', help='json file with parameters for Fmax and Fmin distributions')
parser.add_argument('--figdir', help='Directory for plots')

group = parser.add_argument_group('additional option arguments')
group.add_argument('--use_simulated', type=int,
                   help='set to 0 or 1 if you want to use simulated distribution (1) or'
                   'not (0). Otherwise program will decide.')
group.add_argument('-fmaxq', '--fmax_query', default='', help='Query string to select good variants for fmax fitting')
group.add_argument('-fminq', '--fmin_query', default='', help='Query string to select good variants for fmin fitting')



def useSimulatedOrActual(variant_table, cutoff):
    # if at least 20 data points have at least 10 counts in that bin, use actual
    # data. This statistics seem reasonable for fitting
    index = variant_table.dG < cutoff
    counts, binedges = np.histogram(variant_table.loc[index].numTests,
                                    np.arange(1, variant_table.numTests.max()))
    if (counts > 10).sum() >= 20:
        use_actual = True
    else:
        use_actual = False
    return use_actual

if __name__=="__main__":
    args = parser.parse_args()
    processing.update_logger(logging, args.log)
        
    # load variant_table
    filter_function = getattr(filterfunctions, args.filterfun)
    grouped = initial_points.dropna(subset=['variant_number']).groupby('variant_number')
    grouped_sub = filter_function(initial_points).dropna(subset=['variant_number']).groupby('variant_number')
    variant_table = grouped.median()
    variant_table.loc[:, 'numTests'] = grouped.size()
    variant_table.loc[:, 'fitFraction'] = grouped_sub.size()/grouped.size()
    variant_table.loc[:, 'pvalue'] = processing.findPvalueFitFraction(variant_table.fitFraction, variant_table.numTests)

    # do things on the good binders
    good_fits = variant_table.pvalue < pvalue_cutoff
    tight_binders = variant_table.loc[good_fits&(variant_table.dG<affinity_cutoff)]

    # save variant table
    variant_table.to_csv(outFile, sep='\t', compression='gzip')   

    # find fmax distribution
    logging.info(('%d out of %d variants pass cutoff')
           %(len(tight_binders), len(variant_table)))
    if len(tight_binders) < 10:
        
        logging.error('Error: need more variants passing cutoffs to fit')
        logging.error("Only saved init file... ")
        sys.exit()

    # find good variants
    plotting.plotFmaxVsKd(variant_table.loc[good_fits], parameters.find_Kd_from_dG(affinity_cutoff))
    plotting.savefig(os.path.join(figDirectory, 'fmax_vs_kd_init.pdf'))
    
    plotting.plotFmaxVsKd(variant_table.loc[good_fits], parameters.find_Kd_from_dG(affinity_cutoff),
                         plot_fmin=True)
    plotting.savefig(os.path.join(figDirectory, 'fmin_vs_kd_init.pdf'))  

    # if use_simulated is not given, decide
    if use_simulated is None:
        use_simulated = not useSimulatedOrActual(variant_table.loc[good_fits], affinity_cutoff)
    if use_simulated:
        logging.info('Using fmaxes drawn randomly from clusters')
    else:
        logging.info('Using median fmaxes of variants')

 

    # find fmax dist object
    fmaxDist = distribution.findParams(tight_binders,
                                use_simulated=use_simulated,
                                table=initial_points)
    if fmaxDist is None:
        logging.info("could not make fmax dist file. Only saved init file... ")
        sys.exit()

    plotting.savefig(os.path.join(figDirectory, 'counts_vs_n_tight_binders.pdf'))
    plotting.plt.close()
    plotting.savefig(os.path.join(figDirectory, 'offset_fmax_vs_n.pdf'))
    plotting.plt.close()
    plotting.savefig(os.path.join(figDirectory, 'stde_fmax_vs_n.pdf'))
    plotting.plt.close()
    plotting.savefig(os.path.join(figDirectory, 'fmax_dist.all.pdf'))
    
    # save
    pickle.dump(fmaxDist, open( fileio.stripExtension(outFile)+'.fmaxdist.p', "wb" ))
        

    
    # generate example distributions
    bounds = [0, distribution.findUpperboundFromFmaxDistObject(fmaxDist)]
    numExamples = tight_binders.numTests.value_counts()
    minNumDist = 20
    while minNumDist >= 5:
        validNs = pd.Series(numExamples.loc[numExamples >= minNumDist].sort_index().index.tolist())
        if len(validNs) < 3:
            minNumDist = minNumDist/2
        else:
            break
        
    if len(validNs) == 0:
        logging.info('Error: no number of measurements has at least %s variants. Not generating example distribution plots.'%minNumDist)
    else:
        if len(validNs) <= 5:
            plotTheseNs = validNs
        else:
            plotTheseNs = validNs.quantile([0, 0.1, 0.5, 0.9, 1], interpolation='nearest')
            
        for n in plotTheseNs:
            plotting.plotAnyN(tight_binders, fmaxDist, n, bounds)
            plotting.plt.xlim(bounds)
            plotting.savefig(os.path.join(figDirectory, 'fmax_dist.n_%d.pdf'%n))

    # plot fraction fit
    plotting.plotFractionFit(variant_table, pvalue_threshold=pvalue_cutoff)
    plotting.savefig(os.path.join(figDirectory, 'fraction_passing_cutoff_in_affinity_bins.pdf'))
    plotting.plt.close()
    plotting.savefig(os.path.join(figDirectory, 'histogram_fraction_fit.pdf'))
    

    

                
    
