#!/usr/bin/env python
""" Find the fmax distribution.

Returns a class describing the fmax.

Parameters:
-----------
variant_table : per-variant DataFrame including columns fmax, dG, and fmin, pvalue

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
import json

from fittinglibs import (fitting, plotting, fileio, processing, distribution, initfits, filterfunctions)



parser = argparse.ArgumentParser(description='Find Fmax Fmin distribution from variant fits')
# processing.add_common_args(parser.add_argument_group('common arguments'),
#                            required_f=False, required_a=False, required_x=False)

parser.add_argument('-vf', help='CPvariant single cluster fit result, filtered and bootstrapped to the variant level')
parser.add_argument('-o', '--output', help='json file with parameters for Fmax and Fmin distributions')
parser.add_argument('--figdir', help='Directory for plots')

group = parser.add_argument_group('additional option arguments')

group.add_argument('--variant_filter', help='query string for good fit quaulity')                
group.add_argument('-fmaxq', '--fmax_filter', default='', help='Query string to select good variants for fmax fitting')
group.add_argument('-fminq', '--fmin_filter', default='', help='Query string to select good variants for fmin fitting')



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


def findFmaxParams(good_variants_table, fit_fmin=False, variant_n_size_cutoff=10, figdir=None):
    """
    Returns a dictionary of fit parameters, mu and sigma.{a, b}
    """
    if fit_fmin:
        var_name = 'fmin'
    else:
        var_name = 'fmax'

    mu = good_variants_table[var_name].median()
    sigma_dict = fitting.fit_sigma_n_fmax(good_variants_table, fit_fmin=fit_fmin, variant_n_size_cutoff=10, return_type='dict')

    fmax_param_dict = {'mu': mu, 'sigma':sigma_dict}

    if not figdir is None:
        plotting.plotFmaxStdVsN(fmax_param_dict, good_variants_table, var_name)
        plotting.savefig(os.path.join(figdir, '%s_stderr_vs_n.pdf'%var_name))
        plotting.plt.close()

    return fmax_param_dict


if __name__=="__main__":
    args = parser.parse_args()
    figdir = args.figdir
    # processing.update_logger(logging, args.log)
        
    filterFmaxVariant = lambda table: table.query(args.variant_filter).query(args.fmax_filter)
    filterFminVariant = lambda table: table.query(args.variant_filter).query(args.fmin_filter)

    # load variant_table
    variant_table = pd.read_csv(args.vf, sep='\t')

    # filter good variants and check
    unstable_good_variants = filterFmaxVariant(variant_table)
    stable_good_variants = filterFminVariant(variant_table)

    print(('%d out of all %d variants pass fmax cutoff')
           %(len(unstable_good_variants), len(variant_table)))
    print(('%d out of all %d variants pass fmin cutoff')
           %(len(stable_good_variants), len(variant_table)))
    if (len(unstable_good_variants) < 10) or (len(stable_good_variants) < 10):
        print('Error: need more variants passing cutoffs to fit')
        print("Only saved init file... ")
        sys.exit()

    # plot fmax fmin filtering
    plotting.plotFmaxVsdG(variant_table, args.variant_filter, args.fmax_filter, T=60.0)
    plotting.savefig(os.path.join(figdir, 'fmax_vs_dG_init.pdf'))
    plotting.plt.close()
    
    plotting.plotFmaxVsdG(variant_table, args.variant_filter, args.fmin_filter, plot_fmin=True, T=20.0)
    plotting.savefig(os.path.join(figdir, 'fmin_vs_dG_init.pdf'))
    plotting.plt.close()
 
    # fit the distribution parameters
    fmax_param = findFmaxParams(unstable_good_variants, figdir=figdir)    
    fmin_param = findFmaxParams(stable_good_variants, fit_fmin=True, figdir=figdir)    

    with open(args.output, 'w') as fh:
        json.dump({'fmax': fmax_param, 'fmin': fmin_param}, fh, indent=4)

    print("Parameter file saved to %s"%args.output)