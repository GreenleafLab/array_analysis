#!/usr/bin/env python
""" Bootstrap the medians of param given single cluster fits.

Parameters:
-----------

"""

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
from joblib import Parallel, delayed
from scikits.bootstrap import bootstrap
import itertools
import warnings
from fittinglibs import fileio, processing

### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description='bootstrap on off rate fits')
parser.add_argument('-cf', '--single_cluster_fits', metavar="CPfitted.pkl",
                   help='CPfitted file with single cluster fits')
parser.add_argument('-p', '--param', nargs='+', default='dH Tm', metavar="[dH | Tm]",
                   help='param to bootstrap ["dG" | "koff" | "kobs"]. Default is "dH Tm" for off rates')
parser.add_argument('-a', '--annotated_clusters', metavar=".CPannot.pkl",
                   help='file with clusters annotated by variant number')

group = parser.add_argument_group('additional option arguments')
group.add_argument('-out', '--out_file', 
                   help='output filename. default is basename of CPfitted filename')
group.add_argument('--n_samples', default=1000, type=int, metavar="N",
                   help='number of times to bootstrap samples. default = 1000.')
group.add_argument('-n', '--numCores', default=20, type=int, metavar="N",
                   help='number of cores. default = 20')
 

def filterFits(table):
    """Filter fits specific to off rates."""
    table = table.astype(float)
    index = ((table.rsq > 0.5)&(table.fmax_stde < table.fmax)&
             (table.koff_stde < table.koff)&
             (table.fmin_stde < table.fmin))
    return table.loc[index]

def bootstrapErrors(params, group, name, n_samples):
    """ Bootstrap columns in group. """
    allbounds = {}
    for param in params:
        # bootstrap medians of e
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                bounds = bootstrap.ci(group.loc[:, param].dropna(), n_samples=n_samples, statfunction=np.median)
            except IndexError:
                bounds = np.ones(2)*np.nan
        allbounds[param+'_lb'] = bounds[0]
        allbounds[param+'_ub'] = bounds[1]
    return pd.Series(allbounds, name=name)

def findPerVariantInfo(annotated_results, param_name):
    """ Group results by variant number and find median param_name, fmin, and fmax. """
    variant_table = processing.findVariantTable(annotated_results,
                                                  parameter=param_name,
                                                  filterFunction=filterFits)    
    return variant_table

def findBootstrappedVariantInfo(annotated_results, variant_table, param_name, n_samples=1000):
    """ Group results by variant number and bootstrap param_name, fmin, and fmax. """
    # group by variant number
    grouped = annotated_results.groupby('variant_number')
    groupDict = {}
    for name, group in grouped:
        groupDict[name] = group
    
    # save a couple params
    variant_table.loc[:, 'numIter'] = n_samples
    variant_table.loc[:, 'rsq'] = grouped.median().loc[:, 'rsq']
    
    # parallelize bootstrapping
    params = ['fmin', param_name, 'fmax']
    bounds = (Parallel(n_jobs=numCores, verbose=10)
            (delayed(bootstrapErrors)(params, group, name, n_samples)
             for name, group in groupDict.items()))
    bounds = pd.concat(bounds, axis=1).transpose()
    variant_table.loc[bounds.index, bounds.columns] = bounds
    variant_table.loc[:, params] = variant_table.loc[:, ['%s_init'%param for param in params]].values
    
    return variant_table


##### SCRIPT #####
if __name__ == '__main__':
    # load files
    args = parser.parse_args()
    outFile = args.out_file
    fittedBindingFilename = args.single_cluster_fits
    annotatedClusterFile  = args.annotated_clusters

    param     = args.param
    numCores  = args.numCores
    n_samples = args.n_samples
    
    # find out file
    if outFile is None:
        outFile = fileio.stripExtension(args.single_cluster_fits)

    # laod data
    print('Loading data..')
    
    # load sing cluster fits and make sure param is in columns
    cluster_data = fileio.loadFile(fittedBindingFilename)
    if not param in cluster_data.columns.tolist():
        print("Error: Param name %s not in single cluster fit file."%param)
        sys.exit()
    
    # load annotations
    try:
        annotated_clusters = fileio.loadFile(annotatedClusterFile)[['clusterID', 'RefSeq']]
    except:
        print("Check if ['clusterID', 'RefSeq'] are in the header of the annotated cluster file %s"%annotatedClusterFile)

    # load annotations for bootstrapping
    annotated_results = pd.merge(left=annotated_clusters, right=cluster_data, on='clusterID').dropna(subset=annotated_clusters.columns.tolist() + [param])
    
    # save
    variant_table = findPerVariantInfo(annotated_results, param)
    variant_table = findBootstrappedVariantInfo(annotated_results, variant_table, param, n_samples=n_samples)
    variant_table.to_csv(outFile + '.CPvariant', sep='\t')
