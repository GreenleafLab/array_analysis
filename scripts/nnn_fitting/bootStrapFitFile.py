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
parser.add_argument('-p', '--param', nargs='+', default='dH Tm', metavar="dH Tm",
                   help='param to bootstrap, excluding fmax and fmin. Default is "dH Tm"')
parser.add_argument('-a', '--annotated_clusters', metavar=".CPannot.pkl",
                   help='file with clusters annotated by variant number')

group = parser.add_argument_group('additional option arguments')
group.add_argument('-out', '--out_file', 
                   help='output filename. default is basename of CPfitted filename')
group.add_argument('-g', '--good_cluster_index',
                   help='output file for cluster indeces that are considered good, default to basename of CPfitted + good_cluster_ind.txt')
group.add_argument('--n_samples', default=1000, type=int, metavar="N",
                   help='number of times to bootstrap samples. default = 1000.')
group.add_argument('-n', '--numCores', default=20, type=int, metavar="N",
                   help='number of cores. default = 20')
 

def filterFits(table):
    """
    Filter single cluster fits.
    Defines the criterion for a success cluster
    """
    query = "rsqr > 0.5 & fmin > -1 & fmin < 2 & fmin_stderr < fmin + 1 & fmax > 0 & fmax < 3 & fmax_stderr < fmax & Tm_stderr < 10 & dH_stderr < 100"
    return table.query(query)

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

def findPerVariantInfo(annotated_results, param_name, variant_col='SEQID'):
    """ Group results by variant number and find median param_name, fmin, and fmax. """
    variant_table = processing.findVariantTable(annotated_results,
                                                  test_stats=param_name,
                                                  variant_col=variant_col,
                                                  filterFunction=filterFits)    
    return variant_table


def findBootstrappedVariantInfo(annotated_results, variant_table, param_name, n_samples=1000, variant_col='SEQID', filter_fits=False):
    """
    Group results by variant number and bootstrap param_name, fmin, and fmax.
    """
    # filter single clusters
    if filter_fits:
        annotated_results = filterFits(annotated_results)

    # group by variant number
    grouped = annotated_results.groupby(variant_col)
    groupDict = {}
    for name, group in grouped:
        groupDict[name] = group
    
    # save a couple params
    variant_table.loc[:, 'numIter'] = n_samples
    variant_table.loc[:, 'rsqr'] = grouped.median().loc[:, 'rsqr']
    
    # parallelize bootstrapping
    params = ['fmin', 'fmax'] + param_name
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
    clusterIndexFile = args.good_cluster_index
    fittedBindingFilename = args.single_cluster_fits
    annotatedClusterFile  = args.annotated_clusters

    param     = args.param
    numCores  = args.numCores
    n_samples = args.n_samples
    
    # Name of the column to group clusters on
    # Might change according to your library file format
    variant_col = 'RefSeq'

    # find out file
    if outFile is None:
        outFile = fileio.stripExtension(args.single_cluster_fits)
    if clusterIndexFile is None:
        clusterIndexFile = fileio.stripExtension(args.single_cluster_fits) + '_good_cluster_ind.txt'

    # load data
    print('Loading data..')
    
    # load single cluster fits and make sure param is in columns
    cluster_data = fileio.loadFile(fittedBindingFilename)
    if not all(p in cluster_data.columns.tolist() for p in param):
        print("Error: Param name %s not in single cluster fit file."%param)
        sys.exit()
    
    # load annotations
    try:
        annotated_clusters = pd.read_table(annotatedClusterFile)[['clusterID', variant_col]]
    except:
        print("Check if ['clusterID', %s] are in the header of the annotated cluster file %s" % (variant_col, annotatedClusterFile))

    # load annotations for bootstrapping
    annotated_results = pd.merge(left=annotated_clusters, right=cluster_data, on='clusterID').dropna(subset=annotated_clusters.columns.tolist() + param).set_index('clusterID')
    print('Clusters annotated')

    # save
    variant_table = findPerVariantInfo(annotated_results, param + ['fmax', 'fmin'], variant_col=variant_col)
    print('Found per variant info')

    variant_table = findBootstrappedVariantInfo(annotated_results, variant_table, param, n_samples=n_samples, variant_col=variant_col, filter_fits=True)
    print('Finished bootstrapping')

    variant_table.to_csv(outFile + '.CPvariant', sep='\t')
    good_ind = filterFits(annotated_results).index.tolist()
    with open(clusterIndexFile, 'w') as fh:
        fh.writelines([l+'\n' for l in good_ind])
    print('Saved to CPvariant file and good cluster index file')