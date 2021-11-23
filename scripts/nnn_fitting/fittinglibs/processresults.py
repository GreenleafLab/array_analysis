#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from scikits.bootstrap import bootstrap
import warnings
import itertools
import ipdb
import os
import scipy.stats as st
import matplotlib as mpl
from fittinglibs import plotting, fitting, distribution
from plotting import fix_axes
import ipdb

def findExtInList(directory, ext):
    """Find files in directory with extension."""
    if os.path.isdir(directory):
        files = os.listdir(directory)
        return [os.path.join(directory, i) for i in files if i.find(ext)>-1 and i.find(ext)==len(i)-len(ext)]
    else:
        print 'No directory named: %s'%directory
        return []

def loadFile(directory, ext):
    """Load particular file with extension in directory."""
    filenames = findExtInList(directory, ext)
    if len(filenames)==1:
        data = fileio.loadFile(filenames[0])
        print 'Loaded file: %s'%filenames[0]
    else:
        data = None
        if len(filenames) > 1:
            print 'More than one file found: %s'%('\t\n'.join(filenames))
        else:
            print 'Could not find extension %s in directory %s'%(ext, directory)
    return data

def initialize(directory):
    """ Find and load all files of normal structure in directory.
    
    Returns loaded file if extension was found, None otherwise.
    
    Parameters:
    -----------
    directory : the directory containing the CPvariant, CPseries, CPfitted, etc files.
    
    Outputs:
    --------
    variant_table, binding_series, cluster_table, annotated_clusters, time_series, tiles, times, concentrations, timedict
    
    variant_table : per variant fits. 'normalized.CPvariant' in directory.
    binding_series : normalized binding series. 'normalized.CPseries.pkl'
    cluster_table : single cluster fits. 'CPfitted.pkl' in directory
    annotated_clusters : annotated clusters. 'CPannot.pkl' in directory
    time_series : binned time series. 'CPtimeseries.pkl' in directory.
    tiles: clusters annotated with tile info. 'CPtiles.pkl' in directory.
    times : single time vector giving corresponding time for time_series columns. 'times' in directory
    concentrations : vector of concentrations for a binding series. 'concentrations.txt' in directory/../
    timedict : dict of times with keys = tiles. 'timeDict.p' in directory.

    """
    variant_table = loadFile(directory, 'normalized.CPvariant')
    binding_series = loadFile(directory, 'normalized.CPseries.pkl')
    cluster_table = loadFile(directory, 'normalized.CPfitted.pkl')
    annotated_clusters = loadFile(directory,'CPannot.pkl')
    time_series = loadFile(directory, 'CPtimeseries.pkl')
    tiles = loadFile(directory, 'CPtiles.pkl')

    # x
    times = loadFile(directory, 'times')
    concentrations = loadFile(os.path.join(directory, '..'), 'concentrations.txt')
    timedict = loadFile(directory, 'timeDict.p')
        
    return variant_table, binding_series, cluster_table, annotated_clusters, time_series, tiles, times, concentrations, timedict

def errorPropagateAverage(sigmas, weights):
    """Propagate error using weights."""
    sigma_out = np.sqrt(np.sum([np.power(weight*sigma/weights.sum(), 2)
                                if weight != 0
                                else 0
                                for sigma, weight in itertools.izip(sigmas, weights)]))
    return sigma_out

def errorPropagateAverageAll(sigmas, weights, index=None):
    """Propagate error using weights."""
    if index is None:
        index = sigma_dGs.index
    subweights = weights.loc[index].astype(float)
    subsigmas = sigmas.loc[index].astype(float)
    
    sigma_out = pd.Series([errorPropagateAverage(subsigmas.loc[i], weights.loc[i])
                           for i in subsigmas.index], index=subsigmas.index)
    
    return sigma_out

def weightedAverage(values, weights):
    """Average using weights."""
    average = (np.sum([value*weight
                       if weight != 0
                       else 0
                       for value, weight in itertools.izip(values, weights)])/
                weights.sum())
    return average

def weightedAverageAll(values, weights, index=None):
    """Return weighted average """
    if index is None:
        index = values.index
    
    subvalues = values.loc[index]
    subweights = weights.loc[index]
    average = pd.Series([weightedAverage(subvalues.loc[i], weights.loc[i])
                           for i in subvalues.index], index=subvalues.index)
    return average 
    
def getValueInTable(series, name):
    """ Given an index name, return value of index if in series, else return nan."""
    if name in series.index.tolist():
        return series.loc[name]
    else:
        return np.nan
    
class perVariant():
    """For a given CPvariant file, plot stuff."""
    def __init__(self, variant_table=None, annotated_clusters=None, binding_series=None, x=None, cluster_table=None, tiles=None, fmaxdist=None):

        self.binding_series = binding_series
        if annotated_clusters is not None:
            self.annotated_clusters = annotated_clusters.loc[:,'variant_number']
        self.variant_table = variant_table
        self.x = x
        self.cluster_table = cluster_table
        self.tiles = tiles
        self.fmaxdist = fmaxdist

    def getVariantBindingSeries(self,variant ):
        """Return binding series for clusters of a particular variant."""
        index = self.annotated_clusters == variant
        return self.binding_series.loc[self.annotated_clusters.index.tolist()].loc[index]
    
    def getVariantTiles(self, variant):
        """Return tile numbers for clusters of a particular variant."""
        index = self.annotated_clusters == variant 
        return self.tiles.loc[self.annotated_clusters.index.tolist()].loc[index]
    
    
    def plotBindingCurve(self, variant, annotate=True, func_kwargs={}, plot_init=False):
        """Plot a binding curve of a particular variant."""
        subSeries = self.getVariantBindingSeries(variant)
        if len(subSeries)==0:
            print 'No fluorescence data associated with variant %s'%str(variant)
            return
        concentrations = self.x
        variant_table = self.variant_table
        
        # plot
        fig = plt.figure(figsize=(3,3))
        ax = fig.add_subplot(111)
        plotting.plotFitCurve(concentrations,
                            subSeries,
                            variant_table.loc[variant],
                            ax=ax, kwargs=func_kwargs)
        if annotate:
            names = ['dG', 'dG_lb', 'dG_ub', 'fmax', 'fmax_lb', 'fmax_ub', 'numTests', 'pvalue', 'rsq', 'flag']
            vec = pd.Series([getValueInTable(variant_table.loc[variant], name) for name in names], index=names)
            ax = annotateBindingCurve(vec, ax)
        
        if plot_init:
            initialPoints = distribution.findInitialPoints(variant_table)
            more_concentrations = np.logspace(np.log10(concentrations[0]/10.), np.log10(concentrations[-1]*10.))
            fit = fitting.bindingCurveObjectiveFunction(fitting.returnParamsFromResults(initialPoints.loc[variant]), more_concentrations, **func_kwargs)
            ax.plot(more_concentrations, fit, 'm--')            

    def plotOffrateCurve(self, variant, annotate=False, numtiles=None, tiles=None, no_errors=False):
        """Plot an off rate curve of a particular variant."""
        subSeries = self.getVariantBindingSeries(variant).dropna(how='all')
        timeDict = self.x
        variant_table = self.variant_table
        fitParameters = pd.DataFrame(columns=['fmax', 'koff', 'fmin'])
        
        if len(subSeries)==0:
            print 'No fluorescence data associated with variant %s'%str(variant)
            return
        
        # find tiles with variant cluster
        clusterTiles = self.getVariantTiles(variant)
        if tiles is None:
            tileNum = clusterTiles.value_counts()
            tileInit = tileNum.idxmax() # start by plotting the tile with the most clusters
            if numtiles is None: # plot as many tiles as there are clusters if num tiles is not specified
                numtiles = len(tileNum)
            otherTiles = tileNum.iloc[1:numtiles].index.tolist()
        else:
            # if tile is specified, just plot that tile.
            tileInit, otherTiles = tiles[0], tiles[1:]

        if np.in1d(clusterTiles, np.hstack([tileInit, otherTiles])).sum()==0:
            print 'No fluorescence data associated with variant %s'%str(variant)
            return 
        
        # initiae plot
        fig = plt.figure(figsize=(3.5,3))
        ax = fig.add_subplot(111)
        
        # plot fit and one tile's data
        if no_errors:
            errors = [np.nan*np.ones(len(timeDict[tileInit]))]*2
        else:
            errors = None

        plotting.plotFitCurve(timeDict[tileInit], subSeries.loc[clusterTiles==tileInit], variant_table.loc[variant], fitParameters, ax=ax, capsize=0,
                 log_axis=False, func=fitting.objectiveFunctionOffRates, fittype='off', errors=errors)
        
        # plot the other tile's data
        for tile in otherTiles:
            # plot fit and one tile's data
            if no_errors:
                errors = [np.nan*np.ones(len(timeDict[tileInit]))]*2
            else:
                errors = None
            plotting.plotDataErrorbars(timeDict[tile], subSeries.loc[clusterTiles==tile], ax, capsize=0, errors=errors)
        
        # annotate
        if annotate:
            names = ['koff', 'koff_lb', 'koff_ub', 'fmax', 'fmax_lb', 'fmax_ub', 'numTests', 'pvalue', 'rsq']
            vec = pd.Series([getValueInTable(variant_table.loc[variant], name) for name in names], index=names)
            ax = annotateOffrate(vec, ax)
        return ax
    
    def plotOffrateCurveBinned(self, variant, annotate=False, no_errors=False):
        """Plot off rate curve as above, but bin the tiles into time bins."""
        subSeries = self.getVariantBindingSeries(variant).dropna(how='all')
        timeDict = self.x
        variant_table = self.variant_table
        fitParameters = pd.DataFrame(columns=['fmax', 'koff', 'fmin'])
        clusterTiles = self.getVariantTiles(variant)
        
        if len(subSeries)==0:
            print 'No fluorescence data associated with variant %s'%str(variant)
            return
        
        time_delta = np.min([(np.array(times[1:]) - np.array(times[:-1])).min() for times in timeDict.values()])
        max_time = np.max([np.max(times) for times in timeDict.values()])
        time_bins = np.arange(0, max_time+time_delta*2, time_delta)
        
        # bin the times
        d  = []
        for key, times in timeDict.items():
            idx = np.digitize(times, time_bins) - 1
            subMat = subSeries.loc[clusterTiles==key]
            subMat.columns=idx
            d.append(subMat)
        newMat = pd.concat(d)
                
        # drop empty bins
        times = np.array([(time_bins[i]+time_bins[i+1])*0.5 for i in newMat.columns.tolist()])
        
        if no_errors:
            errors = [np.nan*np.ones(len(times))]*2
        else:
            errors = None
        # initiae plot
        fig = plt.figure(figsize=(3.5,3))
        ax = fig.add_subplot(111)
        plotting.plotFitCurve(times, newMat, variant_table.loc[variant], fitParameters, ax=ax, capsize=0,
                 log_axis=False, func=fitting.objectiveFunctionOffRates, fittype='off', errors=errors)
            
    def plotClusterOffrates(self, cluster=None, variant=None, idx=None):
        """Plot an off rate curve of a particular cluster."""
        times = self.x
        
        if cluster is not None:
            fluorescence = self.binding_series.loc[cluster]
        else:
            if variant is None:
                print 'Error: need to define either variant or cluster!'
                return
            subSeries = self.getVariantBindingSeries(variant)
            if idx is not None:
                fluorescence = subSeries.iloc[idx]
            else:
                fluorescence = subSeries.iloc[0]
        cluster = fluorescence.name
        cluster_table = self.cluster_table
        if cluster in cluster_table.index:
            results = cluster_table.loc[cluster]
        else:
            results = pd.Series(index=cluster_table.columns)
        
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        
        fitParameters = pd.DataFrame(columns=['fmax', 'koff', 'fmin'])
        plotting.plotFitCurve(times, fluorescence, results, fitParameters, ax=ax,
                 log_axis=False, func=fitting.objectiveFunctionOffRates, fittype='off')

    def plotClusterBinding(self, variant, cluster=None, idx=None):
        """Plot a binding curve of a particular cluster."""
        concentrations = self.x
        
        if cluster is not None:
            fluorescence = subSeries.loc[cluster]
        else:
            subSeries = self.getVariantBindingSeries(variant)
            if idx is not None:
                fluorescence = subSeries.iloc[idx]
            else:
                fluorescence = subSeries.iloc[0]
         
        cluster = fluorescence.name
        cluster_table = self.cluster_table
        
        # plot
        fig = plt.figure(figsize=(3,3))
        ax = fig.add_subplot(111)
        plotting.plotFitCurve(concentrations,
                            fluorescence,
                            cluster_table.loc[cluster],
                            ax=ax)
     
    def plotBootstrappedDist(self, variant, param, log_axis=False):
        """Plot the distribution of a param."""
        variant_table = self.variant_table
        cluster_table = self.cluster_table
        subSeries = self.getVariantBindingSeries(variant)
        
        params = cluster_table.loc[subSeries.index, param]
        
        # make bootstrapped dist
        if log_axis:
            vec = np.log10(params.dropna())
            med = np.log10(variant_table.loc[variant, param])
            ub = np.log10(variant_table.loc[variant, param+'_ub'])
            lb = np.log10(variant_table.loc[variant, param+'_lb'])
            xlabel = 'log '+param
        else:
            vec = params.dropna()
            med = variant_table.loc[variant, param]
            ub = variant_table.loc[variant, param+'_ub']
            lb = variant_table.loc[variant, param+'_lb']
            xlabel = param
        plt.figure(figsize=(4,3))
        sns.distplot(vec, color='r', kde=False)
        plt.axvline(med, color='0.5', linestyle='-')
        plt.axvline(lb, color='0.5', linestyle=':')
        plt.axvline(ub, color='0.5', linestyle=':')
        plt.xlabel(xlabel)
        fix_axes(plt.gca())
        plt.tight_layout()
        
    def plotFractionFit(self):
        """Plot the fraction fit."""
        variant_table = self.variant_table
        pvalue_cutoff = 0.01
        # plot
        binwidth=0.01
        bins=np.arange(0,1+binwidth, binwidth)
        plt.figure(figsize=(4, 3.5))
        plt.hist(variant_table.loc[variant_table.pvalue <= pvalue_cutoff].fitFraction.values,
                 alpha=0.5, color='red', bins=bins, label='passing cutoff')
        plt.hist(variant_table.loc[variant_table.pvalue > pvalue_cutoff].fitFraction.values,
                 alpha=0.5, color='grey', bins=bins,  label='fails cutoff')
        plt.ylabel('number of variants')
        plt.xlabel('fraction fit')
        plt.legend(loc='upper left')
        plt.tight_layout()
        fix_axes(plt.gca())
    
    def plotErrorByNumberofMeasurements(self, xlim=[1,100], ylim=[0,1], ax=None, color='r', marker='>'):
        """ Plot how the ci changes with number of measruements. """
        variant_table = self.variant_table
        variant_table.loc[:, 'ci_width'] = (variant_table.dG_ub - variant_table.dG_lb)/2
        x, y, yerr = returnFractionGroupedBy(variant_table, 'numTests', 'ci_width')
        if ax is None:
            fig = plt.figure(figsize=(4,3))
            ax = fig.add_subplot(111)
            plt.subplots_adjust(left=0.15, bottom=0.15, top=0.95, right=0.95)
        ax.scatter(x, y, c=color, edgecolor='k', linewidth=0.5, marker=marker, s=5)
        ax.errorbar(x, y, yerr, fmt='-', elinewidth=1, capsize=0, capthick=1,
                    color=color, linestyle='', ecolor='k', linewidth=0.5)
        plt.xlabel('number of measurements')
        plt.ylabel('average error (kcal/mol)')
        plt.xlim(xlim)
        plt.ylim(ylim)
        fix_axes(ax)

    def plotErrorByDeltaGBin(self, binedges=np.arange(-13, -6.1, 0.1), ylim=[0,0.8], min_n=5, xlim=None, ax=None, color='r', marker='>'):
        """ Plot how the ci changes with dG. """
        variant_table = self.variant_table

        variant_table.loc[:, 'ci_width'] = (variant_table.dG_ub - variant_table.dG_lb)/2
        variant_table.loc[:, 'dG_bin'] = np.digitize(variant_table.dG, binedges)
        
        x, y, yerr = returnFractionGroupedBy(variant_table.loc[variant_table.numTests >= min_n], 'dG_bin', 'ci_width')
        binstart = binedges[0]
        binwidth = (binedges[1] - binedges[0])
        bincenters = np.arange(binedges[0]-binwidth*.5, binedges[-1]+2*binwidth, binwidth)
        x = bincenters[np.array(x).astype(int)]
        if ax is None:
            fig = plt.figure(figsize=(3,3))
            ax = fig.add_subplot(111)
            plt.subplots_adjust(left=0.2, bottom=0.2, top=0.95, right=0.95)
        ax.scatter(x, y, c=color, edgecolor='k', linewidth=0.5, marker=marker, s=20)
        ax.errorbar(x, y, yerr, fmt='-', elinewidth=1, capsize=0, capthick=1,
                    color=color, linestyle='', ecolor='k', linewidth=0.5)


        plt.xlabel('$\Delta G$ (kcal/mol)')
        plt.ylabel('average error (kcal/mol)')
        plt.ylim(ylim)
        plt.xlim(xlim)
        fix_axes(ax)

    def plotNumberOfMeasurments(self, xlim=[1,100], ax=None):
        """ Plot how the ci changes with dG. """
        variant_table = self.variant_table
        if ax is None:
            fig = plt.figure(figsize=(4,3))
            ax = fig.add_subplot(111)
            plt.subplots_adjust(left=0.15, bottom=0.15, top=0.95, right=0.95)
        sns.distplot(variant_table.numTests, bins=np.arange(*xlim), kde=False, hist_kws={'histtype':'stepfilled', 'linewidth':1}, color='0.5')
        plt.xlabel('number of measurements')
        plt.ylabel('number of variants')
        plt.xlim(xlim)
        fix_axes(ax)

    def plotBothErrorAndHistogram(self, xlim=[1,100], ylim_ci=None, plot_n=None, marker='>'):
        fig = plt.figure(figsize=(4,4))
        gs = gridspec.GridSpec(2,1, left=0.2, bottom=0.15)
    
        ax1 = fig.add_subplot(gs[0])
        self.plotNumberOfMeasurments(xlim=xlim, ax=ax1)
        
        ax2 = fig.add_subplot(gs[1])
        self.plotErrorByNumberofMeasurements(xlim=xlim, ylim=ylim_ci, marker=marker, ax=ax2)
        
        if plot_n is not None:
            for ax in [ax1, ax2]:
                ax.axvline(plot_n, color='k', linewidth=0.5)
        
    
    def getResultsFromVariantTable(self, other_variant_table=None):
        """Return results format for only one variant table. """
        parameters = fitting.fittingParameters()
        if other_variant_table is None:
            other_variant_table = self.variant_table.copy()
            other_variant_table.numTests = 0
        variant_tables = [self.variant_table, other_variant_table]
        offset = 0
        return getResultsFromVariantTables(variant_tables, offset)

    def fitVariant(self, variant, func, kwargs={}):
        """Fit a variant to the objective function 'func'."""
        concentrations = self.x
        fluorescence = self.getVariantBindingSeries(variant).median().rename(variant)       
        fitParameters = fitting.getInitialFitParameters(concentrations)
        if func == objfunctions.bindingCurveNonspecificTermNonlinear:
            fitParameters.loc[:, 'dG_ns'] = [-np.inf, fitParameters.loc['upperbound', 'dG'], np.inf]
        fitParameters.loc['initial', 'fmax'] = fluorescence.max()
        results = fitting.perCluster(concentrations, fluorescence, fitParameters,
                                  change_params=False, func=func)
        return results

    def redoBootstrapped(self, variant, fmaxDistObject, weighted_fit=True, fmin_float=False, func_kwargs={}, original_slope=0.0006):
        """Redo the bootstrapping step here."""
        concentrations = self.x
        variant_table = self.variant_table
        subSeries = self.getVariantBindingSeries(variant)
        fluorescenceMat = pd.concat([self.annotated_clusters, self.binding_series.astype(float)], axis=1).sort('variant_number')
        
        # get parameters
        parameters = fitting.fittingParameters(concentrations)

    def redoBootstrapped(self, variant, fmaxDistObject, weighted_fit=True, fmin_float=False, func_kwargs={}, original_slope=0.0006):
        """Redo the bootstrapping step here."""
        concentrations = self.x
        variant_table = self.variant_table
        subSeries = self.getVariantBindingSeries(variant)
        fluorescenceMat = pd.concat([self.annotated_clusters, self.binding_series.astype(float)], axis=1).sort('variant_number')
        
        # get parameters
        parameters = fitting.fittingParameters(concentrations)

        # initiate fitParameters and initial points
        initialPoints = distribution.findInitialPoints(variant_table)
        fmin_fixed = distribution.returnFminFromFluorescence(initialPoints, fluorescenceMat, parameters.mindG)
        fitParameters = fitting.getInitialFitParametersVary(concentrations)
        fitParameters.loc['vary', 'fmin'] = fmin_float
        fitParameters.loc['initial', 'fmin'] = fmin_fixed
        fitParameters.loc['initial', 'fmax'] = fmaxDistObject.getDist(1).stats(moments='m')
        
        fittype = 'binding'
        if 'fit_slope' in func_kwargs.keys():
            if func_kwargs['fit_slope']:
                fitParameters.loc[:, 'slope'] = [0, 0.0006, np.inf, True]
                fittype = 'binding_linear'
        print fitParameters
        
        # fit
        results = fitting.perVariant(concentrations, subSeries,
                                        fitParameters,
                                        fmaxDistObject,
                                        initial_points=initialPoints.loc[variant],
                                        weighted_fit=weighted_fit,
                                        kwargs=func_kwargs)
        # plot results
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        plotting.plotFitCurve(concentrations,
                            subSeries,
                            results,
                            fittype=fittype,
                            kwargs=func_kwargs,
                            ax=ax)
        names = ['dG', 'dG_lb', 'dG_ub', 'fmax', 'fmax_lb', 'fmax_ub',  'fmin', 'fmin_lb', 'fmin_ub', 'slope', 'numTests', 'pvalue', 'rsq', 'flag']
        vec = pd.Series([getValueInTable(results, name) for name in names], index=names)
        annotationText = ['dG= %4.2f (%4.2f, %4.2f)'%(vec.dG, vec.dG_lb, vec.dG_ub),
                      'fmax= %4.2f (%4.2f, %4.2f)'%(vec.fmax, vec.fmax_lb, vec.fmax_ub),
                      'fmin= %4.2f (%4.2f, %4.2f)'%(vec.fmin, vec.fmin_lb, vec.fmin_ub),
                      'Nclusters= %4.0f'%len(subSeries),
                      'average Rsq= %4.2f'%vec.rsq,
                      'slope= %4.1e'%vec.slope,
                      'weighted= %d'%weighted_fit,
                      'fmin float= %d'%fmin_float]
        ax.annotate('\n'.join(annotationText), xy=(0.05, .95), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top', fontsize=9)
        more_concentrations = np.logspace(np.log10(concentrations[0]/10.), np.log10(concentrations[-1]*10.))
        try:
            fit = fitting.bindingCurveObjectiveFunction(fitting.returnParamsFromResults(initialPoints.loc[variant]), more_concentrations, slope=original_slope)
            ax.plot(more_concentrations, fit, 'm--')
        except:
            pass
        return results
    
    def plotResiduals(self, func_kwargs={}, bounds=None):
        """For a given range of dG bounds, plot the residuals in the final and initial data."""
        variant_series = pd.concat([self.annotated_clusters, self.binding_series.astype(float)], axis=1).groupby('variant_number').median()
        variant_table = self.variant_table
        initial_points = distribution.findInitialPoints(variant_table)
        concentrations = self.x
        
        resMany = []
        if bounds is None:
            binedges = variant_table.loc[variant_table.pvalue<0.01].dG_init.quantile(np.linspace(0, 1, 20))
        else:
            binedges = bounds
        for lb, ub in zip(binedges.iloc[:-1], binedges.iloc[1:]):
        # find variants within bounds
            variants = variant_table.loc[(variant_table.dG_init < ub)&(variant_table.dG_init >= lb)].index.tolist()
            
            func = fitting.bindingCurveObjectiveFunction
            resMat = pd.melt(pd.concat([func(fitting.returnParamsFromResults( variant_table.loc[variant]), concentrations, **func_kwargs) - variant_series.loc[variant] for variant in variants], axis=1).transpose(), value_name='final').set_index('variable')
            resMatInit = pd.melt(pd.concat([func(fitting.returnParamsFromResults( initial_points.loc[variant]), concentrations, **func_kwargs) - variant_series.loc[variant] for variant in variants], axis=1).transpose(), value_name='init').set_index('variable')
            resAll = pd.melt(pd.concat([resMat, resMatInit], axis=1).reset_index(), id_vars=['variable'], var_name='fit')
            resAll.loc[:, 'dGbin'] = '%4.1f_%4.1f'%(lb, ub)
            resMany.append(resAll)
        resMany = pd.concat(resMany)
        
        sns.factorplot(x="variable", y="value", data=resMany, hue="fit", kind='box', boxprops={'linewidth':1}, whiskerprops={'linewidth':1}, capprops={'linewidth':1})
        plt.title('final')
        plt.axhline(0, color='r', linewidth=0.5)

        return resMany


       
class perFlow():
    """Class for plots that combine off rate and variant data."""
    def __init__(self, affinityData, offRate):
        self.affinityData = affinityData
        self.offRate = offRate
        self.all_variants = pd.concat([affinityData.variant_table, offRate.variant_table], axis=1).index
    
    def getGoodVariants(self, ):
        variants = (pd.concat([self.affinityData.variant_table.pvalue < 0.01,
                               self.offRate.variant_table.pvalue < 0.1,
                               self.offRate.variant_table.fmax_lb>0.2], axis=1)).all(axis=1)
        return variants
    
    def plotDeltaGDoubleDagger(self, variant=None, dG_cutoff=None, plot_on=False, params=['koff', 'dG'], variants=None):
        parameters = fitting.fittingParameters()
        koff = self.offRate.variant_table.loc[self.all_variants, params[0]].astype( float)
        dG = self.affinityData.variant_table.loc[self.all_variants, params[1]].astype(float)

        if variant is None:
            variant = 34936
        
        # find dG predicted from off and on rates
        dG_dagger = parameters.find_dG_from_Kd(koff.astype(float))
        kds = parameters.find_Kd_from_dG(dG)
        dG_dagger_on = parameters.find_dG_from_Kd((koff/kds).astype(float))
        
        # find the variants to plot based on goodness of fit and dG cutoff if given
        if variants is None:
            variants = self.getGoodVariants()
            if dG_cutoff is not None:
                variants = variants&(dG < dG_cutoff)
            
        # deced which y to plot based on user input
        x = (dG - dG.loc[variant]).loc[variants]
        if plot_on:
            y = -(dG_dagger_on - dG_dagger_on.loc[variant]).loc[variants]
        else:
            y = (dG_dagger - dG_dagger.loc[variant]).loc[variants]
        
        # plot
        fig = plt.figure(figsize=(3,3));
        ax = fig.add_subplot(111, aspect='equal')
        
        xlim = [min(x.min(), y.min()), max(x.max(), y.max())]
        im = plt.hexbin(x, y,  extent=xlim+xlim, gridsize=100, cmap='Spectral_r', mincnt=1)
        #sns.kdeplot(x, z,  cmap="Blues", shade=True, shade_lowest=False)
        slope, intercept, r_value, p_value, std_err = st.linregress(x,y)
        
        # offset

        num_variants = 100
        offset = (x - y).mean()
        
        xlim = np.array(ax.get_xlim())
        plt.plot(xlim, xlim*slope + intercept, 'k--', linewidth=1)
        if plot_on:
            plt.plot(xlim, [y.mean()]*2, 'r', linewidth=1)
        else:
            plt.plot(xlim, xlim, 'r', linewidth=1)
        plt.xlabel('$\Delta \Delta G$')
        plt.ylabel('$\Delta \Delta G_{off}\dagger$')
        
        #plt.colorbar(im)
        plt.tight_layout()
        fix_axes(ax)

        annotationText = ['slope = %4.2f '%(slope),
                          'intercept = %4.2f'%intercept,
                          'pvalue = %4.1e'%p_value
                          ]
        ax.annotate('\n'.join(annotationText), xy=(.05, .95), xycoords='axes fraction',
                    horizontalalignment='left', verticalalignment='top')
        residuals = y - (x*slope+intercept)
        return variants, residuals
    
    def plotKdVersusKoff(self, ):
        parameters = fitting.fittingParameters()
        results_off = self.offRate.variant_table.loc[self.all_variants]
        dG = self.affinityData.variant_table.dG.loc[self.all_variants]
        
        plt.figure(figsize=(4,4))
        plt.hexbin(dG, results_off.koff, yscale='log', mincnt=1, cmap='Spectral_r')
        plt.xlabel('$\Delta G$ (kcal/mol)')
        plt.ylabel('$k_{off}$ (s)')
        fix_axes(plt.gca())
        plt.tight_layout()


    def plotEquilibrationTimes(self, concentration, wait_time, initial=1E-9):
        parameters = fitting.fittingParameters()
        variants = self.getGoodVariants()
        koff = self.offRate.variant_table.loc[variants].koff
        dG = self.affinityData.variant_table.loc[variants].dG
        kds = parameters.find_Kd_from_dG(dG.astype(float))*1E-9
        
        # for meshgrid
        kobs_bounds = [-5, 0]
        koff_bounds = [-4, -2]
        xx, yy = np.meshgrid(np.linspace(*koff_bounds, num=200),
                             np.linspace(*kobs_bounds, num=200))
        
        # for each kobs = x(i,j)+y(i,j), plot fraction equilibrated
        labels = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1-1E-12])
        min_kobs = kobs_bounds[0] + np.log10(concentration/initial)
        plt.figure(figsize=(4,4))
        cs =plt.contour(xx, yy, 100*fraction_equilibrated(np.power(10, xx)+np.power(10, yy),
                                                      wait_time),
                        labels*100, colors='k', linewidths=1)
        plt.clabel(cs, inline=1, fontsize=10, fmt='%1.0f')
        
        plt.hexbin(np.log10(koff), np.log10(concentration*koff/kds), mincnt=1,
                   cmap='Spectral_r', extent=koff_bounds+[min_kobs, min_kobs+2], gridsize=150)
        plt.xlabel('log$(k_{off})$')
        plt.ylabel('log$([$flow$] k_{on})$')
        plt.title('%.2e'%(concentration*1E9))
        plt.ylim(min_kobs, min_kobs+2)
        ax = fix_axes(plt.gca())
        plt.tight_layout()
    
class compareFlow():
    def __init__(self, affinityData1, affinityData2):
        self.all_variants = np.intersect1d(*[np.array(affinityData.variant_table.index.tolist())
                                             for affinityData in [affinityData1, affinityData2]])
        self.expt1 = affinityData1
        self.expt2 = affinityData2
        #self.all_variants = pd.concat([affinityData1.variant_table, affinityData2.variant_table], axis=1).index

    def getGoodVariants(self, ):
        variants = (pd.concat([self.expt1.variant_table.pvalue < 0.01,
                               self.expt2.variant_table.pvalue < 0.01], axis=1)).all(axis=1)
        return variants

    def compareParam(self, param, log_axes=False, filter_pvalue=False, min_n=0,
                     max_dG=None, variants=None):
        """ Compare measured values for two experiments.
        
        Can apply different cutoffs. Default is to use all values that are not
        NaN in both. Can also apply 'pvalue' cutoff, min_n measurements, maxdG,
        or you can provide list of variants in which case these override any cutoffs.
        
        """
        x = self.expt1.variant_table.loc[self.all_variants, param]
        y = self.expt2.variant_table.loc[self.all_variants, param]
        
        
        if variants is None:
            variants = np.logical_not(pd.concat([x, y], axis=1).isnull()).all(axis=1)
            if filter_pvalue:
                variants = variants&self.getGoodVariants()

            if min_n > 0:
                variants = (variants&
                            pd.concat([self.expt1.variant_table.numTests >= min_n,
                                      self.expt2.variant_table.numTests >= min_n], axis=1).all(axis=1))
            if max_dG is not None:
                variants = (variants&
                            pd.concat([self.expt1.variant_table.dG <= max_dG,
                                       self.expt2.variant_table.dG <= max_dG], axis=1).all(axis=1))
        
        x = x.loc[variants]
        y = y.loc[variants]
            
            
        plt.figure(figsize=(4,4))
        if log_axes:
            plt.hexbin(x, y, cmap='Spectral_r', mincnt=1, yscale='log', xscale='log')
        else:
            plt.hexbin(x, y, cmap='Spectral_r', mincnt=1,)
        plt.xlabel('expt1 %s'%param)
        plt.ylabel('expt2 %s'%param)
        ax = fix_axes(plt.gca())
        plt.tight_layout()
        
        # plot linear fit
        xlim = np.array(ax.get_xlim())
        
        if log_axes:
            slope, intercept, r_value, p_value, std_err = st.linregress(np.log(x),np.log(y))
            plt.plot(xlim, np.power(xlim, slope)*np.exp(intercept), 'k', linewidth=1)
            plt.plot(xlim, xlim/np.mean(x)*np.mean(y), 'r:', linewidth=1)
        else:
            slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
            plt.plot(xlim, xlim*slope+intercept, 'k', linewidth=1)
            plt.plot(xlim, xlim-np.mean(x)+np.mean(y), 'r:', linewidth=1)
        
        
        return slope, r_value**2
    
    def findCombinedTable(self, offset=None, binedges=None):
        """ Combine two variant tables with some relevant info to reproducibilty. """
        if offset is None:
            offset = 0
        variant_tables = [self.expt1.variant_table.loc[self.all_variants], self.expt2.variant_table.loc[self.all_variants]]
        eminus = pd.concat([(table.dG - table.dG_lb) for table in variant_tables],
                         axis=1, keys=['rep1', 'rep2'])
        eplus  = pd.concat([(table.dG_ub - table.dG) for table in variant_tables],
                             axis=1, keys=['rep1', 'rep2'])
        
        combined = pd.concat([table.dG for table in variant_tables], axis=1,
            keys=['rep1', 'rep2'])
        
        combined.loc[:, 'difference'] = combined.rep2 - combined.rep1
        combined.loc[:, 'eplus']      = np.sqrt((eplus**2).sum(axis=1))
        combined.loc[:, 'eminus']     = np.sqrt((eminus**2).sum(axis=1))
        combined.loc[:, 'within_bound'] = (
            (combined.difference - offset - combined.eminus <= 0)&
            (combined.difference - offset + combined.eplus  >= 0))
        
        # deltaG bins
        if binedges is None:
            binedges = np.arange(-12, -5, 0.5)
        combined.loc[:, 'rep1_bin'] = np.digitize(combined.loc[:, 'rep1'], binedges)
        combined.loc[:, 'rep2_bin'] = np.digitize(combined.loc[:, 'rep2'], binedges)
        combined.loc[:, 'rep1_n'] = variant_tables[0].numTests
        combined.loc[:, 'rep2_n'] = variant_tables[1].numTests
        
        return combined

    def findOptimalOffset(self):
        """ Find the offset that minimizes the fraction different between replicates. """
        offsets = np.linspace(-1, 1, 100)
        number = pd.Series(index=offsets)
        total_number = float(len(self.findCombinedTable()))
        for offset in offsets:
            number.loc[offset] = self.findCombinedTable(offset=offset).within_bound.sum()/total_number
        
        plt.figure(figsize=(4,3));
        plt.plot(offsets, number)
        plt.xlabel('offset')
        plt.ylabel('number within bound')
        plt.tight_layout()
        fix_axes(plt.gca())
        return number.idxmax()
    

    def plotFractionNotDifferentByN(self, offset=0, xlim=[0,25], rep=1):
        """ Plot fraction not different between replicates. """ 
        
        # get data to plot
        combined = self.findCombinedTable(offset=offset)
        x, y, yerr = returnFractionGroupedBy(combined, 'rep%d_n'%rep, 'within_bound')
        
        fig = plt.figure(figsize=(4,3));
        ax = fig.add_subplot(111)
        
        ax.scatter(x, y, c='r', edgecolor='k', linewidth=0.5, marker='>')
        ax.errorbar(x, y, yerr, fmt='-', elinewidth=1, capsize=2, capthick=1,
                    color='r', linestyle='', ecolor='k', linewidth=1)
        majorLocator   = mpl.ticker.MultipleLocator(5)
        majorFormatter = mpl.ticker.FormatStrFormatter('%d')
        minorLocator   = mpl.ticker.MultipleLocator(1)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        plt.xlabel('number of measurements')
        plt.ylabel('fraction not different from replicate')
        plt.ylim(0, 1)
        plt.xlim(xlim)
        fix_axes(ax)
        plt.tight_layout()
    
    
    def plotFractionNotDifferentByDeltaG(self, offset=0, binedges=np.linspace(-12, -6, 12), min_n=0):
        combined = self.findCombinedTable(offset=offset, binedges=binedges)
        x, y, yerr = returnFractionGroupedBy(combined, 'rep1_bin', 'within_bound')
    
        fig = plt.figure(figsize=(4,3));
        ax = fig.add_subplot(111)
        
        ax.scatter(x, y, c='r', edgecolor='k', linewidth=0.5, marker='o')
        ax.errorbar(x, y, yerr, fmt='-', elinewidth=1, capsize=2, capthick=1,
                    color='r', linestyle='', ecolor='k', linewidth=1)
        binlabels = binedges
        
        plt.xlabel('$\Delta G$')
        plt.ylabel('fraction not different from replicate')
        plt.ylim(0, 1)
        plt.xticks(np.arange(1, len(binedges)+1)-0.5, ['%.2f'%i for i in binedges], rotation=90)
        plt.xlim([-0.5, len(binedges)+0.5])
        fix_axes(ax)
        plt.tight_layout()
    
    def getResults(self, offset):
        """ Combine variant tables to make results table."""
        variant_tables = [self.expt1.variant_table.loc[self.all_variants], self.expt2.variant_table.loc[self.all_variants]]
        return getResultsFromVariantTables(variant_tables, offset) 
    
            
        
class compareFluor():
    def __init__(self, fluor1, fluor2, ):
        self.fluor1 = fluor1
        self.fluor2 = fluor2
        clusters = pd.concat([fluor1, fluor2], axis=1).index

        
        # get rid of Nans and signal less than 50
        self.x = fluor1.loc[clusters]
        self.y = fluor2.loc[clusters]

        
    def getGoodClusters(self, min_fluorescence=None):
        if min_fluorescence is None:
            min_fluorescence = 50
        x = self.x
        y = self.y
        index = np.logical_not(pd.concat([x,y], axis=1).isnull().any(axis=1))
        index.loc[(x < min_fluorescence)|(y<min_fluorescence)] = False
        return index

    def plotHex(self, xlim=[0, 800], min_fluorescence=None):
        x = self.x
        y = self.y
        
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111, aspect='equal')
        im = ax.hexbin(x, y, cmap='Spectral_r', mincnt=1, bins='log', extent=xlim+xlim)
        plt.colorbar(im)
        #slope, intercept, r_value, p_value, std_err = st.linregress(x.loc[index],y.loc[index])
        
        ax.set_xlim(xlim)
        ax.set_ylim(xlim)
        
        index = self.getGoodClusters(min_fluorescence=min_fluorescence)
        #plt.plot(xlim, xlim*slope + intercept, 'k')
        slope = (y.loc[index]/x.loc[index]).median()
        plt.plot(xlim, np.array(xlim)*slope, 'k:')
        plt.plot(xlim, xlim, color='0.5', alpha=0.5 )
        plt.xlabel('signal in image 1')
        plt.ylabel('signal in image 2')
        annotation_text = ('slope (median)=%4.2f'%(slope))
        plt.annotate(annotation_text, xy=(.05, .95), xycoords='axes fraction',
                        horizontalalignment='left', verticalalignment='top')
        plt.tight_layout()
        fix_axes(ax)
        return x, y, index, 
    
    def bootstrapSlope(self, min_fluorescence=None, log_axis=False):
        x = self.x
        y = self.y
        index = self.getGoodClusters(min_fluorescence=min_fluorescence)
        
        if log_axis:
            slopes = np.log(y.loc[index]/x.loc[index])
            bins = np.linspace(-3, 3, 100)
            xlim = [-3, 3]
        else:
            slopes = y.loc[index]/x.loc[index]
            bins = np.linspace(0, 5, 100)
            xlim = [0, 2]       
        
        plt.figure(figsize=(4,4));
        sns.distplot(slopes, bins=bins, kde_kws={'clip':xlim},
                     hist_kws={'histtype':'stepfilled'}, color='0.5')
        
        plt.axvline(slopes.median(), color='k', linewidth=1)
        
        # find probability distribution max
        digitized = np.digitize(slopes, bins)
        mode = st.mode(digitized).mode[0]
        max_prob = (bins[mode] + bins[mode-1])*0.5
        plt.axvline(max_prob, color='r', linewidth=1, linestyle=':', label='max probability')

        annotation_text = ('slope (median)=%4.3f\nslope (max prob)=%4.5f'%(slopes.median(), max_prob))
        plt.annotate(annotation_text, xy=(.05, .95), xycoords='axes fraction',
                        horizontalalignment='left', verticalalignment='top')
        
        #lb, ub = bootstrap.ci(slopes, statfunction=np.median, n_samples=1000)
        #for bound in [lb, ub]:
        #    plt.axvline(bound, color='k', linewidth=1, linestyle=':')
        plt.xlabel('slope')
        plt.ylabel('probability distribution')
        plt.xlim(xlim)
        plt.tight_layout()
        fix_axes(plt.gca())
        return slopes
    
    def plotSlopeVersusSignal(self, min_fluorescence=None):
        x = self.x
        y = self.y
        index = self.getGoodClusters(min_fluorescence=min_fluorescence)
        
        slopes = y.loc[index]/x.loc[index]
        mean_signal = (x + y).loc[index]/2.
        
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111,)
        im = ax.hexbin(mean_signal, np.log2(slopes), cmap='Spectral_r', mincnt=1, bins='log')
        plt.xlabel('mean signal per cluster')
        plt.ylabel('signal in image2/image1')
        plt.tight_layout()
        fix_axes(plt.gca())
    
    def findAlpha(self, n):
        x = self.x
        y = self.y
        index = (x>0)&(y>0)
        #vec = np.exp(np.log(y/x)/n).loc[index]
        vec = np.power(y/x, 1/float(n)).loc[index]
        alpha = vec.median()
        lb, ub = bootstrap.ci(vec, statfunction=np.median, n_samples=1000)
        xlim = [0.9 ,1.1]
        bins = np.arange(0.898, 1.1, 0.001)
        plt.figure(figsize=(4,3));
        sns.distplot(vec, bins=bins, hist_kws={'histtype':'stepfilled'}, kde_kws={'clip':xlim});
        plt.axvline(alpha, color='k', linewidth=0.5);
        plt.axvline(lb, color='k', linestyle=':', linewidth=0.5);
        plt.axvline(ub, color='k', linestyle=':', linewidth=0.5);
        plt.xlabel('photobleach fraction per image');
        plt.ylabel('probability');
        plt.xlim(xlim)
        fix_axes(plt.gca());
        plt.tight_layout()
        return alpha, lb, ub

class manyFlows():
    """Class to store many different perVariant classes"""
    def __init__(self, perVariants, names=None):
        self.flowData = perVariants
        if names is None:
            names = np.arange(len(perVariants))
        self.names = names
        
    def plotAllInitVsFinal(self, param='dG', colorby='fmax_init', limits=[-12.5, -6], order=None, **kwargs):
        """For each of the variant tables, plot the dG init versus dG final"""
        alldata = []
        for name, data in self.flowData.items():
            subdata = data.variant_table.loc[:, ['%s_init'%param, param]]
            subdata.loc[:, 'expt'] = name
            subdata.loc[:, 'colorby'] = data.variant_table.loc[:, colorby]
            alldata.append(subdata)
        alldata = pd.concat(alldata)      

        if order is None: order = self.names

        cmap = sns.diverging_palette(220, 10, center="dark", as_cmap=True)
        g = sns.FacetGrid(alldata, col="expt", col_order=order, col_wrap=4, xlim=limits, ylim=limits)
        g.map(my_scatterfun, '%s_init'%param, param, "colorby", cmap=cmap,
              edgecolor='none', marker='.', **kwargs)
            

def my_scatterfun(x, y, c, **kwargs):
    plt.scatter(x, y, c=c, **kwargs)

def fraction_equilibrated(kobs, time_waited):
    return 1-np.exp(-kobs*time_waited)
 

    
def _makebootstrappeddist(vec, n_samples=1000, statfunction=np.median):
    indices = np.random.choice(vec.index, [n_samples, len(vec)])
    bootstrapped_vec = []
    for index in indices:
        bootstrapped_vec.append(statfunction(vec.loc[index]))
    return np.array(bootstrapped_vec)

                            
def returnFractionGroupedBy(mat, param_in, param_out):
    """ Given data matrix, return mean of param_out binned by param_in. """
    grouped = mat.groupby(param_in)[param_out]
    y = grouped.mean()
    x = y.index.tolist()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yerr = np.array([np.abs(bootstrap.ci(group, method='pi', n_samples=100) - y.loc[name])
                         for name, group in grouped]).transpose()
    return x, y, yerr   

def getResultsFromVariantTables(variant_tables, offset=0, param='dG', include_fmax=False):
    """combine variant tables to form a results table."""
    # make flags
    # 0x1 data from rep 1
    # 0x2 data from rep 2
    # 0x4 data not measured in rep 1
    # 0x8 data not measured in rep 2
    parameters = fitting.fittingParameters()
    

    values = pd.concat([table.loc[:, param] for table in variant_tables],
                       axis=1, keys=['0', '1']).astype(float)
    fmax = pd.concat([table.loc[:, 'fmax'] for table in variant_tables],
                       axis=1, keys=['0', '1']).astype(float)    
    numTests = pd.concat([table.numTests for table in variant_tables],
                       axis=1, keys=['0', '1'])
    numTests.fillna(0, inplace=True)
    
    # correct for offset between 1st and 2nd measurement
    values.iloc[:, 1] -= offset
    
    # find error measurements
    eminus = pd.concat([(table.loc[:, param] - table.loc[:, '%s_lb'%param]) for table in variant_tables],
                         axis=1, keys=values.columns)
    eplus  = pd.concat([(table.loc[:, '%s_ub'%param] - table.loc[:, param]) for table in variant_tables],
                         axis=1, keys=values.columns)
    variance = (((eminus + eplus)/2)**2).astype(float)
    weights = pd.DataFrame(data=1, index=variance.index, columns=variance.columns).astype(float)
    index = (variance > 0).all(axis=1)
    weights.loc[index] = 1/variance.loc[index]
    
    # final variant table
    if include_fmax:
        cols =  ['numTests1', 'numTests2', 'weights1', 'weights2', 'numTests', param, 'eminus', 'eplus', 'fmax', 'flag']
    else:
        cols =  ['numTests1', 'numTests2', 'weights1', 'weights2', 'numTests', param, 'eminus', 'eplus', 'flag']
        
    results = pd.DataFrame(index=pd.concat(variant_tables, axis=1).index, columns=cols)
    results.loc[:, ['numTests1', 'numTests2', 'numTests', 'flag']] = 0
    
    # if one of the measurements has less than 5 clusters, use only the other measurement
    indexes = [(numTests.iloc[:, 0] < 5)&(numTests.iloc[:, 1] >= 5),
               (numTests.iloc[:, 0] >= 5)&(numTests.iloc[:, 1] < 5),
               (numTests >= 5).all(axis=1)]
    weights.loc[indexes[0], '0'] = 0
    weights.loc[indexes[1], '1'] = 0

    for i, index in enumerate(indexes):
        results.loc[index, param]     = weightedAverageAll(values, weights, index=index)
        results.loc[index, 'eminus'] = errorPropagateAverageAll(eminus, weights, index=index)
        results.loc[index, 'eplus']  = errorPropagateAverageAll(eplus, weights, index=index)
        if include_fmax:
            results.loc[index, 'fmax']  = weightedAverageAll(fmax, weights, index=index)
        
        results.loc[index, 'flag']  += np.power(2, i)
        results.loc[index, ['numTests1', 'numTests2']] = numTests.loc[index].values
        results.loc[index, ['weights1', 'weights2']]   = weights.loc[index].values
        results.loc[index, 'numTests'] = weightedAverageAll(numTests, weights, index=index)
    return results

def annotateBindingCurve(vec, ax):
    """Given a pd.Series of things ot annotate, annotate them onto ax"""
    annotationText = ['dG= %4.2f (%4.2f, %4.2f)'%(vec.dG,
                                                    vec.dG_lb,
                                                    vec.dG_ub),
                      'fmax= %4.2f (%4.2f, %4.2f)'%(vec.fmax,
                                                    vec.fmax_lb,
                                                    vec.fmax_ub),
                      'Nclusters= %4.0f'%vec.numTests,
                      'pvalue= %.1e'%vec.pvalue,
                      'average Rsq= %4.2f'%vec.rsq,
                      'fmax enforced=%s'%str(vec.flag)
                      ]
    ax.annotate('\n'.join(annotationText), xy=(0.05, .95), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top', fontsize=9)
    return ax

def annotateOffrate(vec, ax):
    """Given a pd.Series of things ot annotate, annotate them onto ax"""
    annotationText = ['koff= %4.2e (%4.2e, %4.2e)'%(vec.koff,
                                                    vec.koff_lb,
                                                    vec.koff_ub),
                      'fmax= %4.2f (%4.2f, %4.2f)'%(vec.fmax,
                                                    vec.fmax_lb,
                                                    vec.fmax_ub),
                      'Nclusters= %4.0f'%vec.numTests,
                      'pvalue= %.1e'%vec.pvalue,
                      'average Rsq= %4.2f'%vec.rsq,
                      ]
    ax.annotate('\n'.join(annotationText), xy=(.05, .95), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top', fontsize=9)

    return ax
