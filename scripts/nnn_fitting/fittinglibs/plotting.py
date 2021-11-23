"""
Sarah Denny
Stanford University

"""

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import warnings
import argparse
import itertools  
import seaborn as sns
import scipy.spatial.distance as ssd
import scipy.stats as st
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib as mpl
from matplotlib import gridspec
from joblib import Parallel, delayed
sns.set_style("white", {'xtick.major.size': 4,  'ytick.major.size': 4,
                        'xtick.minor.size': 2,  'ytick.minor.size': 2,
                        'lines.linewidth': 1})
from  fittinglibs import objfunctions
import fittinglibs.fitting as fitting

def fix_axes(ax):
    ax.tick_params(which='minor', top='off', right='off')
    ax.tick_params(top='off', right='off', pad=2, labelsize=10, labelcolor='k')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ax

def savefig(filename, **kwargs):
    """save fig"""
    plt.savefig(filename, **kwargs)

def annotate_axes(s, ax=None, **kwargs):
    """annotate the current figure with string 's' """
    if ax is None:
        ax = plt.gca()
    if 'xy' not in list(kwargs.keys()):
        kwargs['xy'] = (0.05, 0.95)
    if 'xycoords' not in list(kwargs.keys()):
        kwargs['xycoords'] = 'axes fraction'
    if 'horizontalalignment' not in list(kwargs.keys()):
        kwargs['horizontalalignment'] = 'left'
    if 'verticalalignment' not in list(kwargs.keys()):
        kwargs['verticalalignment'] = 'top'
    if 'fontsize' not in list(kwargs.keys()):
        kwargs['fontsize'] = 10
    ax.annotate(s, **kwargs)    

def get_c(x, y, distance_threshold=None):
    """ Given two arrays x and y, return the number of points within a certain distance of each point."""
    if distance_threshold is None:
        distance_threshold = min(x.std(), y.std())/10
    distance_mat = ssd.squareform(ssd.pdist(pd.concat([x, y], axis=1)))
    c = ((distance_mat < distance_threshold).sum(axis=1) - 1)/2
    c = (c-c.min())/(c.max()-c.min()).astype(float)
    return c


def my_smoothed_scatterplot(x,y, distance_threshold=None, color=None,**kwargs):
    """ given x and y, plot a scatterplot with color according to density."""
    c = get_c(x, y, distance_threshold=distance_threshold)
    if 'cmap' not in list(kwargs.keys()):
        if color is not None:
            cmap = sns.dark_palette(color, as_cmap=True)
        else:
            cmap = None
        plt.scatter(x, y, c=c, cmap=cmap,
                edgecolors='none', marker='.', rasterized=True, **kwargs)
    else:
        plt.scatter(x, y, c=c, 
                edgecolors='none', marker='.', rasterized=True, **kwargs)

    return



def plotDataErrorbars(x, subSeries, ax=None, capsize=2, errors=None):
    """ Find errorbars on set of cluster fluorescence and plot. """
    
    # set errors to NaN unless successfully find them later on with bootstrapping
    if errors is None:
        use_default=False
        default_errors = [np.ones(len(x))*np.nan]*2 
    else:
        use_default = True
        default_errors = errors
    
    # if subseries is a dataframe, find errors along columns. if verctor, no errors will be found.
    if len(subSeries.shape) == 1:
        fluorescence = subSeries
        numTests = np.array([1 for col in subSeries])
    else:
        fluorescence = subSeries.median()
        numTests = np.array([len(subSeries.loc[:, col].dropna()) for col in subSeries])
    
    # option to use only default errors provdided for quicker runtime
    if not use_default:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            eminus, eplus = fitting.findErrorBarsBindingCurve(subSeries)
    else:
        eminus, eplus = default_errors
    
    # if ax is given, plot to that. else
    if ax is None:
        ax = plt.gca()
    
    # plot binding points
    ax.errorbar(x, fluorescence,
                 yerr=[eminus, eplus], fmt='.', elinewidth=1,
                 capsize=capsize, capthick=1, color='k', linewidth=1)
    return ax

def plotFit(x, params, func, ax=None, fit_kwargs=None):
    """ given x values, lmfit params,  and fitting function, calculate fit and plot to current axes."""
    if fit_kwargs is None:
        fit_kwargs = {}
    fit = func(params, x, **fit_kwargs)
    
    if ax is None:
        ax = plt.gca()
    ax.plot(x, fit, 'r')
    return ax

def plotFitBounds(x, params_lb, params_ub, func, ax=None, fit_kwargs=None):
    """ given x values, lmfit params,  and fitting function, calculate fit and plot to current axes."""
    if fit_kwargs is None:
        fit_kwargs = {}
    ub = func(params_ub, x, **fit_kwargs)
    lb = func(params_lb, x, **fit_kwargs)
    
    # plot upper and lower bounds
    if ax is None:
        ax = plt.gca()
    ax.fill_between(x, lb, ub, color='0.5',
                         label='95% conf int', alpha=0.5)
    return ax

def plotFitCurve(x, subSeries, results, param_names=None, ax=None, log_axis=True, capsize = 2,
                 func=objfunctions.binding_curve, fittype='binding', kwargs=None, errors=None):
    if kwargs is None:
        kwargs = {}

    # these are useful definitions for three commonly used fitting functions
    further_process=True
    if fittype == 'binding':
        param_names_tmp = ['fmax', 'dG', 'fmin']
        ub_vec = ['_ub', '_lb', '']
        lb_vec = ['_lb', '_ub', '']
        log_axis = True
        xlabel = 'concentration (nM)'
    elif fittype == 'off':
        param_names_tmp = ['fmax', 'koff', 'fmin']
        ub_vec = ['_ub', '_lb', '_ub']
        lb_vec = ['_lb', '_ub', '_lb']
        log_axis = False
        xlabel = 'time (s)'
    elif fittype == 'on':
        param_names_tmp = ['fmax', 'kobs', 'fmin']
        ub_vec = ['_ub', '_ub', '_ub']
        lb_vec = ['_lb', '_lb', '_lb']
        log_axis=False
        xlabel = 'time (s)'
    elif fittype == 'binding_linear':
        param_names_tmp = ['fmax', 'dG', 'fmin', 'slope']
        ub_vec = ['_ub', '_ub', '_ub', '_ub']
        lb_vec = ['_lb', '_lb', '_lb', '_lb']
        log_axis = True
        xlabel = 'concentration (nM)'
    else:
        further_process = False
        xlabel = ''
        if param_names is None:
            print('Need to define param names or fittype. fittype %s not recognized! Exiting.'%(fittype))
            sys.exit()
        
    # allow custom definition of param_names with fitParameters 
    if param_names is None:
        param_names = param_names_tmp
    
    # get params for fit function
    params = fitting.returnParamsFromResults(results, param_names)
    
    # gerenate x values for fit function
    if ax is None:
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
    
    # generate x values for fit function
    x = np.array(x)
    if log_axis:
        more_x = np.logspace(np.log10(x.min()/10), np.log10(x.max()*2), 100)
    else:
        more_x = np.linspace(x.min(), x.max(), 100)

    # plot the data
    plotDataErrorbars(x, subSeries, ax, capsize=capsize, errors=errors)

    # plot fit
    plotFit(more_x, params, func, ax=ax, fit_kwargs=kwargs)
    
    if further_process:
        # plot upper and lower bounds
        all_param_names = [['%s%s'%(param, s) for param, s in zip(param_names, vec)]
                           for vec in [ub_vec, lb_vec]]
        if np.all(np.in1d(all_param_names, results.index.tolist())):
            params_ub = fitting.returnParamsFromResultsBounds(results, param_names, ub_vec)
            params_lb = fitting.returnParamsFromResultsBounds(results, param_names, lb_vec)
            plotFitBounds(more_x, params_lb, params_ub, func, ax=ax, fit_kwargs=kwargs)

    # format
    ylim = ax.get_ylim()
    xlim = more_x[[0,-1]]
    plt.ylim(0, ylim[1])
    plt.xlim(xlim)
    plt.xlabel(xlabel)
    if log_axis:
        plt.xscale('log')
    else:
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ylabel('normalized fluorescence')
    fix_axes(ax)
    plt.tight_layout()
    return



def plotFmaxVsKd(variant_table, cutoff, subset=None, kde_plot=False,
                 plot_fmin=False, xlim=None, ylim=None):
    if subset is None:
        if kde_plot: subset = True
        else: subset = False
    
    parameters = fitting.fittingParameters()
    
    kds = parameters.find_Kd_from_dG(variant_table.dG)
    if plot_fmin:
        fmax = variant_table.fmin
        ylabel_text = 'min'
    else:
        fmax = variant_table.fmax
        ylabel_text = 'max'
        
    # find extent in x
    if xlim is None:
        log_kds = np.log10(kds)
        kds_bounds = [np.floor(log_kds.min()), log_kds.median() + log_kds.std()*3]
    else:
        kds_bounds = [np.log10(i) for i in xlim]
    
    #find extent in y
    if ylim is None:
        fmax_bounds = [0, fmax.median() + 3*fmax.std()]
    else:
        fmax_bounds = ylim
    
    # initiate plot
    plt.figure(figsize=(3,3))
    ax = plt.gca()
    if subset:
        x = kds.iloc[::100]
        y = fmax.iloc[::100]
    else:
        x = kds
        y = fmax
    
    if kde_plot:   
        sns.kdeplot(np.log10(x),
                    y,
                    shade=True, shade_lowest=False, n_levels=20, clip=[kds_bounds, fmax_bounds],
                    cmap="binary")
        xticks = ax.get_xticks()
        ax.set_xticklabels(['$10^%d$'%x for x in xticks])
        cutoff = np.log10(cutoff)
    else:
        ax.hexbin(x,
                  y, xscale='log',
                  extent=np.hstack([kds_bounds, fmax_bounds]),
                  cmap="Spectral_r", mincnt=1)
    
    fix_axes(ax)
    plt.xlabel('$K_d$ (nM)')
    
    plt.ylabel('initial $f_{%s}$'%ylabel_text)
    ylim=ax.get_ylim()
    plt.plot([cutoff]*2, ylim, 'r:', label='cutoff for 95% bound')
    plt.tight_layout()

def plotFmaxStdeVersusN(fmaxDist, stds_object, maxn, ax=None):
    # plot
    x = stds_object.index.tolist()
    y = stds_object.loc[:, 'std']
    params = fmaxDist.params
    x_fit = np.arange(1, maxn)
    y_fit = fmaxDist.sigma_by_n_fit(params, x_fit)
    if ax is None:
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        marker='o'
        color='k'
        linestyle='-'
        linecolor = 'c'
    else:
        marker='.'
        color='0.5'
        linestyle=':'
        linecolor =color     
    ax.scatter(x, y, s=10, marker=marker, color=color)
    ax.plot(x_fit, y_fit, linestyle=linestyle, color=linecolor)
    plt.xlabel('number of measurements')
    plt.ylabel('standard deviation of median fmax')
    plt.xlim(0, maxn)
    fix_axes(ax)
    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)
    return ax

def plotFmaxOffsetVersusN(fmaxDist, stds_object, maxn, ax=None):
    x = stds_object.index.tolist()
    y = stds_object.offset.values
    yerr = stds_object.offset_stde.values

    if ax is None:
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        marker='o'
        color='k'
        linestyle='-'
        linecolor = 'c'
    else:
        marker='.'
        color='0.5'
        linestyle=':'
        linecolor =color    

    plt.errorbar(x, y, yerr=yerr, fmt='.', color=color, marker=marker, markersize=3)
    plt.axhline(0, color=linecolor, linestyle=linestyle)
    #plt.plot(x, y_smoothed, 'c')
    plt.xlabel('number of measurements')
    plt.ylabel('offset')
    plt.xlim(0, maxn)
    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)
    fix_axes(plt.gca())
    return

def plotNumberVersusN(n_tests, maxn):
    plt.figure(figsize=(4, 3))
    sns.distplot(n_tests, bins=np.arange(maxn+1), kde=False, color='grey',
                 hist_kws={'histtype':'stepfilled'})
    plt.xlabel('number of measurements')
    plt.ylabel('count')
    plt.xlim(0, maxn)
    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)
    fix_axes(plt.gca())
    

def plotFmaxInit(variant_table):
    fmax_subset = variant_table.loc[~variant_table.flag.astype(bool)].fmax
    bounds = [fmax_subset.median()-3*fmax_subset.std(), fmax_subset.median() + 3*fmax_subset.std()]
    parameters = fitting.fittingParameters()
    detection_limit = -6.93
    cmap = sns.diverging_palette(220, 20, center="dark", as_cmap=True)
    index = variant_table.loc[variant_table.numClusters >= 5].index
    
    x = parameters.find_Kd_from_dG(variant_table.loc[index].dG_init.astype(float))
    y = parameters.find_Kd_from_dG(variant_table.loc[index].dG.astype(float))
    
    xlim = [1, 1E5]
    fig = plt.figure(figsize=(4.5,3.75))
    ax = fig.add_subplot(111, aspect='equal')
    im = ax.scatter(x, y, marker='.', alpha=0.5,
                    c=variant_table.loc[index].fmax_init, vmin=bounds[0], vmax=bounds[1], cmap=cmap, linewidth=0)
    plt.plot(xlim, xlim, 'c:', linewidth=1)
    plt.plot([detection_limit]*2, xlim, 'r:', linewidth=1)
    plt.plot(xlim, [detection_limit]*2, 'r:', linewidth=1)
    plt.xlim(xlim); plt.xlabel('$K_d$ initial (kcal/mol)')
    plt.ylim(xlim); plt.ylabel('$K_d$ final (kcal/mol)')
    plt.colorbar(im, label='fmax initial')
    ax.tick_params(top='off', right='off')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()
    return

def plotBoundFluorescence(signal, bounds):
    """ Plot histogram of all RNA fluorescence and bounds imposed on distribution."""
    lowerbound, upperbound  = bounds
    binwidth = (upperbound - lowerbound)/50.
    plt.figure(figsize=(4,3))
    sns.distplot(signal.dropna(), bins = np.arange(signal.min(), signal.max()+binwidth, binwidth), color='seagreen')
    ax = plt.gca()
    ylim = ax.get_ylim()
    plt.plot([lowerbound]*2, ylim, 'k:')
    plt.plot([upperbound]*2, ylim, 'k:')
    plt.xlim(0, upperbound + 2*signal.std())
    plt.xlabel('all cluster fluorescence')
    plt.ylabel('probability density')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    fix_axes(ax)
    plt.tight_layout()
    
def plotGammaFunction(vec, func, results=None, params=None, bounds=None):
    """ Take vector and fit and plot distribution. """
    # plot pdf
    if bounds is None:
        bounds = np.percentile(vec, [0,100])
    plotDist(vec, bounds)
    more_x = np.linspace(*bounds, num=100)

    if results is None and params is None:
        print("need to define either results or params to plot fit")
    else:
        if params is None:
            params = fitting.returnParamsFromResults(results, param_names=['mean', 'std', 'offset'])
        plt.plot(more_x, func(params, more_x, return_pdf=True), 'r')
        
def plotAnyN(tight_binders, fmaxDistObject, n, bounds):
    """Plot a particular std given the number of measurements."""
    x = np.linspace(*bounds, num=100)
    plotDist(tight_binders.loc[tight_binders.numTests==n].fmax, bounds)
    plt.plot(x, fmaxDistObject.getDist(n).pdf(x), 'r')
    plt.tight_layout()
    
    
def plotDist(vec, bounds):
    """Plot the disttribution of fmax initial."""
    plt.figure(figsize=(3.5,3))
    sns.distplot(vec, hist_kws={'histtype':'stepfilled'}, color='0.5', kde_kws={'clip':bounds})
    ax = fix_axes(plt.gca())
    xlim = ax.get_xlim()
    plt.xlim(0, xlim[1])
    plt.xlabel('initial $f_{max}$')
    plt.ylabel('probability density')
    plt.subplots_adjust(bottom=0.2, left=0.2, top=0.95, right=0.95)
    return

def plotErrorInBins(variant_table, xdelta=None):
    parameters = fitting.fittingParameters()
    variant_table = variant_table.astype(float)
    errors = variant_table.dG_ub - variant_table.dG_lb
    plotErrorBars(parameters.find_Kd_from_dG(variant_table.dG),
                  variant_table.numTests, errors, xdelta=xdelta)
    return

def plotPercentErrorInBins(variant_table, xdelta=None):
    parameters = fitting.fittingParameters()
    variant_table = variant_table.astype(float)
    errors = ((parameters.find_Kd_from_dG(variant_table.dG_ub) -
               parameters.find_Kd_from_dG(variant_table.dG_lb))/
               parameters.find_Kd_from_dG(variant_table.dG))*100
    ylim = [0, 250]
    yticks = np.arange(0, 250, 50)
    ylabel = 'percent error on Kd'
    plotErrorBars(parameters.find_Kd_from_dG(variant_table.dG),
                  variant_table.numTests, errors,
                  ylim=ylim, yticks=yticks,
                  ylabel=ylabel, xdelta=xdelta)

def plotErrorBars(kds, numTests, errors, ylim=None, yticks=None, ylabel=None, xdelta=None):

    binedges = np.power(10., [0, 1, 2, 3, 4, 5])
    binned_Kds = np.digitize(kds, binedges)


    
    if ylim is None:
        ylim = [0, 1.5]
    if yticks is None:
        yticks = np.arange(0, 1.5, 0.25)
    if ylabel is None:
        ylabel = 'confidence interval width (kcal/mol)'
    if xdelta is None:
        xdelta = 5

    numbers = np.unique(numTests)
    xticks = np.arange(0, len(numbers), xdelta)[1:]
    xticklabels = ['%d'%n  for n in xticks ]

    fig = plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(len(binedges)-1, 1)
    for i, bin_idx in enumerate(np.arange(1, len(binedges))):
        index = binned_Kds==bin_idx
        ax = fig.add_subplot(gs[i])
        sns.barplot(x=numTests.loc[index],
                    y=errors.loc[index],
                    ax=ax,
                    order=numbers,
                    color="r",
                    edgecolor="r",
                    error_kw={'elinewidth':0.5})
        ax.set_xticks(xticks)
        if bin_idx == len(binedges)-1:
            ax.set_xticklabels(xticklabels)
            ax.set_xlabel('number of tests')
        else:
            ax.set_xticklabels('')
            ax.set_ylabel('')
            
        ax.set_ylim(ylim)
        ax.set_yticks(yticks)

        ax.tick_params(top='off', right='off')
        ax.annotate('$10^%d \leq K_d \leq 10^%d$ nM'%(np.log10(binedges[i]),
                                                   np.log10(binedges[i+1])),
                    xy=(0.95, 0.95),
                    xycoords='axes fraction',
                    horizontalalignment='right', verticalalignment='top',
                    fontsize=12)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    plt.annotate(ylabel, rotation=90,
                 xy=(0.05, 0.5),
                 xycoords='figure fraction',
                 horizontalalignment='left', verticalalignment='center',
                 fontsize=12)
    return


        
def plotNumberInBins(variant_table, xdelta=None):
    if xdelta is None:
        xdelta=5
    parameters = fitting.fittingParameters()
    kds = parameters.find_Kd_from_dG(variant_table.dG.astype(float))
    numTests = variant_table.numTests
    
    binedges = np.power(10., [0, 1, 2, 3, 4, 5])
    binned_Kds = np.digitize(kds, binedges)

    numbers = np.unique(numTests)
    xticks = np.arange(0, len(numbers), xdelta)[1:]
    xticklabels = ['%d'%n  for n in xticks ]
    ylabel = 'number of variants'
        
    fig = plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(len(binedges)-1, 1)
    for i, bin_idx in enumerate(np.arange(1, len(binedges))):
        index = binned_Kds==bin_idx
        ax = fig.add_subplot(gs[i])
        plt.hist(numTests.loc[index].values,
                 bins=numbers,
                 color=sns.xkcd_rgb['charcoal'],
                 rwidth=0.8,
                 alpha=0.5)
        ax.set_xticks(xticks)
        ax.set_xlim(0, xticks[-1]+1)
        if bin_idx == len(binedges)-1:
            ax.set_xticklabels(xticklabels)
            ax.set_xlabel('number of tests')
        else:
            ax.set_xticklabels('')
            ax.set_ylabel('')
        #
        #ylim = ax.get_ylim()
        ##ax.set_ylim(ylim)
        #delta = np.ceil(ylim[1]/100.)
        #ax.set_yticks(np.arange(0, delta*100, np.around(delta/4.)*100))

        ax.tick_params(top='off', right='off')
        ax.annotate('$10^%d \leq K_d \leq 10^%d$ nM'%(np.log10(binedges[i]),
                                                   np.log10(binedges[i+1])),
                    xy=(0.95, 0.95),
                    xycoords='axes fraction',
                    horizontalalignment='right', verticalalignment='top',
                    fontsize=12)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    plt.annotate(ylabel, rotation=90,
                 xy=(0.05, 0.5),
                 xycoords='figure fraction',
                 horizontalalignment='left', verticalalignment='center',
                 fontsize=12)
    return

def plotNumberTotal(variant_table, binedges=None, variant_table2=None):
    bins = np.arange(120)
    fig = plt.figure(figsize=(3.5,2.5))
    hist, binedges, patches = plt.hist(variant_table.numTests.astype(float).values,
                            bins=bins, histtype='stepfilled', alpha=0.5, color=sns.xkcd_rgb['navy blue'])
    plt.plot((binedges[:-1]+binedges[1:])*0.5, hist, linewidth=1, alpha=0.1, color=sns.xkcd_rgb['navy blue'])
    #sns.distplot(variant_table.numTests.astype(float).values, bins=bins,
    #            hist_kws={'histtype':'stepfilled'}, color='grey')
    plt.xlabel('# tests')
    plt.ylabel('# variants')
    plt.tight_layout()
    ax = plt.gca()
    ax.tick_params(top='off', right='off')

    if variant_table2 is not None:
        hist, binedges, patches = plt.hist(variant_table2.numTests.astype(float).values,
                                           bins=bins, histtype='stepfilled', alpha=0.5, color=sns.xkcd_rgb['dark cyan'])
        plt.plot((binedges[:-1]+binedges[1:])*0.5, hist, linewidth=1, alpha=0.1, color=sns.xkcd_rgb['dark cyan'])

    ylim = ax.get_ylim()
    plt.plot([5]*2, ylim, 'k--', linewidth=1, alpha=0.5)

def plotFractionFit(variant_table, binedges=np.arange(-12, -6, 0.5), param='dG', pvalue_threshold=0.05):
    # plot
    binwidth=0.01
    bins=np.arange(0,1+binwidth, binwidth)
    plt.figure(figsize=(4, 3.5))
    plt.hist(variant_table.loc[variant_table.pvalue <= pvalue_threshold].fitFraction.values,
             alpha=0.5, color='red', bins=bins, label='passing cutoff')
    plt.hist(variant_table.loc[variant_table.pvalue > pvalue_threshold].fitFraction.values,
             alpha=0.5, color='grey', bins=bins,  label='fails cutoff')
    plt.ylabel('number of variants')
    plt.xlabel('fraction fit')
    plt.legend(loc='upper left')
    plt.tight_layout()
    fix_axes(plt.gca())    

    subtable = pd.DataFrame(index=variant_table.index,
                            columns=['binned_dGs', 'pvalueFilter'],
                            data=np.column_stack([np.digitize(variant_table.loc[:, param], binedges),
                                                  variant_table.pvalue <= pvalue_threshold]))
    g = sns.factorplot(x="binned_dGs", y="pvalueFilter", data=subtable,
                order=np.unique(subtable.binned_dGs),
                color="r", kind='bar');
    g.set(ylim=(0, 1.1), ylabel='fraction pass pvalue cutoff');
    g.set_xticklabels(['<%4.1f'%binedges[0]] +
        ['%3.1f:%3.1f'%(i, j) for i, j in zip(binedges[:-1], binedges[1:])]+
        ['>%4.1f'%binedges[-1]], rotation=90)
    g.set(xticks=np.arange(len(binedges)+1))
    g.fig.subplots_adjust(hspace=.2, bottom=0.35)
    fix_axes(plt.gca())


