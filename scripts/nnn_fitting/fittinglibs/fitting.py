from lmfit import minimize, Parameters, report_fit, conf_interval
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scikits.bootstrap import bootstrap
from statsmodels.distributions.empirical_distribution import ECDF
import warnings
import itertools
import scipy.stats as st
import copy
import datetime
from fittinglibs import variables
from variables import fittingParameters

def getFitParam(param, concentrations=None, init_val=None, vary=None, ub=None, lb=None):
    """For a given fit parameter, return reasonable lowerbound, initial guess, and upperbound.
    
    Different inputs may be provided.
    - For param='dG', need to provide the concentrations, (i.e. dG_params = {'concentrations':[1,2,4,8]} so that reasonable bounds can be found on the dG.
    - can specify 'vary', which will be a bool by default set to true, in addition to the lb, initial, and ub.
    - For param='fmin', may want to provide the 
    """
    fitParam = pd.Series(index=['lowerbound', 'initial', 'upperbound'], name=param)
    if param=='dG':
        if concentrations is None:
            print 'must specify concentrations to find initial parameters for dG'
            return fitParam
        parameters = fittingParameters(concentrations=concentrations)
        fitParam.loc['lowerbound'] = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.99, concentrations[0]))
        fitParam.loc['initial'] = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.5, concentrations[-1]))
        fitParam.loc['upperbound'] = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.01, concentrations[-1]))
    elif param=='dGns':
        parameters = fittingParameters(concentrations=concentrations)
        fitParam.loc['lowerbound'] = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.99, concentrations[0]))
        fitParam.loc['initial'] = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.01, concentrations[-1]))
        fitParam.loc['upperbound'] = np.inf
        
    elif param=='fmin':
        fitParam.loc[:] = [0, 0, np.inf]
    elif param=='fmax':
        fitParam.loc[:] = [0, np.nan, np.inf]
    elif param=='slope':
        fitParam.loc[:] = [0, 3E-4, np.inf]

    else:
        print 'param %s not recognized.'%param
    
    # change vary
    if vary is not None:
        fitParam.loc['vary'] = vary
        
    # change init param
    if init_val is not None:
        fitParam.loc['initial'] = init_val
        
    # change ub
    if ub is not None:
        fitParam.loc['upperbound'] = ub
    if lb is not None:
        fitParam.loc['lowerbound'] = lb        
    return fitParam
    
def getInitialFitParameters(concentrations):
    """ Return fitParameters object with minimal constraints.
    
    Input: concentrations
    Uses concencetration to provide constraints on dG
    """
    parameters = fittingParameters(concentrations=concentrations)
    
    fitParameters = pd.DataFrame(index=['lowerbound', 'initial', 'upperbound'],
                                 columns=['fmax', 'dG', 'fmin'])
    
    # find fmin
    param = 'fmin'
    fitParameters.loc[:, param] = [0, 0, np.inf]

    # find fmax
    param = 'fmax'
    fitParameters.loc[:, param] = [0, np.nan, np.inf]
    
    # find dG
    fitParameters.loc[:, 'dG'] = [parameters.find_dG_from_Kd(
        parameters.find_Kd_from_frac_bound_concentration(frac_bound, concentration))
             for frac_bound, concentration in itertools.izip(
                                    [0.99, 0.5, 0.01],
                                    [concentrations[0], concentrations[-1], concentrations[-1]])]
 
    return fitParameters

def getInitialFitParametersVary(concentrations):
    """ Return initial fit parameters from single cluster fits.
    
    Add a row 'vary' that indicates whether parameter should vary or not. """
    fitParameters = getInitialFitParameters(concentrations)
    return pd.concat([fitParameters.astype(object),
                      pd.DataFrame(True, columns=fitParameters.columns,
                                         index=['vary'], dtype=bool)])

    
def convertFitParametersToParams(fitParameters):
    """ Return lmfit params structure starting with descriptive dataframe. """
    param_names = fitParameters.columns.tolist()
    # store fit parameters in class for fitting
    params = Parameters()
    for param in param_names:
        if 'vary' in fitParameters.loc[:, param].index.tolist():
            vary = fitParameters.loc['vary', param]
        else:
            vary = True
        if 'lowerbound' in fitParameters.loc[:, param].index.tolist():
            lowerbound = fitParameters.loc['lowerbound', param]
        else:
            lowerbound = None
        if 'upperbound' in fitParameters.loc[:, param].index.tolist():
            upperbound = fitParameters.loc['upperbound', param]
        else:
            upperbound = np.inf
        params.add(param, value=fitParameters.loc['initial', param],
                   min = lowerbound,
                   max = upperbound,
                   vary= vary)
    return params

def getWeightsFromError(errors):
    """Given errors=[eminus, eplus], find weights based on inverse."""
    eminus, eplus = errors
    weights = 1/(eminus+eplus)
    if np.isnan(weights).any():
        weights = None
    return weights

def getWeightsFromBindingSeries(ys, num_cutoff=5):
    """Given binding series, find errors and then find weights."""
    if len(ys) >= num_cutoff:
        eminus, eplus = findErrorBarsBindingCurve(ys)
        weights = getWeightsFromError([eminus, eplus])
    else:
        weights = None
    return weights
        

def fitSingleCurve(x, y, fitParameters, func,
                          weights=None, do_not_fit=False, kwargs={}, min_kws={'maxfev':100}):
    """ Fit an objective function to data, weighted by errors. """

    # fit parameters
    params = convertFitParametersToParams(fitParameters)
    param_names = fitParameters.columns.tolist()
    
    # initiate output structure  
    final_params = pd.Series(index=(param_names +
                                    ['%s_stde'%param for param in param_names] +
                                    ['rsq', 'exit_flag', 'rmse']))
    
    # make sure fluorescence doesn't have NaN terms
    index = np.array(np.isfinite(y))
    
    # find the number of free parameters
    num_free_parameters = len(param_names)
    if 'vary' in fitParameters:
        num_free_parameters = fitParameters['vary'].sum()
    
    # don't fit if there are not enough finite entries    
    if index.sum() <= num_free_parameters:
        do_not_fit = True

    # return here if you don't want to actually fit
    if do_not_fit:
        final_params.loc['exit_flag'] = -1
        return final_params
    
    # add the arguments to the kwargs dict
    kwargs = kwargs.copy()
    kwargs.update({'data':y, 'weights':weights, 'index':index}) 

    # do the fit
    results = minimize(func, params,
                       args=(x,),
                       kws=kwargs, **min_kws)

    
    # find rqs
    rsq, rmse = get_rsq_rmse(y, results.residual)
    
    # save params in structure
    for param in param_names:
        final_params.loc[param] = params[param].value
        final_params.loc['%s_stde'%param] = params[param].stderr
    final_params.loc['rsq'] = rsq
    final_params.loc['exit_flag'] = results.ier
    final_params.loc['rmse'] = rmse
    
    return final_params

def get_rsq_rmse(y, residuals):
    """Return the rsq and the rmse given the y values and the residuals """
    ss_total = np.sum((y - y.mean())**2)
    ss_error = np.sum(residuals**2)
    rsq = 1-ss_error/ss_total
    rmse = np.sqrt(ss_error)
    
    return rsq, rmse

def findErrorBarsBindingCurve(subSeries, min_error=0):
    """ Return bootstrapped confidence intervals on columns of an input data matrix.
    
    Assuming rows represent replicate measurments, i.e. clusters. """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        eminus=[]
        eplus = [] 
        for i in subSeries:
            vec = subSeries.loc[:, i].dropna()
            success = True
            if len(vec) > 1:
                try:
                    bounds = bootstrap.ci(vec, np.median, n_samples=1000)
                except IndexError:
                    success = False
            else:
                success = False
                
            if success:
                eminus.append(vec.median() - bounds[0])
                eplus.append(bounds[1] - vec.median())
            else:
                eminus.append(np.nan)
                eplus.append(np.nan)
                
        
        eminus = pd.Series([max(min_error, e) for e in eminus], index=subSeries.columns)
        eplus = pd.Series([max(min_error, e) for e in eplus], index=subSeries.columns)
    return eminus, eplus

def enforceFmaxDistribution(median_fluorescence, fmaxDist, verbose=None, cutoff=None):
    """ Decide whether to enforce fmax distribution (on binding curves) or let it float.
    
    Cutoff is whether the last point of the (median) fluorescence is above the
    lower bound for fmax. """
    
    if verbose is None:
        verbose = False
    
    if cutoff is None:
        cutoff = 0.01 # only 2.5% of distribution falls beneath this value
    
    lowerbound = fmaxDist.ppf(cutoff)
    
    median_fluorescence = median_fluorescence.astype(float).values
    if median_fluorescence[-1] < lowerbound:
        redoFitFmax = True
        if verbose:
            print (('last concentration is below lb for fmax (%4.2f out of '
                   '%4.2f (%d%%). Doing bootstrapped fit with fmax'
                   'samples from dist')
                %(median_fluorescence[-1], lowerbound,
                  median_fluorescence[-1]*100/lowerbound))
    else:
        redoFitFmax = False
        if verbose:
            print (('last concentration is above lb for fmax (%4.2f out of %4.2f '+
                   '(%d%%). Proceeding by varying fmax')
                %(median_fluorescence[-1], lowerbound,
                  median_fluorescence[-1]*100/lowerbound))
    return redoFitFmax

def getClusterIndices(subSeries, n_samples=100, enforce_fmax=False, verbose=False):
    """Based on indices, find a list of sets of lcusters to do bootstrapping on."""
    # find number of samples to bootstrap
    numTests = len(subSeries)
    if numTests <10 and np.power(numTests, numTests) <= n_samples and not enforce_fmax:
        # then do all possible permutations
        if verbose:
            print ('Doing all possible %d product of indices'
                   %np.power(numTests, numTests))
        indices = [list(i) for i in itertools.product(*[subSeries.index]*numTests)]
    else:
        # do at most 'n_samples' number of iterations
        if verbose:
            print ('making %4.0f randomly selected (with replacement) '
                   'bootstrapped median binding curves')%n_samples
        indices = np.random.choice(subSeries.index,
                                   size=(n_samples, len(subSeries)), replace=True)
    return indices


def bootstrapCurves(x, subSeries, fitParameters, fmaxDist, func,
                    weighted_fit=True, verbose=False, n_samples=100,
                    enforce_fmax=None, min_error=0, func_kwargs={}):
    """ Bootstrap fit of a model to multiple measurements of a single molecular variant. """

    # if last point in binding series is below fmax constraints, do by method B
    median_fluorescence = subSeries.median()
        
    if enforce_fmax is None:
        # if enforce_fmax is not set, decide based on median fluorescence in last binding point
        enforce_fmax = enforceFmaxDistribution(median_fluorescence,
                                               fmaxDist, verbose=verbose)

            
    if enforce_fmax and (fmaxDist is None):
        print ('Error: if you wish to enforce fmax, need to define "fmaxDist"\n'
               'which is a instance of a normal distribution with mean and sigma\n'
               'defining the expected distribution of fmax')
        sys.exit()
        
    # estimate weights to use in weighted least squares fitting
    eminus = eplus = np.ones(len(x))*np.nan
    if weighted_fit:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                eminus, eplus = findErrorBarsBindingCurve(subSeries, min_error)
        except:
            pass
        
    # find number of samples to bootstrap
    indices = getClusterIndices(subSeries, n_samples, enforce_fmax)

    # Enforce fmax if initially told to and cutoff was not met
    fitParameters = fitParameters.copy()
    if enforce_fmax:
        # make sure fmax does not vary and find random variates
        # of fmax distribution
        if 'vary' not in fitParameters.index:
            fitParameters.loc['vary'] = True
        fitParameters.loc['vary', 'fmax'] = False
        fmaxes = fmaxDist.rvs(n_samples)
    
    # proceed with bootstrapping
    singles = {}
    for i, clusters in enumerate(indices):
        if verbose:
            if i%(n_samples/10.)==0:
                print 'working on %d out of %d, %d%%'%(i, n_samples, i/float(n_samples)*100)
        if enforce_fmax:
            fitParameters.loc['initial', 'fmax'] = fmaxes[i]
        
        # find median fluorescence
        fluorescence = subSeries.loc[clusters].median()
        
        # only fit if at least 3 measurements
        index = np.isfinite(fluorescence)
        if index.sum() <= 3:
            do_not_fit = True
        else:
            do_not_fit = False
        
        # do an iteration of fitting
        singles[i] = fitSingleCurve(x[index.values],
                                    fluorescence.loc[index],
                                    fitParameters,
                                    errors=[eminus[index.values], eplus[index.values]],
                                    func=func,
                                    do_not_fit=do_not_fit, kwargs=func_kwargs)
    # concatenate all resulting iterations
    singles = pd.concat(singles, axis=1).transpose()
    results = findProcessedSingles(singles, param_names)
    results.loc['rsq'] = findRsq(concentrations, median_fluorescence, convertFitParametersToParams(fitParameters), func)
    results.loc['numIter'] = (singles.exit_flag > 0).sum()
    results.loc['flag'] = enforce_fmax
    return results, singles

def findRsq(x, y, params, func):
    """find rsq of a fit."""
    ss_total = np.sum((y - y.mean())**2)
    ss_error = np.sum((y - func(params, x))**2)
    return 1-ss_error/ss_total


def findProcessedSingles(singles, param_names):
    """Given the output of the singles, find the upper and lower bounds."""
    # save results
    data = pd.concat({param:singles.loc[:, param].quantile([0.5, 0.025, 0.975]) for param in param_names}, axis=1)
    data.index = ['', '_lb', '_ub']
    results = data.stack().swaplevel(0,1).sort_index()
    results.index = [''.join(s) for s in results.index.tolist()]
    return results


def fitSetClusters(fitParams, ySeries, print_bool=True):
    """ Fit a set of curves. """
    singles = {}
    times = {}
    for i, idx in enumerate(ySeries.index.tolist()):
        # track progress
        if print_bool:
            num_steps = max(min(100, (int(len(ySeries)/100.))), 1)
            if (i+1)%num_steps == 0:
                print ('working on %d out of %d iterations (%d%%)'
                       %(i+1, len(ySeries), 100*(i+1)/
                         float(len(ySeries))))
                sys.stdout.flush()
        # fit single cluster
        y = ySeries.loc[idx]
        
        # find init time
        t0 = datetime.datetime.now()
        
        # fit
        singles[idx] = perCluster(fitParams, y)
        
        # find time diff
        t1 = datetime.datetime.now()
        times[idx] = (t1 - t0).total_seconds()
    return pd.concat(singles).unstack()

def perCluster(fitParams, y, plot=False):
    """ Fit a single binding curve. """
    fitParams.fit_curve(y)
    if plot:
        fitParams.plot_fit()
    return fitParams.results

def perVariant(variantParams, variant, n_samples=100, enforce_fmax=None, weighted_fit=True, min_error=0, func_kwargs={}):
    """ Fit a variant to objective function by bootstrapping median fluorescence. """

    t0 = datetime.datetime.now()
    variantParams.fit_set_binding_curves(variant, n_samples=n_samples, weighted_fit=weighted_fit,
                                         enforce_fmax=enforce_fmax)
    t1 = datetime.datetime.now()
    time_diff = (t1 - t0).total_seconds()
    return variantParams.results

def fitSetVariants(variantParams, variants=None,  n_samples=100, enforce_fmax=None, weighted_fit=True, min_error=0, func_kwargs={}, print_bool=True, return_time=False):
    """Fit a set of variants to objective function by bootstrapping median fluorescence."""
    t0 = datetime.datetime.now()
    results = variantParams.fit_binding_curves_all(variants, n_samples=n_samples,
                                                   weighted_fit=weighted_fit,
                                                   enforce_fmax=enforce_fmax,
                                                    print_bool=print_bool,
                                         return_results=True)
    t1 = datetime.datetime.now()
    time_diff = (t1 - t0).total_seconds()
    if return_time:
        return results, time_diff
    return results


def plotFitDistributions(results, singles, fitParameters):
    """ Plot a distribtion of fit parameters. """
    for param in fitParameters.columns.tolist():
    
        plt.figure(figsize=(4,3))
        sns.distplot(singles.loc[:, param].dropna().values,
                     hist_kws={'histtype':'stepfilled'}, color='b')
        ax = plt.gca()
        ylim = ax.get_ylim()
        plt.plot([results.loc[param]]*2, ylim, 'k--', alpha=0.5)
        plt.plot([results.loc['%s_lb'%param]]*2, ylim, 'k:', alpha=0.5)
        plt.plot([results.loc['%s_ub'%param]]*2, ylim, 'k:', alpha=0.5)
        plt.ylabel('prob density')
        plt.xlabel(param)
        plt.tight_layout()
    return

def returnParamsFromResults(final_params, param_names=None):
    """ Given results, convert to lmfit params structure for fitting. """
    if param_names is None:
        param_names = ['fmax', 'dG', 'fmin']
    params = Parameters()
    for param in param_names:
        params.add(param, value=final_params.loc[param])
    return params

def returnParamsFromResultsBounds(final_params, param_names, ub_vec):
    params_ub = Parameters()
    for param in ['%s%s'%(param, suffix) for param, suffix in
                  itertools.izip(param_names, ub_vec)]:
        name = param.split('_')[0]
        params_ub.add(name, value=final_params.loc[param])
    return params_ub

def errorPropagationKdFromKoffKobs(koff, kobs, c, sigma_koff, sigma_kobs):
    koff = koff.astype(float)
    kobs = kobs.astype(float)
    sigma_koff = sigma_koff.astype(float)
    sigma_kobs = sigma_kobs.astype(float)
    sigma_kd = np.sqrt(
        (( c*kobs/(kobs-koff)**2)*sigma_koff)**2 +
        ((-c*koff/(kobs-koff)**2)*sigma_kobs)**2)
    return sigma_kd

def errorProgagationKdFromdG(dG, sigma_dG):
    dG = dG.astype(float)
    sigma_dG = sigma_dG.astype(float)
    parameters = fittingParameters()
    sigma_kd = parameters.find_Kd_from_dG(dG)/parameters.RT*sigma_dG
    return sigma_kd

def returnResultsFromParams(params, results, y):
    """Given a set of params (i.e. from lmfit) return parsed."""
        # find rqs
    ss_total = np.nansum((y - y.mean())**2)
    ss_error = np.nansum((results.residual)**2)
    rsq = 1-ss_error/ss_total
    rmse = np.sqrt(ss_error)
    
    param_names = params.keys()
    index = (param_names + ['%s_stde'%param for param in param_names] +
             ['rsq', 'exit_flag', 'rmse'])
    final_params = pd.Series(index=index)

    # save params in structure
    for param in param_names:
        final_params.loc[param] = params[param].value
        final_params.loc['%s_stde'%param] = params[param].stderr
    final_params.loc['rsq'] = rsq
    final_params.loc['exit_flag'] = results.ier
    final_params.loc['rmse'] = rmse
    return final_params
