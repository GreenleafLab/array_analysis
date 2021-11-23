from lmfit import minimize, Parameters, report_fit, conf_interval
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scikits.bootstrap import bootstrap
from statsmodels.distributions.empirical_distribution import ECDF
import warnings
import itertools
import scipy.stats as st
import fittinglibs.plotting as plotting

class fmaxDistAny():
    __module__ = 'fittinglibs.%s'%os.path.splitext(os.path.basename(__file__))[0]
    # for fitting stde of fmaxes
    def __init__(self, params=None):
        self.params = params    
        
    def sigma_by_n_fit(self, params, x, y=None, weights=None):
        # objective for 1/sqrt(n) fit to stand error
        parvals = params.valuesdict()
        sigma = parvals['sigma']
        c   = parvals['c']
        fit = sigma/np.sqrt(x) + c
        if y is None:
            return fit
        elif weights is None:
            return y-fit
        else:
            return (y-fit)*weights
    
    def getDist(self, n, do_gamma=None):

        if self.params is None:
            print 'Error: define popts'
            return
        params = self.params
        
        sigma = self.sigma_by_n_fit(params, n)
        mean = params.valuesdict()['median']
        
        return self.find_fmax_bounds(mean, sigma,
                                     alpha=None,
                                     return_dist=True,
                                     do_gamma=do_gamma)

    def find_fmax_bounds(self, mean, sigma, alpha=None, return_dist=None, do_gamma=None):
        if alpha is None: alpha = 0.99
        if return_dist is None: return_dist = False
        if do_gamma is None:
            do_gamma = True
            
        self.do_gamma = do_gamma
        if do_gamma:    
            return self._find_fmax_bounds_given_gamma(mean, sigma, alpha, return_dist)
        else:
            return self._find_fmax_bounds_given_norm(mean, sigma, alpha, return_dist)
   
    def _get_percentiles_given_alpha(self, alpha):
        # first, middle, and last percentiles that cover fraction (alpha) of data 
        return (1-alpha)/2., 0.5, (1+alpha)/2
    
    def _find_fmax_bounds_given_gamma(self, mean, sigma, alpha, return_dist):
        # distribution is a gamma distribution
        k, theta = returnGammaParams(mean, sigma)
        dist = st.gamma(k, scale=theta)
        return self._return_bounds(dist, alpha, return_dist)

    def _find_fmax_bounds_given_norm(self, mean, sigma, alpha, return_dist):
        # distribution is a normal distribution
        dist = st.norm(loc=mean, scale=sigma)
        return self._return_bounds(dist, alpha, return_dist)
        
    def _return_bounds(self, dist, alpha, return_dist):
        # return upper and lower bounds or distribution
        if return_dist:
            return dist
        else:
            percentiles = self._get_percentiles_given_alpha(alpha)
            return dist.ppf(percentiles)

def getFractionOfData(n_test_counts, fraction_of_data):
    """Get a subset of the date? Don't actually remember."""
    for n in n_test_counts.index:
        if n_test_counts.loc[:n].sum()/float(n_test_counts.sum()) < fraction_of_data:
            pass
        else:
            min_n = n
            break
    return min_n

def fitGammaDistribution(vec, plot=None, set_mean=None, set_offset=None, initial_mean=None):
    """ Fit the CDF of a vector to the gamma distribution. """
    
    if plot is None:
        plot = False
    
    if initial_mean is None:
        # set initial mean for fit
        initial_mean = vec.median() 
    
    # fit the cdf of vec to find parameters that best represent the distribution
    ecdf = ECDF(vec)
    x, y = ecdf.x, ecdf.y
    param_names = ['mean', 'std', 'offset']
    params = Parameters()
    
    # distribution has mean (first moment)
    if set_mean is None:
        params.add('mean', value=initial_mean, vary=True, min=0, max=np.inf)
    else:
        params.add('mean', value=set_mean, vary=False)

    # distribution has standard deviation (sqrt of second moment)
    params.add('std',  value=vec.std(), vary=True, min=0, max=np.inf)

    # gamma distribution can have offset (probably centered at 0)
    if set_offset is None:
        params.add('offset', value=0, vary=True)
    else:
        params.add('offset', value=set_offset, vary=False)
        
    # go through with least squares minimization
    results = minimize(gammaObjective, params,
                       args=(x,),
                       kws={'data':y})

    # save params and rsq
    index = (param_names + ['%s_stde'%param for param in param_names] +
             ['rsq', 'exit_flag', 'rmse'])
    final_params = pd.Series(index=index)
    
    # save params in structure
    for param in param_names:
        final_params.loc[param] = params[param].value
        final_params.loc['%s_stde'%param] = params[param].stderr

    # find rqs
    ss_total = np.sum((y - y.mean())**2)
    ss_error = np.sum((results.residual)**2)
    rsq = 1-ss_error/ss_total
    rmse = np.sqrt(ss_error)
    final_params.loc['rsq'] = rsq
    final_params.loc['exit_flag'] = results.ier
    final_params.loc['rmse'] = rmse
    
    if plot:
        # plot pdf 
        plotting.plotGammaFunction(vec, gammaObjective, params=params)
    return final_params

def gammaObjective(params, x, data=None, weights=None, return_pdf=None):
    """Objective function for a gamma cdf"""
    if return_pdf is None:
        return_pdf is False
    parvals = params.valuesdict()
    mean = parvals['mean']
    std  = parvals['std']
    mu = parvals['offset']
    
    k, theta = returnGammaParams(mean, std)

    if return_pdf:
        cdf = st.gamma.pdf(x, k, scale=theta, loc=mu)
    else:
        cdf = st.gamma.cdf(x, k, scale=theta, loc=mu)
    
    if data is None:
        return cdf
    elif weights is None:
        return cdf - data
    else:
        return (cdf - data)*weights
    

def returnGammaParams(mean, std):
    """Return shape and scale parameters from moments of distrubtion mean, std"""
    k = (mean/std)**2
    theta = (std)**2/mean
    return k, theta
        
def fitSigmaDist(x, y, weights=None, set_c=None, at_n=None):
    """Fit the relationship between the stde of a distribution and the number of measurements."""
    # fit how sigmas scale with number of measurements
    fmaxDist = fmaxDistAny()
    params = Parameters()
    params.add('sigma', value=y.max(), min=0)
    if set_c is None:
        params.add('c',     value=y.min(),   min=0)
    else:
        params.add('c', expr="%4.5f-sigma/sqrt(%4.5f)"%(set_c, at_n))
    minimize(fmaxDist.sigma_by_n_fit, params,
                                   args=(x,),
                                   kws={'y':y,
                                        'weights':weights},
                                   xtol=1E-6, ftol=1E-6, maxfev=10000)
    return params

def findMinStd(fmaxes, n_tests, mean_fmax, fraction_of_data=0.5):
    """ Find a minimum standard deviation. """
    n_test_counts = n_tests.value_counts().sort_index()
    # find number of tests that contains 50% of the variants
    min_n = getFractionOfData(n_test_counts, fraction_of_data)
    
    # fit to find the "floor" of std below which represents experimental noise
    fmaxHighNFit = fitGammaDistribution(fmaxes.loc[n_tests>=min_n],
                                        set_mean=mean_fmax, set_offset=0)
    min_std = fmaxHighNFit.loc['std']
    at_n = np.average(n_test_counts.loc[min_n:].index.tolist(),
                      weights=n_test_counts.loc[min_n:].values)
    return min_std, at_n

def findStdParams(fmaxes, n_tests, mean_fmax, min_std, at_n, min_num_to_fit=4, single_std=False):
    """ Find the relationship between number of tests and std. """
    if single_std:
        # don't fit. return min_std if not None, else fit all points.
        stds_actual = pd.concat({1:fitGammaDistribution(fmaxes, set_mean=mean_fmax)}).unstack()
        if min_std is None:
            min_std = stds_actual.loc[1, 'std']
        params = Parameters()
        eps = min_std/1E6 # a small value
        params.add('sigma', value=eps)
        params.add('c', value=min_std)
        params.add('median', value=mean_fmax)
        return params, stds_actual

    n_test_counts = n_tests.value_counts().sort_index()
    all_ns = n_test_counts.loc[n_test_counts>=min_num_to_fit].index.tolist()
    
    stds_actual = {}
    for n in all_ns:
        stds_actual[n] = fitGammaDistribution(fmaxes.loc[n_tests==n],
                                              set_mean=mean_fmax,)
    stds_actual = pd.concat(stds_actual).unstack()
    stds_actual.dropna(inplace=True)

    x = stds_actual.index.tolist()
    y = stds_actual.loc[:, 'std']
    weights_fit = np.sqrt(n_test_counts.loc[stds_actual.index])

    params = fitSigmaDist(x, y,
                          weights=weights_fit,
                          set_c=min_std, at_n=at_n)
    params.add('median', value=mean_fmax)
    return params, stds_actual
    
def findParams(tight_binders, use_simulated=None, table=None, single_std=False):
    """ Initialize, find the fmax distribution object and plot. """
    if use_simulated is None:
        use_simulated = False
    if use_simulated and (table is None):
        print "error: need to give table of all cluster fits to do simulated data"
        return
    
    mean_fmax, bounds, loose_bounds = getFmaxMeanAndBounds(tight_binders)
    fmaxes_data, n_tests_data = getFmaxesToFit(tight_binders, bounds=bounds)
    if len(fmaxes_data) < 3:
        print "error: not enough fmaxes to fit distribution"
        return
    
    if use_simulated:
        # find min std of distribution anyways
        min_std, at_n = findMinStd(fmaxes_data, n_tests_data, mean_fmax)
        # if at_n is None then we can skip fitting.
        fmaxes, n_tests = getFmaxesToFitSimulated(table, tight_binders.index, bounds=bounds)
        if np.around(at_n)==1: single_std = True
    else:
        min_std = None; at_n = None
        fmaxes, n_tests = fmaxes_data, n_tests_data

    # fit relationship of std with number of measurements
    params, stds_actual = findStdParams(fmaxes, n_tests, mean_fmax, min_std, at_n, single_std=single_std)
    fmaxDist = fmaxDistAny(params=params)

    # plot
    maxn = getFractionOfData(n_tests_data.value_counts().sort_index(), 0.995)
    ax = plotting.plotFmaxStdeVersusN(fmaxDist, stds_actual, maxn)
    # if use_simulated, plot actual as well
    if use_simulated:
        try:
            params1, stds_not_sim = findStdParams(fmaxes_data, n_tests_data, mean_fmax, None, None)
            fmaxDist1 = fmaxDistAny(params=params1)
            plotting.plotFmaxStdeVersusN(fmaxDist1, stds_not_sim, maxn, ax=ax)
        except:
            pass
        
    # plot offsets
    plotting.plotFmaxOffsetVersusN(fmaxDist, stds_actual, maxn)
    # if use_simulated, plot actual as well
    if use_simulated:
        try:
            params1, stds_not_sim = findStdParams(fmaxes_data, n_tests_data, mean_fmax, None, None)
            fmaxDist1 = fmaxDistAny(params=params1)
            plotting.plotFmaxOffsetVersusN(fmaxDist1, stds_not_sim, maxn, ax=ax)
        except:
            pass  

    # plot numbers
    plotting.plotNumberVersusN(n_tests_data, maxn)
    
    return fmaxDist

def resultsFromFmaxDist(fmaxDist, n):
    """ Return results from fmaxDist. """
    mean, var = fmaxDist.getDist(n).stats(moments='mv')
    return pd.Series({'std':np.sqrt(var), 'mean':mean, 'offset':0})

def findUpperboundFromFmaxDistObject(fmaxDist):
    """ Return loose bounds on fmax from fmaxDist. """
    return fmaxDist.params['median'] + 3*fmaxDist.params['sigma']


def returnFminFromFits(variant_table, cutoff):
    """ Return the estimated fixed fmin based on affinity and fits. """
    return variant_table.loc[variant_table.dG_init> cutoff].fmin_init.median()

def findInitialPoints(variant_table, param_names=['fmax', 'dG', 'fmin']):
    """ Return initial points with different column names. """
    initialPoints = variant_table.loc[:, ['%s_init'%param_name for param_name in param_names] + ['numTests']]
    initialPoints.columns = param_names + ['numTests']  
    return initialPoints

def returnFminFromFluorescence(initialPoints, fluorescenceMat, cutoff):
    """ Return the estimated fixed fmin based on affinity and fluroescence. """
    # if cutoff is not given, use parameters
    initial_dG = initialPoints.loc[:, 'dG']
    firstBindingPoint = getMedianFirstBindingPoint(fluorescenceMat)
    return firstBindingPoint.loc[initial_dG.index].loc[initial_dG > cutoff].median()

def getMedianFirstBindingPoint(table):
    """ Return the median fluoresence in first binding point of each variant. """
    return table.groupby('variant_number').median().iloc[:, 0]

def getFmaxMeanAndBounds(tight_binders, cutoff=1E-12):
    """ Return median fmaxes of variants. """
    # find defined mean shared by all variants by fitting all
    fmaxes = tight_binders.fmax
    fmaxAllFit = fitGammaDistribution(fmaxes, set_offset=0, initial_mean=fmaxes.median())
    
    # use fit to also define upper and lower bound of expected values
    mean_fmax = fmaxAllFit.loc['mean']
    std_fmax =  fmaxAllFit.loc['std']
    lowerbound, median, upperbound = (fmaxDistAny().
                                      find_fmax_bounds(mean_fmax, std_fmax, alpha=1-cutoff))
    loose_bounds = [0, 2*upperbound]
    plotting.plotGammaFunction(fmaxes, gammaObjective, fmaxAllFit, bounds=loose_bounds)
    plt.axvline(lowerbound, color='k', linestyle=':')
    plt.axvline(upperbound, color='k', linestyle=':')
    return mean_fmax, [lowerbound, upperbound], loose_bounds
    

def getFmaxesToFit(tight_binders, bounds=[0, np.inf]):
    """ Return fmax initial fits that fall within bounds. """
    fmaxes = tight_binders.fmax
    
    # find those within bounds
    index = (fmaxes>=bounds[0])&(fmaxes<=bounds[1])
    n_tests = tight_binders.loc[index].numTests
    return fmaxes.loc[index], n_tests

def getFmaxesToFitSimulated(all_clusters, good_variants, bounds=[0, np.inf], n_subset=np.arange(1,15)):
    """ Return simulated median fmaxes from randomly sampled clusters.
    
    Returns a distribution of fmaxes for increasing number of samples.

    Parameters:
    -----------
    all_clusters : per-cluster DataFrame giving initial fits. Columns variant
        number and fmax.
    good_variants : list of variant numbers that represent variants passing cutoffs.
    bounds : only include median values of fmax within bounds.
    n_subset : list of number of tests to simulate.
    
    Returns:
    --------
    DataFrame of median fmaxes
    DataFrame of number of tests
    
    """
        
    # find those clusters associated with good variants and use those it fmaxes
    good_clusters = pd.Series(np.in1d(all_clusters.variant_number,
                                      good_variants),
                              index=all_clusters.index)
    fmax_clusters = all_clusters.loc[good_clusters, 'fmax']
    
    # first assign each cluster to a variant. different assignments for each n
    fmax_df = pd.concat([fmax_clusters,
                         pd.DataFrame(index=fmax_clusters.index, columns=n_subset)], axis=1)
    for n in n_subset:
        # make at most 1000 or max_num_variants independent subsamples of n samples each 
        num_variants = min(1000, np.floor(len(fmax_clusters)/float(n)))
        n_samples = int(n*num_variants)
        clusters = np.random.choice(fmax_clusters.index, n_samples, replace=False)
        fmax_df.loc[clusters, n] = np.tile(np.arange(num_variants), n)
    
    # now find median of fmaxes assigned to each variant
    fmaxes = []
    n_tests = []
    for n in n_subset:
        
        # only take values within bounds
        vec = fmax_df.groupby(n)['fmax'].median()
        index = (vec>=bounds[0])&(vec<=bounds[1])
        
        # make vector of tests
        n_vec = vec.loc[index].copy()
        n_vec.loc[:] = n
        
        # save
        fmaxes.append(vec.loc[index])
        n_tests.append(n_vec)
    
    return pd.concat(fmaxes, ignore_index=True), pd.concat(n_tests, ignore_index=True)
