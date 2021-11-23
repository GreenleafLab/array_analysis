from lmfit import minimize, Parameters, report_fit, conf_interval
import numpy as np
import pandas as pd
import sys
import warnings
import itertools
import scipy.stats as st
import copy
#from sklearn import metrics
from fittinglibs import objfunctions, variables, fitting, plotting

class FitParams():
    """Class with attributes objective function, initial params, upper/lowerbounds on params."""
    def __init__(self, func_name, x=None, y=None, init_kws={}, fit_kws={}, before_fit_ops=[]):
        self.func_name = func_name
        self.func = getattr(objfunctions, func_name)
        self.param_names = self.func(None, None, return_param_names=True)
        self.x = x
        self.y = y
        self.fit_kws = fit_kws  # kwargs pairs taat will be passed to fittings.fitSingleCurve
        self.init_kws = init_kws # dict of {param_name:{'initial':val}} etc

        # these operations will be performed before fitting
        self.before_fit_ops = before_fit_ops # formatted like [('fmax', 'initial', lambda y: y.max)]
        
        if x is not None:        
            # go through and get default initiation point for each fit param
            fit_parameters = {}
            for param_name in self.param_names:
                fit_parameters[param_name] = getattr(self, 'get_init_%s'%param_name)()
                # TODO add functions for koff etc.
            self.fit_parameters = fit_parameters
            self.update_init_params(**init_kws)
        
        
  
    def update_init_params(self, **init_kws):
        """Given kwargs, update the fit_parameter values."""
        self.fit_parameters = _update_init_params(self.fit_parameters, **init_kws)
    

    
        
    def get_init_params(self, y=None, fit_parameters=None):
        """Return lmfit Parameters class."""
        if fit_parameters is not None:
            return _get_init_params(fit_parameters)
        
        fit_parameters = self.fit_parameters
        # TODO: code below is repeated in fit_curve...way to better interface?
        if self.before_fit_ops:
            for (param_name, key, operation) in self.before_fit_ops:
                fit_parameters = _update_init_params(fit_parameters, **{param_name:{key:operation(y)}})
        return _get_init_params(fit_parameters)
    
    def fit_curve(self, y, fit_parameters=None, weights=None, min_kws={'xtol':1E-6, 'ftol':1E-6, 'maxfev':100}, return_results=False):
        """fit a single curve to y values."""
        if fit_parameters is None:
            fit_parameters = self.fit_parameters
        
        # change fit parameters according to list of operations 
        if self.before_fit_ops:
            for (param_name, key, operation) in self.before_fit_ops:
                if fit_parameters[param_name]['vary']:
                    fit_parameters = _update_init_params(fit_parameters, **{param_name:{key:operation(y)}})
        
        # weight fit if weights are given
        results = fitting.fitSingleCurve(self.x, y, _convert_to_expected_fit_parameters(fit_parameters), self.func, weights=weights, kwargs=self.fit_kws, min_kws=min_kws)
        
        if return_results:
            return results
        
        self.y = y
        self.results = results

    def get_params_from_results(self, results=None):
        """Convert between params and results"""
        if results is None:
            results = getattr(self, 'results')
        return _get_params_from_results(results, self.param_names)
    
    def plot_fit(self, y=None, results=None, plot_fit=True, **kwargs):
        """Assuming fit_curve has already been run, plot results."""
        x = self.x
        if y is None:
            y = getattr(self, 'y')
        if results is None:
            try:
                results = getattr(self, 'results')
            except AttributeError:
                plot_fit = False
        
        
        # generate x values for fit function
        x = np.array(x)
        more_x = np.logspace(np.log10(x.min()/10), np.log10(x.max()*2), 100)
    
        # plot the data
        line, = plotting.plt.plot(x, y, 'o', **kwargs)
        if plot_fit:
            if 'color' in kwargs.keys(): kwargs.pop('color')
            params = _get_params_from_results(results, self.param_names)
            plotting.plt.plot(more_x, self.func(params, more_x), color=line.get_color(), **kwargs)


    def plot_initfit(self, y=None, fit_parameters=None, **kwargs):
        """Assuming fit_curve has already bee run, plot results."""
        x = self.x
        params = self.get_init_params(y, fit_parameters=fit_parameters)
        
        # generate x values for fit function
        x = np.array(x)
        more_x = np.logspace(np.log10(x.min()/10), np.log10(x.max()*2), 100)
    
        # plot the data
        plotting.plt.plot(more_x, self.func(params, more_x), '.--', **kwargs)
    
        
    def get_init_dG(self):
        """Given the initial concentrations, find an initial dG value."""
        parameters = variables.fittingParameters()
        num_x = len(self.x)
        lb   = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.99, self.x[0]))
        init = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.5,  self.x[-1]))
        ub   = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.01, self.x[-1]))
        return {'lowerbound':lb, 'initial':init, 'upperbound':ub, 'vary':True}
        
    def get_init_fmax(self):
        """Find initial fmax values."""
        return {'lowerbound':0, 'initial':np.nan, 'upperbound':np.inf, 'vary':True}

    def get_init_fmin(self):
        """Find initial fmin values."""
        return {'lowerbound':0, 'initial':1E-4, 'upperbound':np.inf, 'vary':True}

    def get_init_dGns(self):
        """Find initial nonspecific dG val values."""
        parameters = variables.fittingParameters()
        lb   = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.99, self.x[0]))
        init = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.1,  self.x[-1]))
        ub   = parameters.find_dG_from_Kd(parameters.find_Kd_from_frac_bound_concentration(0.0001, self.x[-1]))
        return {'lowerbound':lb, 'initial':init, 'upperbound':ub, 'vary':True}
        
class MoreFitParams():
    """Add more attributes to enable fitting of curves with specific distribution for fmax in particular."""
    def __init__(self, fitParams, initial_points=None, binding_series_dict=None, binding_series=None, annotated_clusters=None, fmax_dist_obj=None, results=None):

        self.fitParams = fitParams 
        self.initial_points = initial_points
        self.binding_series_dict = binding_series_dict
        self.binding_series = binding_series
        self.annotated_clusters = annotated_clusters
        self.fmax_dist_obj = fmax_dist_obj
        self.results_all = results
        #self.variants = 
        if initial_points is not None and binding_series_dict is not None:
            self.variants = list(set(binding_series_dict.groupby(level=0).first().index.tolist()).union(initial_points.index.tolist()))     
        # make sure there are values for initial points for all variants

    def get_ys(self, idx):
        """Find the set of y values associated with a variant."""
        if self.binding_series_dict is None:
            if self.binding_series is None or self.annotated_clusters is None:
                raise IndexError('Need to define either binding_series_dict or binding series and annotated clusters')
            index = self.annotated_clusters.loc[self.annotated_clusters.variant_number==idx].index
            return self.binding_series.loc[index]
        else:
            try:
                return self.binding_series_dict.loc[idx]
            except KeyError:
                return pd.DataFrame(columns=self.binding_series_dict.columns)
                

    def fit_set_binding_curves(self, idx, enforce_fmax=None, weighted_fit=True, n_samples=100, return_results=False, use_initial=False):
        """Fit a set of y values to a curve."""
        x = self.fitParams.x
        ys = self.get_ys(idx)        
        y = ys.median()
        num = len(ys)
        fmax_dist = self.fmax_dist_obj.getDist(num)
        fit_parameters = copy.deepcopy(self.fitParams.fit_parameters)
        param_names = self.fitParams.param_names
        initial_points = self.initial_points.loc[idx]
        
        # update parameters based on initial points
        if use_initial:
            init_kws = {}
            for param_name, param_dict in fit_parameters.items():
                new_init = initial_points.loc[param_name]
                if param_dict['vary'] and not pd.isnull(new_init):
                    init_kws[param_name] = {'initial':new_init}              
            fit_parameters = _update_init_params(fit_parameters, **init_kws)
        
        # decide whether to enforce fmax or not
        lb, ub = fmax_dist.interval(0.99)
        if enforce_fmax is None:
            enforce_fmax = (fitting.enforceFmaxDistribution(y, fmax_dist) )
            
        # if enforce fmax, find set of fmaxes to use
        if enforce_fmax:
            fmaxes = fmax_dist.rvs(n_samples)

        # find error on ys and use this to weight if weighted fit option is given.
        if weighted_fit:
            weights = fitting.getWeightsFromBindingSeries(ys) # note here is a cutoff
        else:
            weights = None

        # get bootstrap indices
        indices = fitting.getClusterIndices(ys, n_samples, enforce_fmax)
        
        # for each set of clusters, fit
        singles = []
        # if ys is empty, don't fit
        if len(ys) == 0 :
            if return_results:
                return None, singles, ys, indices
            else:
                return
            
        for i, index in enumerate(indices):
            if enforce_fmax:
                fit_parameters_fmax = _update_init_params(fit_parameters, fmax={'initial':fmaxes[i], 'vary':False})
            else:
                fit_parameters_fmax = _update_init_params(fit_parameters, fmax={'lowerbound':lb, 'upperbound':ub})
            #single = fitting.fitSingleCurve(x, ys.loc[index].median(), _convert_to_expected_fit_parameters(fit_parameters), func, weights=weights, kwargs=fit_kws, min_kws=min_kws)
            #singles.append(single)
            results = self.fitParams.fit_curve(ys.loc[index].median(), fit_parameters=fit_parameters_fmax, weights=weights, return_results=True)
            singles.append(results)
        singles = pd.concat(singles, axis=1).transpose()

        # find upper and lower bounds from results
        sub_param_names = [key  for key in param_names if fit_parameters[key]['vary']] # because you want to include ub and lb of fmax
        results = fitting.findProcessedSingles(singles, sub_param_names)
        for key in [param_name for param_name in param_names if param_name not in sub_param_names]:
            results.loc[key] = singles.loc[:, key].median()
        for key in param_names:
            results.loc['%s_init'%key] = initial_points.loc[key]
        
        # find rsq of median
        ypred = pd.Series(self.fitParams.func(_get_params_from_results(results, param_names), x), index=y.index)
        ypred_init = pd.Series(self.fitParams.func(_get_params_from_results(initial_points, param_names), x), index=y.index)

        results.loc['rsq']      = get_r2_score(y, ypred)
        results.loc['rsq_init'] = get_r2_score(y, ypred_init)
        results.loc['rmse']     = get_rmse(y, ypred)
        results.loc['rmse_init'] = get_rmse(y, ypred_init)
        results.loc['num_iter'] = (singles.exit_flag > 0).sum()
        results.loc['num_tests'] = len(ys)
        results.loc['fmax_enforced'] = enforce_fmax
        order = [s for s in results.index.tolist() if s.find('_init')>-1] + [s for s in results.index.tolist() if s.find('_init')==-1]
        results = results.loc[order]
        
        if return_results:
            return results, singles, ys, indices
        
        self.results = results
        self.ys = ys
    
    def fit_binding_curves_all(self, variants=None, enforce_fmax=None, weighted_fit=True, n_samples=100, print_bool=True, return_results=False):
        """Fit a set of variants to curves with bootstrapping method."""
        if variants is None:
            variants = self.variants
        results = {}
        for i, idx in enumerate(variants):

            # track progress
            if print_bool:
                num_steps = max(min(100, (int(len(variants)/100.))), 1)
                if (i+1)%num_steps == 0:
                    print ('working on %d out of %d iterations (%d%%)'
                           %(i+1, len(variants), 100*(i+1)/
                             float(len(variants))))
                    sys.stdout.flush() 
            
            # fit individual curves
            vec = self.fit_set_binding_curves(idx,
                                                       enforce_fmax=enforce_fmax,
                                                       weighted_fit=weighted_fit,
                                                       n_samples=n_samples,
                                                       return_results=True)[0]
            if vec is not None:
                # don't include None
                results[idx] = vec

        results = pd.concat(results).unstack()
        if return_results:
            return results
        
        else:
            self.results_all = results
    
    def plot_specific_binding_curve(self, variant, annotate=True, **kwargs):
        """Plot the fit of a particalr binding curve assuming reuslts_all is set."""
        results_all = getattr(self, 'results_all')
        ys = self.get_ys(variant)
        self.plot_binding_curves(ys, results_all.loc[variant], **kwargs)
        if annotate:
            self.annotate_curve(results=results_all.loc[variant])
    
    
    def plot_binding_curves(self, ys=None, results=None, **kwargs):
        """Plot the resulting fits."""

        if ys is None:
            ys = getattr(self, 'ys')
        if results is None:
            results = getattr(self, 'results')

        x = self.fitParams.x
        y = ys.median()
        func = self.fitParams.func
        param_names = self.fitParams.param_names
 
        params = _get_params_from_results(results, param_names)
        
        # generate x values for fit function
        x = np.array(x)
        more_x = np.logspace(np.log10(x.min()/10), np.log10(x.max()*2), 100)
    
        # plot the data
        line, = plotting.plt.plot(x, y, 'o', **kwargs)
        if 'color' in kwargs.keys(): kwargs.pop('color')
        
        # plot the errorbars
        eminus, eplus = fitting.findErrorBarsBindingCurve(ys)
        plotting.plt.errorbar(x, y, yerr=[eminus, eplus], fmt=',', color=line.get_color(), linewidth=0.5, **kwargs)
        
        # plot the fit
        plotting.plt.plot(more_x, func(params, more_x), color=line.get_color(), **kwargs)
        
        # plot the bounds
        if self.fitParams.func_name=='binding_curve' or self.fitParams.func_name=='binding_curve_nonlinear':
            param_names_ub = [param + '_lb' if param=='dG' else param + '_ub' if param=='fmax' else param for param in param_names]
            param_names_lb = [param + '_ub' if param=='dG' else param + '_lb' if param=='fmax' else param for param in param_names]
            if all([s in results.index.tolist() for s in param_names_ub + param_names_lb]):
                params_ub = _get_params_from_results(results, param_names_ub, param_names)
                params_lb = _get_params_from_results(results, param_names_lb, param_names)
                plotting.plt.fill_between(more_x, func(params_lb, more_x), func(params_ub, more_x), alpha=0.5, color='0.5')
    
    def plot_init_binding_curves(self, results=None, **kwargs):
        """Given the initial parameters, plot the resulting fit."""

        if results is None:
            results = getattr(self, 'results')

        x = self.fitParams.x
        func = self.fitParams.func
        param_names = self.fitParams.param_names
 
        params = _get_params_from_results(results, ['%s_init'%s for s in param_names], param_names)        
        more_x = np.logspace(np.log10(x.min()/10), np.log10(x.max()*2), 100)
        plotting.plt.plot(more_x, func(params, more_x), **kwargs)
        
    def annotate_curve(self, results=None, other_params=['rsq', 'num_tests'], other_param_formats={'rsq':'%4.3f', 'num_tests':'%d'}, **kwargs):
        """annotate the current plot with information about the fit."""
        if results is None:
            results = self.results
        param_names = self.fitParams.param_names
        annotate_list = []
        for param_name in param_names:
            s = '%s = %4.2f'%(param_name, results.loc[param_name])
            if param_name + '_lb' in results.index.tolist() and param_name + '_ub' in results.index.tolist():
                s = s + ' (%4.2f, %4.2f)'%(results.loc[param_name + '_lb'], results.loc[param_name + '_ub'])
            annotate_list.append(s)
        for idx in other_params:
            if idx in results.index.tolist():
                annotate_list.append(('%s = ' + other_param_formats[idx])%(idx, results.loc[idx]))
                
        plotting.annotate_axes('\n'.join(annotate_list))
     
        

def _convert_to_expected_fit_parameters(fit_parameters):
    """Assuming fit parameters are in dict format, convert to old format (dataframe)."""
    return pd.concat({key:pd.Series(val) for key, val in fit_parameters.items()}).unstack(level=0)

def _get_params_from_results(results, param_names, rename_param_names=None):
    """Convert between params and results."""
    if rename_param_names is None: rename_param_names = param_names
    params = Parameters()
    for param, rename_param in zip(param_names, rename_param_names):
        params.add(rename_param, value=results.loc[param])
    return params

def _get_init_params(fit_parameters):
    """Return lmfit Parameters class."""
    params = Parameters()
    for param_name, init_dict in fit_parameters.items():
        params.add(param_name, value=init_dict['initial'],
               min = init_dict['lowerbound'],
               max = init_dict['upperbound'],
               vary= init_dict['vary'])
    return params

def _update_init_params(fit_parameters, **init_kws):
    """Return an updated fit parameters dict bsed on new vals in init_kws."""
    fit_parameters_new = copy.deepcopy(fit_parameters)
    for param_name, param_dict in fit_parameters_new.items():
        if param_name in init_kws.keys():
            for key, val in init_kws[param_name].items():
                if val is not None:
                    param_dict[key] = val
        fit_parameters_new[param_name] = param_dict
    return fit_parameters_new
        
def process_new_params(args ):
    """Return processed values such that each input argument is a list of the same length."""
    args = copy.copy(args)
    names = ['params_init', 'params_vary', 'params_lb', 'params_ub']
    
    # if params_name is defined, make sure all attributes are a list of the same length
    if args.params_name is not None:
        for name in names:
            val = getattr(args, name)
            if val is None:
                setattr(args, name, [None]*len(args.params_name))
            elif len(val) != len(args.params_name):
                raise IndexError('must have same number of values in --params_name and --%s'%name)
    elif any([getattr(args, name) for name in names]):
        print 'Warning: no params are changed. Need to define --params_name if you want to effect changes in the fit params.'
    
    # if there is nothing to change, set all to be a list of length 1
    if args.params_name is None:
        for name in names + ['params_name']:
            setattr(args, name, [None])
    return args

def get_r2_score(y, ypred):
    """Find the rsq, leaving out NaNs."""
    index = pd.concat([y, ypred], axis=1).dropna().index.tolist()
    return st.pearsonr(y.loc[index], ypred.loc[index])[0]**2

def get_rmse(y, ypred):
    """Find the rmse."""
    return np.sqrt(((y-ypred)**2).sum()/len(y))
        