import numpy as np
import pandas as pd
from fittinglibs import fitting
from fittinglibs.variables import fittingParameters
from math import factorial

def rates_off(params, times, data=None, weights=None, index=None, bleach_fraction=1, image_ns=None, return_param_names=False):
    """ Return fit value, residuals, or weighted residuals of off rate objective function. """
    if return_param_names:
        return ['fmax', 'koff', 'fmin']
    
    # process inputs
    if index is None:
        # use all values
        index = np.ones(len(times)).astype(bool)
    if image_ns is None:
        # assume every image is the sequential image order
        image_ns = np.arange(len(times))
        
    parvals = params.valuesdict()
    fmax = parvals['fmax']
    koff = parvals['koff']
    fmin = parvals['fmin']
    fracbound = (fmin +
                 (fmax - fmin)*np.exp(-koff*times)*
                 np.power(bleach_fraction,image_ns))

    # return fit value of data is not given
    if data is None:
        return fracbound[index]
    
    # return residuals if data is given
    elif weights is None:
        return (fracbound - data)[index]
    
    # return weighted residuals if data is given
    else:
        return ((fracbound - data)*weights)[index]  
    
    
def rates_on(params, times, data=None, weights=None, index=None,  bleach_fraction=1, image_ns=None, return_param_names=False):
    """ Return fit value, residuals, or weighted residuals of on rate objective function. """
    if return_param_names:
        return ['fmax', 'kobs', 'fmin']    
    if index is None:
        index = np.ones(len(times)).astype(bool)
    if image_ns is None:
        image_ns = np.arange(len(times))
        
    parvals = params.valuesdict()
    fmax = parvals['fmax']
    kobs = parvals['kobs']
    fmin = parvals['fmin']
    fracbound = fmin + (fmax*(1 - np.exp(-kobs*times)*np.power(bleach_fraction,image_ns)));

    # return fit value of data is not given
    if data is None:
        return fracbound[index]
    
    # return residuals if data is given
    elif weights is None:
        return (fracbound - data)[index]
    
    # return weighted residuals if data is given
    else:
        return ((fracbound - data)*weights)[index]  
        
def binding_curve(params, concentrations, data=None, weights=None, index=None, return_param_names=False):
    """  Return fit value, residuals, or weighted residuals of a binding curve.
    
    Hill coefficient 1. """
    if return_param_names:
        return ['fmax', 'dG', 'fmin']
    
    if index is None:
        index = np.ones(len(concentrations)).astype(bool)
        
    parameters = fittingParameters()
    
    parvals = params.valuesdict()
    fmax = parvals['fmax']
    dG   = parvals['dG']
    fmin = parvals['fmin']

    fracbound = (fmin + fmax*concentrations/
                 (concentrations + np.exp(dG/parameters.RT)/
                  parameters.concentration_units))
    
    # return fit value of data is not given
    if data is None:
        return fracbound[index]
    
    # return residuals if data is given
    elif weights is None:
        return (fracbound - data)[index]
    
    # return weighted residuals if data is given
    else:
        return ((fracbound - data)*weights)[index]
    
def binding_curve_linear(params, concentrations, data=None, weights=None, index=None, return_param_names=False):
    """  Return fit value, residuals, or weighted residuals of a binding curve.
    
    Hill coefficient 1. """
    if return_param_names:
        return ['fmax', 'dG', 'fmin', 'slope']
    
    if index is None:
        index = np.ones(len(concentrations)).astype(bool)
        
    parameters = fittingParameters()
    
    parvals = params.valuesdict()
    fmax = parvals['fmax']
    dG   = parvals['dG']
    fmin = parvals['fmin']
    slope = parvals['slope']
    
    fracbound = (fmin + fmax*concentrations/
                 (concentrations + np.exp(dG/parameters.RT)/
                  parameters.concentration_units)) + slope*concentrations
    
    # return fit value of data is not given
    if data is None:
        return fracbound[index]
    
    # return residuals if data is given
    elif weights is None:
        return (fracbound - data)[index]
    
    # return weighted residuals if data is given
    else:
        return ((fracbound - data)*weights)[index]
    
def powerlaw(params, x, y=None, weights=None, index=None, return_param_names=False):
    """"""
    if return_param_names:
        return ['c', 'exponent', 'amplitude']
    if index is None: index = np.ones(len(x)).astype(bool)
    parvals = params.valuesdict()
    c = parvals['c']
    k = parvals['exponent']
    A = parvals['amplitude']

    y_pred = A*np.power(x, k) + c
    if y is None:
        return y_pred[index]
    elif weights is None:
        return (y - y_pred)[index]
    else:
        return ((y - y_pred)*weights)[index]
    
def exponential(params, x, y=None, weights=None, return_param_names=False):
    """"""
    if return_param_names:
        return ['c', 'exponent', 'amplitude']
    parvals = params.valuesdict()
    c = parvals['c']
    k = parvals['exponent']
    A = parvals['amplitude']

    y_pred = A*np.exp(k*x) + c
    if y is None:
        return y_pred
    elif weights is None:
        return (y - y_pred)
    else:
        return (y - y_pred)*weights

def poisson(params, x, y=None, weights=None, return_param_names=False):
    """"""
    if return_param_names:
        return ['lambda_param']
    
    parvals = params.valuesdict()
    lambda_param = parvals['lambda_param']


    y_pred = (np.power(lambda_param, x)*np.exp(-lambda_param)/np.array([factorial(i) for i in x])).astype(float)
    if y is None:
        return y_pred
    elif weights is None:
        return (y - y_pred)
    else:
        return (y - y_pred)*weights

def powerexp(params, x, y=None, weights=None, index=None):
    """"""
    if index is None: index = np.ones(len(x)).astype(bool)
    parvals = params.valuesdict()
    c = parvals['c']
    k = parvals['base']
    A = parvals['amplitude']

    y_pred = A*np.power(k, x) + c
    if y is None:
        return y_pred[index]
    elif weights is None:
        return (y - y_pred)[index]
    else:
        return ((y - y_pred)*weights)[index]
    
    
def binding_curve_nonlinear(params, concentrations, data=None, weights=None, index=None, return_param_names=False):
    """  Return fit value, residuals, or weighted residuals of a binding curve with nonlinear, nonspecific term.
    
    Hill coefficient 1. """
    if return_param_names:
        return ['fmax', 'dG', 'fmin', 'dGns']
    if index is None:
        index = np.ones(len(concentrations)).astype(bool)
        
    parameters = fittingParameters()
    
    parvals = params.valuesdict()
    fmax = parvals['fmax']
    dG   = parvals['dG']
    fmin = parvals['fmin']
    dG_ns = parvals['dGns']

    kd = np.exp(dG/parameters.RT)/parameters.concentration_units
    kd_ns = np.exp(dG_ns/parameters.RT)/parameters.concentration_units
    fracbound = fmin + fmax*(concentrations/(kd + concentrations))*(1 + concentrations/(kd_ns + concentrations))

    # return fit value of data is not given
    if data is None:
        return fracbound[index]
    
    # return residuals if data is given
    elif weights is None:
        return (fracbound - data)[index]
    
    # return weighted residuals if data is given
    else:
        return ((fracbound - data)*weights)[index]

   
def processFuncInputs(func_name, x, params_to_change=None, params_init=None, params_lb=None, params_ub=None, params_vary=None):
    """Return FitParameters structure given user input."""
    
    # find what the param names should be given the func_name. If func_name not recognized, assume you've specified all param_names.
    if func_name == 'binding_curve':
        param_names = ['fmin', 'dG', 'fmax']
    elif func_name == 'binding_curve_linear':
        param_names = ['fmin', 'dG', 'fmax', 'slope']
    elif func_name == 'rates_on':
        param_names = ['fmin', 'kobs', 'fmax']
    elif func_name == 'rates_off':
        param_names = ['fmin', 'koff', 'fmax']
    elif func_name == 'binding_curve_nonlinear':
        param_names = ['fmin', 'dG', 'fmax', 'dGns']
    elif func_name == 'melt_curve':
        param_names = ['fmin', 'dH', 'Tm', 'fmax']
    else:
        print "Function %s not recognized. Must have default param names for each func in objfunctions!."%func_name
        sys.exit()
    
    # change the params specified by "params_to_change", to values specified in other arguments.
    if params_to_change is None: params_to_change = []
    param_changes = {}    
    if params_init is None: params_init=[None]*len(param_names)
    if params_ub is None: params_ub=[None]*len(param_names)   
    if params_lb is None: params_lb=[None]*len(param_names)
    if params_vary is None: params_vary=[None]*len(param_names)          
    for i, param in enumerate(params_to_change):
        param_changes[param] = [params_lb[i], params_init[i], params_ub[i], params_vary[i]]

    # Go through and either use defaults defined in fitting.getFitParam, or use the ones specified above.                                                 
    fitParameters = []
    for i, param_name in enumerate(param_names):
        # add default initial values
        if param_name in param_changes.keys():
            lowerbound, init_val, upperbound , vary = param_changes[param]
        else:
            lowerbound, init_val, upperbound , vary = None, None, None, None
        # load defaults    
        fitParameters.append(fitting.getFitParam(param_name, concentrations=x, init_val=init_val, vary=vary, ub=upperbound, lb=lowerbound))
    fitParameters = pd.concat(fitParameters, axis=1)
    return fitParameters

  
