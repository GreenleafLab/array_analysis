import numpy as np
import lmfit
import scipy.stats as st

class MeltCurveModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        self.T = kwargs['T']
        self.f_margin = kwargs['f_margin']

        self.kB = 0.0019872

        self.nan_policy = 'omit'
        
        def melt_curve(T, fmin, fmax, dH, Tm):
            frac_unfolded = (fmin + (fmax - fmin)/
                 (1 + np.exp((dH/self.kB) * (1/Tm - 1/T))))
            return frac_unfolded
        
        super(MeltCurveModel, self).__init__(melt_curve, *args, **kwargs)

    def guess(self, data, **kwargs):
        def guess_Tm(data):
            idmin = np.argmin(np.abs(data-0.5))
            return self.T[idmin]
        
        params = self.make_params()
        def pset(param, value, minimum=-np.inf, maximum=np.inf):
            params["%s%s" % (self.prefix, param)].set(value=value, min=minimum, max=maximum)

        pset("fmin", 0, minimum=-self.f_margin, maximum=self.f_margin)
        pset("fmax", 1, minimum=1-self.f_margin, maximum=1+self.f_margin)
        if 'dH' in kwargs.keys():
            pset("dH", kwargs["dH"], maximum=0)
        else:
            pset("dH", -40, maximum=0)
            
        if 'Tm' in kwargs.keys():
            pset("Tm", kwargs["Tm"])
        else:
            pset("Tm", guess_Tm(data))
        
        return lmfit.models.update_param_vals(params, self.prefix, **kwargs)


class GammaModel(lmfit.Model):
    def __init__(self, var_name, *args, **kwargs):
        self.var_name = var_name
        self.nan_policy = 'omit'
        
        def gamma(x, mean, std, offset):
            k, theta = self.returnGammaParams(mean, std)
            cdf = st.gamma.cdf(x, k, scale=theta, loc=offset)
            return cdf
        
        super(GammaModel, self).__init__(gamma, *args, **kwargs)


    @staticmethod
    def returnGammaParams(mean, std):
        """Return shape and scale parameters from moments of distrubtion mean, std"""
        k = (mean/std)**2
        theta = (std)**2/mean
        return k, theta

    def getParams(self, **kwargs):
        params = self.make_params()
        def pset(param, value, minimum=-np.inf, maximum=np.inf, vary=True):
            params["%s%s" % (self.prefix, param)].set(value=value, min=minimum, max=maximum, vary=vary)

        pset("mean", 1)
        pset("std", 1, minimum=0)

        if ("set_offset" in kwargs.keys()):
            pset("offset", kwargs["set_offset"], vary=False)
        else:
            pset("offset", 0)

        return lmfit.models.update_param_vals(params, self.prefix, **kwargs)


class SigmaNModel(lmfit.Model):
    def __init__(self, *args, **kwargs):

        self.nan_policy = 'omit'
        
        def sigma_n(n, a, b):
            return a / np.sqrt(n) + b
        
        super(SigmaNModel, self).__init__(sigma_n, *args, **kwargs)


    def guess(self, **kwargs):
        params = self.make_params()
        def pset(param, value, minimum=-np.inf, maximum=np.inf, vary=True):
            params["%s%s" % (self.prefix, param)].set(value=value, min=minimum, max=maximum, vary=vary)

        pset("a", 4, minimum=0)
        pset("b", 0)

        return lmfit.models.update_param_vals(params, self.prefix, **kwargs)