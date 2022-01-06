import numpy as np
import lmfit

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