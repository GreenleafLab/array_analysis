import numpy as np
import pandas as pd

class fittingParameters():
    """
    stores some parameters and functions
    """
    def __init__(self, concentrations=None, params=None, fitParameters=None,
                 default_errors=None):

        
        # save the units of concentration given in the binding series
        self.concentration_units = 1E-9 # i.e. nM
        self.RT = 0.582
        
        # When constraining the upper and lower bounds of dG, say you only think
        # can fit binding curves if at most it is 99% bound in the first
        # point of the binding series. This defines 'frac_bound_lowerbound'.
        # 'frac_bound_upperbound' is the minimum binding at the last point of the
        # binding series that you think you can still fit.
        self.frac_bound_upperbound = 0.01
        self.frac_bound_lowerbound = 0.99
        self.frac_bound_initial = 0.5
        
        # assume that fluorsecnce in last binding point of tightest binders
        # is on average at least 25% bound. May want to lower if doing
        # different point for binding point. 
        self.saturation_level   = 0.25
        
        # also add other things
        self.cutoff_kd = 5000
        self.cutoff_dG = self.find_dG_from_Kd(self.cutoff_kd)


    def find_dG_from_Kd(self, Kd):
        return self.RT*np.log(Kd*self.concentration_units)

    def find_Kd_from_dG(self, dG):
        return np.exp(dG/self.RT)/self.concentration_units
    
    def find_Kd_from_frac_bound_concentration(self, frac_bound, concentration):
        return concentration/float(frac_bound) - concentration

    def find_dG_from_frac_bound(self, frac_bound, concentration):
        return self.find_dG_from_Kd(self.find_Kd_from_frac_bound_concentration(frac_bound, concentration))

class fittingParametersMelt():
    """
    Stores parameters and conversion functions for melt curves
    """
    def __init__(self) -> None:
        self.energy_unit = r'kcal\mol'
        # set to positive as a convenient conversion factor
        self.absolute_zero = 273.15
        # Bolzemann constant in kcal/(mol.K)
        self.kB = 0.0019872

        # bounds for Tm
        self.Tm_lb = self.find_Kalvin_from_celsius(0)
        self.Tm_ub = self.find_Kalvin_from_celsius(150)

        # fmin fmax bound margin
        self.f_margin = 0.25

    def find_celsius_from_Kalvin(self, temperature_K):
        return temperature_K - self.absolute_zero

    def find_Kalvin_from_celsius(self, temperature_c):
        return temperature_c + self.absolute_zero

    def find_dH_from_frac_unfolded_and_Tm(self, frac_unfolded, Tm, temperature):
        dH = self.kB / (1/Tm - 1/temperature) * np.log(1/frac_unfolded - 1)
        return dH