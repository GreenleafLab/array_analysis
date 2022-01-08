"""
Filter functions for deciding which (single clusters or variants) are good fits.
"""

import numpy as np
import pandas as pd

def default_filter(table):
    """Filter the fit parameters. """
    index = (table.rsq > 0.5)&(table.dG_stde.astype(float) < 1)&(table.fmax_stde.astype(float)<table.fmax.astype(float))
    return table.loc[index]

def filterSingleClusterFitsNNN(table):
    """
    Filter single cluster fits.
    Defines the criterion for a success cluster
    """
    query = "rsqr > 0.5 & fmin > -1 & fmin < 2 & fmin_stderr < fmin + 1 & fmax > 0 & fmax < 3 & fmax_stderr < fmax & Tm_stderr < 10 & dH_stderr < 100"
    return table.query(query)
