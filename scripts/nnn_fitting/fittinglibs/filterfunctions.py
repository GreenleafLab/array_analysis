import numpy as np
import pandas as pd

def default_filter(table):
    """Filter the fit parameters. """
    index = (table.rsq > 0.5)&(table.dG_stde.astype(float) < 1)&(table.fmax_stde.astype(float)<table.fmax.astype(float))
    return table.loc[index]