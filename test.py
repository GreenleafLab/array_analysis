import pandas as pd
from array_fitting_tools.bin_py3.fittinglibs import processing

signal = processing.getSignalFromCPFluor("./data/CPfluor/test.CPfluor")
signal.to_csv("./data/CPseries/test.CPseries")
