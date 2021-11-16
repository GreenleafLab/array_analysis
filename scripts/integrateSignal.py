#!/usr/bin/env python
"""
Wrapper's of Sarah Denny's codebase for snakemake workflow for the NNN project.
Yuxi Ke, Sep 2021
"""
# import pandas as pd
from .array_tools.bin_py3.fittinglibs import processing

signal = processing.getSignalFromCPFluor(snakemake.input[0])
signal.to_csv(snakemake.output[0])

# def makeCPseriesTilefile(CPfluorfile, CPseriesfile):
#     """
#     Convert CPfluor files to CPseries one by one.
#     Args:
#         CPfluorfile - str
#         CPseriesfile - str
#         (optional)CPinfo - Dict[str: str], metadata. 
#             keys="channel", "condition", "tile"
#     """
#     signal = processing.getSignalFromCPFluor(CPfluorfile)
#     signal.to_csv(CPseriesfile)


# def mergeCPSeiresTiles(CPseriesfiles, CPseriesfilename):
#     """
#     Merge all tiles into one.
#     Args:
#         CPseriesfiles - List[str] 
#     """
#     signals = []
#     for fn in CPseriesfiles:
#         signals.append(pd.Series(pd.read_csv(fn)), name=idx)

#     signal = pd.concat(signals, axis=1)
#     signal.to_csv(CPseriesfilename)

# def makeCPseriesfromCPfluors(CPfluorfiles, CPseriesfile):
#     """
#     Merge all tiles into one CPseries file
#     Args:
#         CPfluorfiles - List[str]
#         CpseriesFile - str
#     """
    
