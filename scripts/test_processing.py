import pandas as pd
from array_tools.bin_py3.fittinglibs import processing

processing.makeCPseriesFile('./cpseriestest.csv', pd.Series(['/scratch/users/kyx/fluor/07_Green_Fiducial/Sequence9_tile1_green_2_600ms_2011.07.26-20.00.22.531.CPfluor']), '/scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/aligned/ConsensusReads_20211012_exact.CPseq')
