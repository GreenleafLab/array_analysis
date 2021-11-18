"""
Test script for `plot_seqs.py`
"""

from plot_seqs import *
import os

CPseqs = [os.path.join('/scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/filtered_tiles', 'ALL_tile%03d_Bottom_filtered.CPseq'%i) for i in range(1,19)]
pngs = [os.path.join('/scratch/groups/wjg/kyx/NNNlib2b_Oct6/fig/fiducial', 'tile%03d_Bottom_fiducial.png'%i) for i in range(1,19)]
plot_fiducial_all_tiles(CPseqs, pngs)
