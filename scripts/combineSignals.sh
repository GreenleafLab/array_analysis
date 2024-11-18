#!/bin/bash
#SBATCH --job-nam=combineSignals
#SBATCH --output=/scratch/groups/wjg/kyx/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,BEGIN
#SBATCH --mem-per-cpu=36G
#SBATCH --cpus-per-task=3
#SBATCH --time=0:30:00


source ~/.bashrc
conda activate py36
python3 array_tools/bin_py3/processData.py -mf /oak/stanford/groups/wjg/kyx/data/rf009/tmp/NNNlib2b_RNA_20240224.map -od /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/ --appendLibData /oak/stanford/groups/wjg/kyx/data/rf009/aligned/ConsensusReads_20240313_exact.CPseq --num_cores 3