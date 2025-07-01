#!/bin/bash
#SBATCH --job-nam=boostrapSingleClusterFit
#SBATCH --output=/scratch/groups/wjg/kyx/NNNlib2b_Oct6/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/NNNlib2b_Oct6/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,START
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=12
#SBATCH --time=02:00:00

source ~/.bashrc
conda activate fitting
python scripts/nnn_fitting/bootStrapFitFile.py -cf /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_single_cluster/NNNlib2b_DNA_20220314.CPfitted.gz -a /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/aligned/ConsensusReads_20220316_exact.CPannot -g /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_single_cluster/NNNlib2b_DNA_20220314_good_cluster_ind.txt                -p dH Tm -vc SEQID --query "rsqr>0.5&fmin>-1&fmin<2&fmin_stderr<fmin+1&fmax>0&fmax<3&fmax_stderr<fmax&Tm_stderr<10&dH_stderr<100" --n_samples 100
