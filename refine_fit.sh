#!/bin/bash
#SBATCH --job-nam=refineFit
#SBATCH --output=/scratch/groups/wjg/kyx/NNNlib2b_Oct6/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/NNNlib2b_Oct6/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,START
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00

source ~/.bashrc
conda activate fitting
python scripts/nnn_fitting/refineVariantFits.py --parallel -b /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/series_normalized/NNNlib2b_DNA_20211022_normalized.pkl -vf /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_single_cluster/NNNlib2b_DNA_20211022.CPvariant -x /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/series_normalized/NNNlib2b_DNA_20211022_xdata.txt        --cluster_annot /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/aligned/ConsensusReads_20220115_exact.CPseq --fmax_fmin /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_fmax_fmin/NNNlib2b_DNA_20211022-fmax_fmin.json -o /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_variant/NNNlib2b_DNA_20211022_v4.CPvariant.gz --figdir /scratch/groups/wjg/kyx/NNNlib2b_Oct6/fig/NNNlib2b_DNA_20211022-fit_refine_variant/        --param dH Tm --variant_col SEQID --n_bootstraps 100 --good_clusters /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_single_cluster/NNNlib2b_DNA_20211022_good_cluster_ind.txt --variant_filter "rsqr>0.5&fmin_init>-1&fmin_init<2&fmax_init>0&fmax_init<2&dH_init>-400"