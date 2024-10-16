#!/bin/bash
#SBATCH --job-nam=refineFit
#SBATCH --output=/scratch/groups/wjg/kyx/NNNlib2b_Nov11/out/%x_%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/NNNlib2b_Nov11/out/%x_%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,START
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00

source ~/.bashrc
conda activate fitting
# python scripts/nnn_fitting/refineVariantFits.py --parallel -b /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/series_normalized/NNNlib2b_DNA_20211022_normalized.pkl -vf /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_single_cluster/NNNlib2b_DNA_20211022.CPvariant -x /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/series_normalized/NNNlib2b_DNA_20211022_xdata.txt        --cluster_annot /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/aligned/ConsensusReads_20220115_exact.CPseq --fmax_fmin /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_fmax_fmin/NNNlib2b_DNA_20211022-fmax_fmin.json -o /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_variant/NNNlib2b_DNA_20211022_v6.1.CPvariant.gz --figdir /scratch/groups/wjg/kyx/NNNlib2b_Oct6/fig/NNNlib2b_DNA_20211022-fit_refine_variant/        --param dH Tm --variant_col SEQID --n_bootstraps 100 --good_clusters /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fitted_single_cluster/NNNlib2b_DNA_20211022_good_cluster_ind.txt --variant_filter "rsqr>0.5&fmin_init>-1&fmin_init<2&fmax_init>0&fmax_init<2&dH_init>-400" --mapfile /home/groups/wjg/kyx/array_analysis/config/nnnlib2b_map.csv
# python scripts/nnn_fitting/refineVariantFits.py --parallel -b /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/series_normalized/NNNlib2b_DNA_20211221_normalized.pkl -vf /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_single_cluster/NNNlib2b_DNA_20211221.CPvariant -x /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/series_normalized/NNNlib2b_DNA_20211221_xdata.txt        --cluster_annot /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/aligned/ConsensusReads_20211118_exact.CPseq --fmax_fmin /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_fmax_fmin/NNNlib2b_DNA_20211221-fmax_fmin.json -o /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_variant/NNNlib2b_DNA_20211221_v6.CPvariant.gz --figdir /scratch/groups/wjg/kyx/NNNlib2b_Nov11/fig/NNNlib2b_DNA_20211221-fit_refine_variant/        --param dH Tm --variant_col SEQID --n_bootstraps 100 --good_clusters /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_single_cluster/NNNlib2b_DNA_20211221_good_cluster_ind.txt --variant_filter "rsqr>0.5&fmin_init>-1&fmin_init<2&fmax_init>0&fmax_init<2&dH_init>-400" --mapfile /home/groups/wjg/kyx/array_analysis/config/nnnlib2b_20211221_map.csv
python3 scripts/nnn_fitting/refineVariantFits.py \
    -b /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/series_normalized/NNNlib2b_DNA_20220314_normalized.pkl \
    -vf /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_single_cluster/NNNlib2b_DNA_20220314.CPvariant \
    -x /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/series_normalized/NNNlib2b_DNA_20220314_xdata.txt \
    --mapfile config/nnnlib2b_20220314_map.csv \
    --cluster_annot /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/aligned/ConsensusReads_20220316_exact.CPannot \
    --fmax_fmin /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_fmax_fmin/NNNlib2b_DNA_20220314-fmax_fmin.json \
    -o /scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/fitted_variant/NNNlib2b_DNA_20220314_v2.CPvariant.gz \
    --figdir /scratch/groups/wjg/kyx/NNNlib2b_Nov11/fig/NNNlib2b_DNA_20220314-fit_refine_variant/ \
    --param dH Tm --variant_col SEQID --n_bootstraps 100 --parallel
