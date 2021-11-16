#!/bin/bash
#SBATCH --job-nam=makeFilteredTiles
#SBATCH --mail-user=kyx@stanford.edu
source ~/.bashrc
conda activate snakemake
snakemake --profile slurm --use-conda
