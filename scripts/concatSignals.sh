#!/bin/bash
#SBATCH --job-nam=concatTilesSignal
#SBATCH --output=/scratch/groups/wjg/kyx/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf,owners
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,BEGIN
#SBATCH --mem-per-cpu=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00


source ~/.bashrc
conda activate py36
python3 scripts/concatTilesSignal.py --tiles /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile1_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile2_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile3_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile4_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile5_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile6_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile7_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile8_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile10_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile11_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile12_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile13_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile14_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile15_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile16_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile17_green_600ms.CPseries.json.csv /oak/stanford/groups/wjg/kyx/data/rf009/series_20240324/CPseries/RNAmelt_tile18_green_600ms.CPseries.json.csv -o /oak/stanford/groups/wjg/kyx/data/rf009/series_merged/NNNlib2b_RNA_20240324_v0.pkl -m config/nnnlib2b_20240324_map.csv