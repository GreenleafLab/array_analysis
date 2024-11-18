#!/bin/bash
#SBATCH --job-nam=registration
#SBATCH --output=/scratch/groups/wjg/kyx/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,BEGIN
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00


source ~/.bashrc
conda activate py36
module load matlab
export MATLABPATH=/home/groups/wjg/kyx/array_analysis/scripts/array_tools/CPscripts/:/home/groups/wjg/kyx/array_analysis/scripts/array_tools/CPlibs/
python3 $GROUP_HOME/kyx/array_analysis/scripts/array_tools/array_data_processing/getRegistrationOffsets.py -id $OAK/kyx/data/rf003/RedFidTest/ -sd $OAK/kyx/data/rf003/FidTileTest/ -gv /home/groups/wjg/kyx/array_analysis/scripts/array_tools/ -f FID -od $OAK/kyx/data/rf003/registration -op registration_offset_test_20240217
