#!/bin/bash
#SBATCH --job-nam=registration
#SBATCH --output=/scratch/groups/wjg/kyx/NNNlib2b_Nov11/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/NNNlib2b_Nov11/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,BEGIN
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=6
#SBATCH --time=1:00:00


source ~/.bashrc
conda activate py36
module load matlab
export MATLABPATH=/home/groups/wjg/kyx/array_analysis/scripts/array_tools/CPscripts/:/home/groups/wjg/kyx/array_analysis/scripts/array_tools/CPlibs/
python3 $GROUP_HOME/kyx/array_analysis/scripts/array_tools/array_data_processing/getRegistrationOffsets.py -id $OAK/kyx/data/rf003/Red09_Fid3/ -sd $OAK/kyx/data/rf003/filtered_tiles/ -gv /home/groups/wjg/kyx/array_analysis/scripts/array_tools/ -f FID -od $OAK/kyx/data/rf003/registration -op registration_offset_20240216
