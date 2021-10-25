#!/bin/bash
#SBATCH --job-name=quantify_images
#################  
#a file for job output, you can check job progress
#SBATCH --output=quantify_tiles.out
#################
# a file for errors from the job
#SBATCH --error=quantify_tiles.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=1:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=wjg,biochem,owners
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-task=16" instead
# Sherlock nodes have 16 cpus. For some reason, you can request more than that on 'owners' partition, but not on others. 
# It can't give you more obviously, but it's interesting that it lets you ask on one partion but not another.
# Note: On Sherlock2, most of the nodes we have access to have 20 cores. 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=24G
#################

#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END

source ~/.bashrc
conda activate ame
module load matlab
export MATLABPATH=scripts/array_tools/CPscripts:scripts/array_tools/CPlibs
python scripts/array_tools/CPscripts/quantifyTilesDownstream.py -id /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/images/Green15_25/ -ftd /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/filtered_tiles_libregion/ -fd /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/fluor/ -rod /scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/roff/ -n 18 -rs LibRegion -sf MiSeq_to_TIRFStation1 -gv scripts/array_tools/ 
