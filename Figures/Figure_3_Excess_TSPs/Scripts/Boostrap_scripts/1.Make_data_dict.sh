#!/usr/bin/env bash
#
#SBATCH -J data_pickle # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1day
#SBATCH --mem 30G
#SBATCH -o /scratch/csm6hg/err/pick.out # Standard output
#SBATCH -e /scratch/csm6hg/err/pick.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Working directory
cd /scratch/csm6hg/moments/

# Modules
module load anaconda/2023.07-py3.11
source activate msprime_env

# Run script
date
python3 1.Make_data_dict.py

# Finish
date