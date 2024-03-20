#!/usr/bin/env bash
#SBATCH -J admixture # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-00:00 # 2 hours
#SBATCH --mem 30G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/admix.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/admix.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load SMC++ & Modules
module load plink

# Working & temp directory
wd="/scratch/csm6hg/daphnia_phylo/admixture"
admix="/home/csm6hg/admixture_linux-1.3.0/admixture"

# Output
out="daphnia.new.sub"

# Move to working directory
cd ${wd}/data_sub

# Go through various K values
echo "K:" ${SLURM_ARRAY_TASK_ID}

# Run admixture
${admix} \
${wd}/${out}.bed \
${SLURM_ARRAY_TASK_ID} \
--cv \
-j5 \
-B100 > \
${wd}/data_sub/log${SLURM_ARRAY_TASK_ID}.sub.out
