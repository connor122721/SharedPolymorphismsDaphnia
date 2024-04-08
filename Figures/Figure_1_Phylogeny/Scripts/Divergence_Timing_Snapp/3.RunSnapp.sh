#!/usr/bin/env bash
#
#SBATCH -J run_beast2 # A single job name for the array
#SBATCH --ntasks-per-node=20 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # 3 days
#SBATCH --mem 100G
#SBATCH -o /project/berglandlab/connor/snapp5/snappRun.out # Standard output
#SBATCH -e /project/berglandlab/connor/snapp5/snappRun.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load bcftools
module load ruby

# Working directory
wd="/project/berglandlab/connor/snapp5"
cd ${wd}

# Beast2 directory
beast2="/home/csm6hg/beast/bin/beast"

# Run Beast2
${beast2} \
-threads 20 \
${wd}/snapp.xml

# Add posterior probabilities
/home/csm6hg/beast/bin/treeannotator \
-burnin 10 \
-heights mean \
${wd}/snapp.trees \
${wd}/snapp.tre
