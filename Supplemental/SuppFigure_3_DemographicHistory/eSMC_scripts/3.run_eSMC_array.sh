#!/usr/bin/env bash
#
#SBATCH -J run_eSMC_chrom
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-24:00 # hours
#SBATCH --mem 15G
#SBATCH -o /project/berglandlab/connor/err/run_eSMC.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/run_eSMC.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules; SLURM_ARRAY_TASK_ID=183
module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

# Message
echo "Start" 
date

# Go to wd
wd="/project/berglandlab/connor"
cd ${wd}

# Samples
paramFile=${wd}/eSMC.1.list

# Extract constants from parameter file
file=$( sed -n ${SLURM_ARRAY_TASK_ID}p $paramFile | cut -f1 )
chrom=$( sed -n ${SLURM_ARRAY_TASK_ID}p $paramFile | cut -f3 )
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p $paramFile | cut -f4 )
echo ${sample} ${chrom}

# Run script
Rscript 2a.run_eSMC.R \
${file} \
${sample} \
${chrom}

# Finish
echo "Finish" ${SLURM_ARRAY_TASK_ID}
date
