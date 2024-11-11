#!/usr/bin/env bash
#SBATCH -J MSMC-daphnia # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 1-00:00 # Running time of 1 day
#SBATCH --mem 15G # Memory request of 100GB
#SBATCH -o /project/berglandlab/connor/msmc/err/msmc2.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/msmc/err/msmc2.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Load Modules
module load anaconda/2020.11-py3.8
source activate msprime_env

# Working directory
wd="/project/berglandlab/connor"
threads=10

# Samples
paramFile=${wd}/msmc/mlg.1.species.country.list

# Extract constants from parameter file
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p $paramFile )

# Progress
echo ${sample} "Number:" ${SLURM_ARRAY_TASK_ID}

# Run MSMC - split
${wd}/msmc2_linux64bit \
-t ${threads} \
--fixedRecombination \
-p 10*1+15*2 \
-I 0-1 \
-o ${wd}/msmc/output/${sample}.msmc \
${wd}/msmc/Scaffold*/*${sample}*.phase.filt.samp
