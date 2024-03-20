#!/usr/bin/env bash
#SBATCH -J run_poppr # A single job name for the array
#SBATCH --ntasks-per-node=5 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # 3 days
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/poppr/err/poppr2.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/poppr/err/poppr2.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Start message
echo "Start"; date

# Load modules
module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj

# Working directory
wd="/project/berglandlab/connor"

# Parameter file
paramFile=${wd}/poppr/model_paramList2

# Cores
cores=5

# Population size array
popList=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# RScript for classifying MLGs
Rscript ${wd}/scripts/3.samples_poppr_group_country.R \
${popList} \
${cores}

# End message
echo "Finish"; date
