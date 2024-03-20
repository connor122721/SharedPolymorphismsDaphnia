#!/usr/bin/env bash
#
#SBATCH -J fixed_differences # A single job name for the array
#SBATCH --ntasks-per-node=15 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # 3 days
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/fixed/fixed.out # Standard output
#SBATCH -e /project/berglandlab/connor/fixed/fixed.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Start message
echo "Start"; date

# Load modules
module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj

# RScript for fixed differences
Rscript /project/berglandlab/connor/scripts/1.fixed_differences_par.R

# End message
echo "Finish"; date
