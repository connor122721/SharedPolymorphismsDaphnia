#!/usr/bin/env bash
#
#SBATCH -J refBias
#SBATCH --ntasks-per-node=15 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # hours
#SBATCH --mem 60G
#SBATCH -o /project/berglandlab/connor/err/refbias.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/refbias.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules to load
module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj

# Start
echo date

# Run script
Rscript /project/berglandlab/connor/scripts/4.Reference_Bias_Calc.R

# Finish
echo "Finish"
echo date
