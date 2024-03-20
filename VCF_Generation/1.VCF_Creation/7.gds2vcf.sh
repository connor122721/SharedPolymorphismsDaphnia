#!/usr/bin/env bash
#SBATCH -J gds2vcf
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # hours
#SBATCH --mem 20G
#SBATCH -o /project/berglandlab/connor/err/vcf2gds.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/vcf2gds.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules to load
module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj

# Start
echo date

# working directory
wd="/project/berglandlab/connor"

# Run script
Rscript ${wd}/scripts/7.gds2vcf.R

# bgzip and tabix
module load tabix
bgzip ${wd}/new_vcf2/daphnia.filt.mlg.genome.11.18.22.vcf
tabix -p vcf ${wd}/new_vcf2/daphnia.filt.mlg.genome.11.18.22.vcf.gz

# Finish
echo "Finish"
echo date
