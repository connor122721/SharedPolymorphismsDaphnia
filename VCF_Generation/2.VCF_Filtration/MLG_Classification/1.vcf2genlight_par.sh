#!/usr/bin/env bash
#
#SBATCH -J vcf2genlight
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 120G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/vcf2genlight_par.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/vcf2genlight_par.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# Convert a VCF to a Genlight object

# Modules to load
module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj

# Run script
Rscript vcf2genlight_par.R
