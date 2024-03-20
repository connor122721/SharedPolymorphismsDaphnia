#!/usr/bin/env bash
#
#SBATCH -J CombineVCFs # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/CombineVCFs.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/CombineVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will merge all final VCFs together.

# Load modules
module load gatk

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER="/scratch/csm6hg/mapping/vcf/lift"
cd $WORKING_FOLDER

# Get sample names per VCF
module load bcftools

# Create VCF combine list
ls $WORKING_FOLDER/*.am2eurolift.vcf.gz > $WORKING_FOLDER/interval_paramList.list

# VCF list location
intervals=$WORKING_FOLDER/interval_paramList.list

# Combined VCF name
output="$WORKING_FOLDER/namPulex.filt.new.am2euro.vcf.gz"

# Gatk parameters
JAVAMEM=100G
CPU=1

# Move to working directory
cd ${WORKING_FOLDER}

# Combine VCFs using MergeVcfs
gatk MergeVcfs \
-I ${intervals} \
-O ${output}

# Finish
echo "Complete" $(date)
#bcftools query -f '%CHROM %POS[\t%GT]\n' namPulex.filt.am2euro.vcf.gz > napulex.snps.genotype.lifted