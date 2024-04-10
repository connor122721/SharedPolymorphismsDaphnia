#!/usr/bin/env bash
#
#SBATCH -J CombineVCFs # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5-00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/CombineVCFs.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/CombineVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will merge all final VCFs together.

# Load modules
module load gatk/4.1.6.0

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER="/scratch/csm6hg/daphnia_phylo/vcf"

# Create VCF combine list
ls $WORKING_FOLDER/*vcf.gz > $WORKING_FOLDER/interval_paramList.list

# VCF list location
intervals=/scratch/csm6hg/daphnia_phylo/vcf/interval_paramList.list

# Combined VCF name
output="/scratch/csm6hg/daphnia_phylo/vcf/raw_vcf/raw.vcf.gz"

# Gatk parameters
JAVAMEM=100G
CPU=1

# Move to working directory
cd $WORKING_FOLDER

# Combine VCFs using MergeVcfs
gatk --java-options "-Xmx${JAVAMEM}" MergeVcfs \
-I $intervals \
-O $output

# Finish
echo "Complete" $(date)
