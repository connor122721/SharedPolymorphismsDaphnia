#!/usr/bin/env bash
#
#SBATCH -J SNP_filtering
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 80G
#SBATCH -o /project/berglandlab/connor/err/FilterVCFs.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/FilterVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# This script will remove Individuals with high missingness

# Modules to load
module load vcftools
module load tabix
module load gcc/9.2.0 htslib/1.10.2

# Working directory
wd=/project/berglandlab/connor

# Individuals to remove - high missingness
Individuals_to_remove=/project/berglandlab/connor/metadata/high_missingness_samples_raw

# Individuals to remove - pooled samples
Individuals_to_remove2=/project/berglandlab/connor/metadata/pooledsamples.list

# Individuals to remove - duplicated samples from Europe
Individuals_to_remove3=/project/berglandlab/connor/metadata/duplicated.blank.Euro.samples.txt

# Input file
IN_GZVCF=${wd}/vcf/combined.filtsnps10bpindels_snps.vcf.gz

# Output file
Out=${wd}/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss

# Move to working directory
cd ${wd}

vcftools \
--gzvcf $IN_GZVCF \
--remove $Individuals_to_remove \
--remove $Individuals_to_remove2 \
--remove $Individuals_to_remove3 \
--recode \
--recode-INFO-all \
--out $Out

bgzip $Out.recode.vcf

# Index filtered vcf
tabix -p vcf $Out.recode.vcf.gz
