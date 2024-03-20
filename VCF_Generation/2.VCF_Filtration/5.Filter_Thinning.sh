#!/usr/bin/env bash
#
#SBATCH -J VCF_thining
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # hours
#SBATCH --mem 80G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/ThinVCFs.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/ThinVCFs.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# This script will thin SNPS every 500 base pairs using vcftools

# Modules to load
module load vcftools
module load tabix
module load gcc/9.2.0 htslib/1.10.2

# Working directory
wd="/project/berglandlab/connor/new_vcf"

# Output directory
out="/scratch/csm6hg/daphnia_phylo/admixture"

# Input file
IN_GZVCF="/project/berglandlab/connor/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz"

# Master Bed file
bed="/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.bed"

# Move to working directory
cd ${wd}

# Thin VCF every 500 SNPs
vcftools \
--gzvcf ${IN_GZVCF} \
--exclude-bed ${bed} \
--thin 500 \
--recode \
--recode-INFO-all \
--stdout | \
bgzip > \
${out}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.thin500.vcf.gz

# Index filtered vcf
tabix \
-p vcf \
${out}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.thin500.vcf.gz

# Finisg
echo "Finish"
