#!/usr/bin/env bash
#SBATCH -J PLINK
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/FilterVCFs.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/FilterVCFs.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# This script will remove Individuals and make input for LD

# Modules to load
module load tabix
module load gsl
module load bcftools

# Working directory
wd="/project/berglandlab/connor"
cd ${wd}/linkage

# Input file
vcf=${wd}/linkage/pulex_filt_hard_wSNPids.vcf.gz

#  Metadata
euro_list=${wd}/linkage/euro_pulex.txt
nam_list=${wd}/linkage/nam_pulex.txt

cat $euro_list | cut -f1 > ${wd}/linkage/euro_pulex.list
cat $nam_list | cut -f1 > ${wd}/linkage/nam_pulex.list

# Extract samples
bcftools view \
-Ov \
--threads 5 \
--samples-file ${wd}/linkage/euro_pulex.list \
${vcf} \
-o ${wd}/linkage/euro_pulex.hard.vcf

bgzip ${wd}/linkage/euro_pulex.hard.vcf
tabix -p vcf ${wd}/linkage/euro_pulex.hard.vcf.gz

# Extract samples
bcftools view \
-Ov \
--threads 5 \
--samples-file ${wd}/linkage/nam_pulex.list \
${vcf} \
-o ${wd}/linkage/nam_pulex.hard.vcf

bgzip ${wd}/linkage/nam_pulex.hard.vcf
tabix -p vcf ${wd}/linkage/nam_pulex.hard.vcf.gz
