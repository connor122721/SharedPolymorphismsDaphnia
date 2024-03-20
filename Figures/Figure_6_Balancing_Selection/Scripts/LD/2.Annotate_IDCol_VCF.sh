#!/usr/bin/env bash
#SBATCH -J ld_cand_snps
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/ld.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/ld.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

#load modules
module load tabix
module load bcftools
module load vcftools

# Working directory
wd="/project/berglandlab/connor"

# VCF
vcf="/project/berglandlab/connor/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz"

# add annotation
#bcftools query \
#-f '%CHROM\t%POS\t%POS\t%CHROM\_%POS\_SNP\n' \
#$vcf > ${wd}/linkage/snps.meta.txt

#Index the annotation
#bgzip ${wd}/linkage/snps.meta.txt
#tabix -s1 -b2 -e2 ${wd}/linkage/snps.meta.txt.gz

# Add annotation
bcftools annotate \
-a ${wd}/linkage/snps.meta.txt.gz \
-c CHROM,FROM,TO,ID \
${vcf} > ${wd}/linkage/pulex_filt_wSNPids.vcf

#tabix
bgzip ${wd}/linkage/pulex_filt_wSNPids.vcf
tabix -p vcf ${wd}/linkage/pulex_filt_wSNPids.vcf.gz

#sanity check
bcftools query \
-f '%CHROM\t%POS\t%ID\n' \
${wd}/linkage/pulex_filt_wSNPids.vcf.gz | head

# Master Bed file
bed="/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.bed"

# Hard filter sites from vcf
vcftools \
--gzvcf ${wd}/linkage/pulex_filt_wSNPids.vcf.gz \
--exclude-bed ${bed} \
--recode \
--recode-INFO-all \
--stdout | \
bgzip > \
${wd}/linkage/pulex_filt_hard_wSNPids.vcf.gz

# Index filtered vcf
tabix \
-p vcf \
${wd}/linkage/pulex_filt_hard_wSNPids.vcf.gz

# Finisg
echo "Finish"
