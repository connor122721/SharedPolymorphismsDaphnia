#!/usr/bin/env bash
#
#SBATCH -J Filt_SNPs
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # hours
#SBATCH --mem 80G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/FiltSNPs2.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/FiltSNPs2.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab_standard

# This script will filter SNPs based on Gatk protocols without a reference SNP panel

# Modules to load
module load gatk/4.1.6.0

# Working directory
wd="/project/berglandlab/connor/new_vcf2"

# Parameters
JAVAMEM=50G
CPU=1

# Move to working directory
cd $wd

echo "Quality filter step"

# Run gatk program
gatk --java-options "-Xmx${JAVAMEM}" VariantFiltration \
-R /scratch/csm6hg/daphnia_phylo/totalHiCwithallbestgapclosed.fa \
-V ${wd}/combined.filtsnps10bpindels_snps_filthighmiss.recode.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O ${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.vcf.gz

echo "Quality filter step"
module load tabix
tabix -p vcf ${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.vcf.gz

# Run gatk program
gatk --java-options "-Xmx${JAVAMEM}" SelectVariants \
-R /scratch/csm6hg/daphnia_phylo/totalHiCwithallbestgapclosed.fa \
-V ${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.vcf.gz \
--exclude-filtered \
-O ${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.vcf.gz

tabix -p vcf ${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.vcf.gz
