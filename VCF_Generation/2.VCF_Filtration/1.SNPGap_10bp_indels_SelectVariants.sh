#!/usr/bin/env bash
#
#SBATCH -J SNP_filtering
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/err/FilterVCFs.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/err/FilterVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will remove SNPS within 10 base pairs of indels by a chromosome basis using bcftools
# Also filters out only SNPs using gatk

# Modules to load
module load bcftools/1.9
module load gatk/4.1.6.0

# Working directory
wd=/scratch/csm6hg/daphnia_phylo/vcf

# Parameters
JAVAMEM=10G
CPU=4

# Intervals to analyze
intervals=/scratch/csm6hg/daphnia_phylo/interval_DBI_paramList

# Chromosome
chrom=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Start
start=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Stop
stop=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Input file
IN_GZVCF=/scratch/csm6hg/daphnia_phylo/vcf/${chrom}.${start}.${stop}.vcf.gz

# Move to working directory
cd ${wd}

# Remove SNPs 10 bp from indels
bcftools filter \
	--SnpGap 10 \
	${IN_GZVCF} \
	--output ${wd}/raw_vcf/${chrom}.${start}.${stop}_filtsnps10bpindels.vcf \
	-O z \
	--threads ${CPU}

module load tabix
module load gcc/9.2.0 htslib/1.10.2

# Index filtered vcf
tabix -p vcf ${wd}/raw_vcf/${chrom}.${start}.${stop}_filtsnps10bpindels.vcf.gz

# Filter only SNPs
gatk --java-options "-Xmx${JAVAMEM}" SelectVariants \
	-V ${wd}/raw_vcf/${chrom}.${start}.${stop}_filtsnps10bpindels.vcf.gz \
	-select-type SNP \
	-R /scratch/csm6hg/daphnia_phylo/totalHiCwithallbestgapclosed.fa \
	-O ${wd}/raw_vcf/${chrom}.${start}.${stop}_filtsnps10bpindels_snps.vcf.gz
