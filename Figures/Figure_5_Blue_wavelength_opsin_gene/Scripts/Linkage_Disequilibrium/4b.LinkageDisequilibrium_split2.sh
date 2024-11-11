#!/usr/bin/env bash
#SBATCH -J ld_cand_snps
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 20G
#SBATCH -o /project/berglandlab/connor/err/ld.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/ld.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will calcualte linkage disequlibrium from focal SNPs

# Modules to load
module load bcftools
module load tabix
module load plink

# Working directory
wd="/project/berglandlab/connor"
cd ${wd}/linkage

# Input file
euro_vcf=${wd}/linkage/euro_pulex.hard.vcf.gz
nam_vcf=${wd}/linkage/nam_pulex.hard.vcf.gz

# SNPs to caluclate R^2
snp=($(cat ${wd}/linkage/all_genes_tail_half))

# Go through each focal SNP
snpi=${snp[${SLURM_ARRAY_TASK_ID}]}
echo $snpi
date

# Rename
snp_name=$( echo $snpi | cut -f7 -d"/" | cut -f1 -d"." )
echo ${snp_name}
echo "Start Europe"

# Redo
#rm ${wd}/linkage/out/*${snp_name}.ld

# Generate the input file in plink format
plink \
--memory 15000 \
--vcf ${euro_vcf} \
--r2 \
--ld-snp-list ${snpi} \
--ld-window-r2 0.0 \
--geno 0.999 \
--maf 0.01 \
--double-id \
--out ${wd}/linkage/out/pulex_euro_${snp_name} \
--allow-extra-chr \
--set-missing-var-ids @:#[b37]

echo "Start North America"

# Generate the input file in plink format
plink \
--memory 15000 \
--vcf ${nam_vcf} \
--r2 \
--ld-snp-list ${snpi} \
--ld-window-r2 0.0 \
--geno 0.999 \
--maf 0.01 \
--double-id \
--out ${wd}/linkage/out/pulex_nam_${snp_name} \
--allow-extra-chr \
--set-missing-var-ids @:#[b37]

# Remove intermediate files
rm ${wd}/linkage/out/*${snp_name}*.log
rm ${wd}/linkage/out/*${snp_name}*.nosex

echo "Finish"
date
