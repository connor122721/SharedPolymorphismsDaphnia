#!/usr/bin/env bash
#
#SBATCH -J Filt_Miss
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # hours
#SBATCH --mem 60G
#SBATCH -o /project/berglandlab/connor/err/FiltMissVCF.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/FiltMissVCF.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab_standard

# This script will filter SNPs based on both:
# 1) a missingness bed file
# 2) a repetetive elements bed file
# 3) a depth, Ns, and chromosomal endpoint bed file

# Modules to load
module load vcftools
module load tabix
module load gcc/9.2.0 htslib/1.10.2

# Working directory
wd="/project/berglandlab/connor"

# Missingness bed file
missbed="/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.bed"

# Out vcf
Out_vcf="${wd}/new_vcf2/daphnia.filt.qual.miss.rep.dep.chr"

# Move to working directory
cd ${wd}

# Message
echo "Missingness filter step"

# Minimum depth
MIN_DEPTH=8
MAX_DEPTH=35

# Run bcftools
vcftools \
--gzvcf ${wd}/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.vcf.gz \
--min-meanDP $MIN_DEPTH \
--max-meanDP $MAX_DEPTH \
--exclude-bed ${missbed} \
--recode \
--recode-INFO-all \
--out ${Out_vcf}

bgzip ${Out_vcf}.recode.vcf

# Index filtered vcf
tabix -p vcf ${Out_vcf}.recode.vcf.gz

echo "Finish"
