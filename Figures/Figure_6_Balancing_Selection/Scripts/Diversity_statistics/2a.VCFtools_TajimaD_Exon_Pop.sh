#!/usr/bin/env bash
#SBATCH -J exon_fst_pi_vcftools
#SBATCH --ntasks-per-node=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # hours
#SBATCH --mem 1G
#SBATCH -o /project/berglandlab/connor/err/exon_vcf.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/exon_vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# start
echo "Start" $(date)

# Exon Functions
exonFunc () {

# Load Modules
module load vcftools
module load tabix

# Working folder is core folder where this pipeline is being run.
wd=/project/berglandlab/connor/new_vcf2

# Input file
IN_GZVCF=${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz

# Master Bed file
bed=/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed

# Meta Exon
exon=($( ls /project/berglandlab/connor/data/exon/* ))

# Move to working directory
cd ${wd}

#i="Scaffold_7757_HRSCAF_8726       3472398 3472582 exon    Daphnia04572-RA"

# Start
echo ${i}
date

# Extract metadata
chrom=$( echo ${i} | cut -f1 -d " " )

# VCF functions
analy="--keep species.nam.TEX.pop.fst.txt \
--TajimaD 1"

# Output name
out_namey="TEX_nam_tajimaD_nomiss"

# Subset VCF and Run VCFTools
tabix -h ${IN_GZVCF} ${chrom} |
vcftools \
--vcf - \
--max-missing 1 \
--chr ${chrom} \
--exclude-bed ${bed} \
`echo ${analy}` \
--out ${wd}/vcftools/${out_namey}

echo "VCF completed" $(date)
