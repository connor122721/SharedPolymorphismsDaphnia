#!/usr/bin/env bash
#
#SBATCH -J snpEFF # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/snpeff.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/snpeff.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# This script will run snpEFF to annotate a VCF.

# Load Modules
module load tabix/0.2.6
module load gcc/9.2.0 htslib/1.10.2

# Working folder is core folder where this pipeline is being run.
wd="/project/berglandlab/connor/new_vcf2"

# Combined VCF name
vcf="combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.vcf.gz"

# Output annotated VCF name
out_vcf="combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf"

# java parameters
JAVAMEM=50G

# Move to working directory
cd ${wd}

# Run snpEFF on raw vcf
java -Xmx${JAVAMEM} -jar \
/home/csm6hg/SNPEFF/snpEff.jar ann \
dpgenome \
${wd}/$vcf \
-o vcf \
-t \
-s ${wd}/SNPeff.summary > \
$out_vcf

# BGzip and tab index annotated vcf
bgzip ${wd}/${out_vcf}
tabix -p vcf ${wd}/${out_vcf}.gz

# Finish
echo "Complete" $(date)
