#!/usr/bin/env bash
#SBATCH -J ADMIXTURE
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/FilterVCFs.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/FilterVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will remove Individuals and make input for admixture

# Modules to load
module load bcftools
module load tabix
module load gsl

echo "Start filtering"
echo date

# Working directory
wd="/scratch/csm6hg/daphnia_phylo/admixture"
cd ${wd}

# Input file
IN_GZVCF=${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.thin500.vcf.gz

# Individuals to keep
Individuals_to_keep=${wd}/pop.clust.sub

# Output file
Out=$( echo "${wd}/daphnia.filt.new.ann.thin500.sub.vcf" )

# Extract samples
bcftools view \
-Ov \
--samples-file $Individuals_to_keep \
$IN_GZVCF \
--threads 5 \
-o $Out

# Load plink
module load plink

# ADDMIXTURE executable
admix="/home/csm6hg/admixture_linux-1.3.0/admixture"

# Generate the input file in plink format
plink \
--vcf ${Out} \
--make-bed \
--maf 0.01 \
--geno 0.999 \
--double-id \
--out daphnia.thin.sub \
--allow-extra-chr \
--set-missing-var-ids @:#[b37]

# Thin VCF
plink \
--bed daphnia.thin.sub.bed \
--bim daphnia.thin.sub.bim \
--fam daphnia.thin.sub.fam \
--make-bed \
--double-id \
--out daphnia.new.sub \
--allow-extra-chr \
--set-missing-var-ids @:#[b37]

# ADMIXTURE does not accept chromosome names that are not human chromosomes
awk '{$1="0";print $0}' daphnia.new.sub.bim > daphnia.new.sub.bim.tmp
mv daphnia.new.sub.bim.tmp daphnia.new.sub.bim

# Finish
echo "done"
echo date
