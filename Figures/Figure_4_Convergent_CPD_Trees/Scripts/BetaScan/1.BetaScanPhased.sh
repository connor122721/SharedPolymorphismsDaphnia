#!/usr/bin/env bash
#SBATCH -J betascan
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # hours
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/betascan.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/betascan.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules to load
module load vcftools
module load tabix
module load gcc/9.2.0 htslib/1.10.2
module load anaconda/2019.10-py2.7

# Working directory
wd="/project/berglandlab/connor"

# glactools executable
glactools="/home/csm6hg/glactools/glactools"

# BetaScan excecutable
beta="/home/csm6hg/BetaScan/BetaScan.py"

# Parameter file
intervals="/project/berglandlab/connor/metadata/D84Agoodscaffstouse.12chromosomes.comma.bed"

# Parameters
chrom=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${intervals} | cut -f1 -d',')
start=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${intervals} | cut -f2 -d',')
stop=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${intervals} | cut -f3 -d',')

### EUROPE ###

# Finish chromosome
echo "Starting chromosome:" ${chrom}

# Allele counts
${glactools} \
vcfm2acf \
--onlyGT \
--fai ${wd}/totalHiCwithallbestgapclosed.fa.fai \
${wd}/for_daniel/${chrom}.euro.phased.vcf > \
${wd}/betascan/${chrom}.snps.genome.euro.acf.gz

# Fold the SFS
${glactools} \
acf2betascan \
--fold \
${wd}/betascan/${chrom}.snps.genome.euro.acf.gz | \
gzip -f > ${wd}/betascan/${chrom}.snps.genome.euro.beta.gz

# Calculate Beta
python ${beta} \
-i ${wd}/betascan/${chrom}.snps.genome.euro.beta.gz \
-fold \
-o ${wd}/betascan/${chrom}.snps.genome.euro.betafin

### NORTH AMERICAN ###

# Allele counts
${glactools} \
vcfm2acf \
--onlyGT \
--fai ${wd}/totalHiCwithallbestgapclosed.fa.fai \
${wd}/for_daniel/${chrom}.nam.phased.vcf > \
${wd}/betascan/${chrom}.snps.genome.nam.acf.gz

# Fold the SFS
${glactools} \
acf2betascan \
--fold \
${wd}/betascan/${chrom}.snps.genome.nam.acf.gz | \
gzip > ${wd}/betascan/${chrom}.snps.genome.nam.beta.gz

# Calculate Beta
python ${beta} \
-i ${wd}/betascan/${chrom}.snps.genome.nam.beta.gz \
-fold \
-o ${wd}/betascan/${chrom}.snps.genome.nam.betafin
