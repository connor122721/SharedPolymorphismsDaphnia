#!/usr/bin/env bash
#
#SBATCH -J split_VCF # A single job name for the array
#SBATCH --ntasks-per-node=1 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 80G
#SBATCH -o /project/berglandlab/connor/err/splitvcf.out # Standard output
#SBATCH -e project/berglandlab/connor/err/splitvcf.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load htslib
module load tabix

# Working directory
wd="/project/berglandlab/connor/new_vcf2"

# Move into VCF directory
cd ${wd}

# Common VCF prefix used
vcf="daphnia.filtered.chr.busco"

bgzip ${wd}/${vcf}.vcf
tabix -p vcf ${wd}/${vcf}.vcf.gz

# Parameter file
intervals="/project/berglandlab/connor/metadata/goodChrom.txt"

# While loop to extract all chromosomes into seperate VCFs
while IFS= read -r line; do

  # Extract chromosome region
  tabix ${wd}/${vcf}.vcf.gz $line -h > $line.${vcf}.vcf;

  # bgzip and tabix individual chromosome
  bgzip ${wd}/${line}.${vcf}.vcf
  tabix -p vcf ${wd}/${line}.${vcf}.vcf.gz

  # Finish chromosome
  echo "Finish chromosome:" ${line}

done < ${intervals}

# End task
echo "Finished job"
