#!/usr/bin/env bash
#
#SBATCH -J liftover_array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # hours
#SBATCH --mem 100G
#SBATCH -o /scratch/csm6hg/err/liftc.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/liftc.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# This script will liftover chromosomal VCFs to Euro pulex

# Load Modules; SLURM_ARRAY_TASK_ID=1
module load tabix
module load picard

# Liftover files
dt="/scratch/csm6hg/mapping/vcf"

# Parameters
JAVAMEM=100G
CPU=1

# Chromosomes to analyze
intervals=/scratch/csm6hg/data/pulex_nam_chroms.txt

# Create liftover folder
if [[ -d "${d}/lift" ]]
then
echo "Working folder exist"
echo "lets move on"
date
else
echo "folder doesnt exist. lets fix that"
mkdir ${d}/lift
date
fi

# Move to working directory
cd ${dt}/lift

# Chromosome
i=$( cat ${intervals} | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f1 )

# Conditional protocol
if [ -f ${i}*am2eurolift.vcf.gz.tbi ]; then 
  echo ${i} ": This file exists! Skipping protocol."
else {

echo ${i} "Not found"

# Liftover using picard
java "-Xmx${JAVAMEM}" -jar $EBROOTPICARD/picard.jar  LiftoverVcf \
I=${dt}/${i}_filtsnps10bpindels_snps_qual.pass.dep.EuroChrms.vcf \
O=${i}.filt.am2eurolift.new.vcf.gz \
CHAIN=/scratch/csm6hg/data/american_to_european_chredit.liftOver \
REJECT=${i}.new.rejected_variants.vcf \
R=/scratch/csm6hg/ref/totalHiCwithallbestgapclosed.fa \
MAX_RECORDS_IN_RAM=50000000 \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true 

tabix -p vcf ${i}.filt.am2eurolift.new.vcf.gz

}
fi