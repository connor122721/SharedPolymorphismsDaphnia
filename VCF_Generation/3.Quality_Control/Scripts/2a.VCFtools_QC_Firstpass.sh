#! /bin/bash
#
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/err/FinQCVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/err/FinQCVCF.%A_%a.err # Standard error

# This script will check various QC statistics from a raw VCF file

# Load Modules
module load vcftools

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/csm6hg/daphnia_phylo/vcf

# Output
PIPELINE=/scratch/csm6hg/daphnia_phylo/vcf/raw_vcf

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
cd $WORKING_FOLDER

# Survey quality of final VCF
analyses=("--depth" \
"--site-mean-depth" \
"--site-quality" \
"--missing-indv" \
"--missing-site"
"--het")

# For loop to run through 6 QC steps
for i in {0..5}

do

echo "Now processing" ${analyses[${i}]}
echo "VCF:" ${chrom}.${start}.${stop}

# Run VCFTools
vcftools \
--gzvcf $IN_GZVCF \
`echo ${analyses[${i}]}` \
--out $PIPELINE/${chrom}.${start}.${stop}

# Finish i
done

echo "VCF QC completed" $(date)
