#!/usr/bin/env bash
#
#SBATCH -J GenotypeGVCFs
#SBATCH --ntasks-per-node=10
#SBATCH -N 1
#SBATCH -t 7-00:00
#SBATCH --mem 80G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/vcf/err/Genotype.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/vcf/err/Genotype.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will conduct genotype calling on the GenomeDBI object

# Load Modules
module load gatk/4.1.6.0

# Parameters
JAVAMEM=80G
CPU=10

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/csm6hg/daphnia_phylo/vcf

# Reference genome
REFERENCE=/scratch/csm6hg/daphnia_phylo/totalHiCwithallbestgapclosed.fa

# Intervals to analyze
intervals=/scratch/csm6hg/daphnia_phylo/interval_DBI_paramList

# This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

# Chromosome
i=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Start
start=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Stop
stop=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Create temp folder
if [[ -d "TEMP_Daphnia_Genotype_${i}_${start}_${stop}" ]]
then
echo "Working TEMP_Daphnia_Genotype folder exist"
echo "Lets move on"
date
else
echo "Folder doesnt exist. Lets fix that"
mkdir $WORKING_FOLDER/TEMP_Daphnia_Genotype_${i}_${start}_${stop}
date
fi

echo ${i}_${start}_${stop} "is being processed" $(date)

# Identify the Genome database to genotyoe
GenomeDB_path=`echo /scratch/csm6hg/daphnia_phylo/DBI/Daphnia_DBI_${i}_${start}_${stop}`

# Genotype call the samples in the DBI merged set
gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
-R $REFERENCE \
-V gendb://$GenomeDB_path \
--tmp-dir=$WORKING_FOLDER/TEMP_Daphnia_Genotype_${i}_${start}_${stop} \
-O $WORKING_FOLDER/${i}.${start}.${stop}.vcf.gz \
--genomicsdb-use-vcf-codec \
-L ${i}:${start}-${stop}

# Remove temp folder
rm -rf $WORKING_FOLDER/vcf/TEMP_Daphnia_Genotype_${i}_${start}_${stop}

echo ${i} "done" $(date)
