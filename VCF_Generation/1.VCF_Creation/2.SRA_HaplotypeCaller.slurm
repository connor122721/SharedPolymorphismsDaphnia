#!/usr/bin/env bash
#
#SBATCH -J gatk_chrom # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00:00 # 3 days
#SBATCH --mem 25G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gatk/4.1.6.0
module load tabix/0.2.6

# Parameters
parameterFile="/scratch/csm6hg/daphnia_phylo/SRA_paramList_1"
wd="/scratch/csm6hg/daphnia_phylo"

# Extract sample name
id=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
samp=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
chrom=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

echo "Haplotype calling -" "Sample:" $SLURM_ARRAY_TASK_ID
echo "Chromosome:" ${chrom}

# Create folder for chromosome
if [[ -d "${wd}/gvcf/${chrom}" ]]
then
	echo "Working chromosome folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir ${wd}/gvcf/${chrom}
	date
fi

# Haplotype Calling
gatk HaplotypeCaller \
-R ${wd}/totalHiCwithallbestgapclosed.fa \
-I ${wd}/final_bam/${samp}_finalmap_RG.bam \
-O ${wd}/gvcf/${chrom}/${samp}.${chrom}.${id}.g.vcf \
-L ${wd}/${chrom}.bed \
-ERC GVCF

# Bgzip
module load gcc/9.2.0 htslib/1.10.2

# Compress and index with Tabix
bgzip ${wd}/gvcf/${chrom}/${samp}.${chrom}.${id}.g.vcf
tabix -p vcf ${wd}/gvcf/${chrom}/${samp}.${chrom}.${id}.g.vcf.gz

# Finish
echo "Complete -" "Sample:" $SLURM_ARRAY_TASK_ID
