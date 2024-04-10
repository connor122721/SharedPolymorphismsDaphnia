#!/usr/bin/env bash
#
#SBATCH -J gatk_chrom # A single job name for the array
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00:00 # 3 days
#SBATCH --mem 50G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gatk/4.1.6.0
module load tabix/0.2.6
module load picard

# Parameters
parameterFile="/scratch/csm6hg/daphnia_phylo/dorthe_paramList"
wd="/scratch/csm6hg/daphnia_phylo"

# Parameters
JAVAMEM=5G
CPU=4

# Extract sample name
samp=$( cat ${parameterFile} | grep ",$SLURM_ARRAY_TASK_ID$" | cut -d',' -f1 )
chrom=$( cat ${parameterFile} | grep ",$SLURM_ARRAY_TASK_ID$" | cut -d',' -f2 )

echo "Haplotype calling -" "Sample:" $SLURM_ARRAY_TASK_ID
echo "Chromosome:" ${chrom}

# Create folder for chromosome
if [[ -d "${wd}/Dorthe_gvcf/${chrom}" ]]
then
	echo "Working chromosome folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir ${wd}/Dorthe_gvcf/${chrom}
	date
fi

echo "Adding read groups -" "Sample:" $SLURM_ARRAY_TASK_ID

# Create forced read group bam
if [[ -f "${wd}/bam_dorthe_fin/${samp}_finalmap_RG.bam" ]]
then
	echo "Bam file does exist"
	date
else
	echo "Bam file doesnt exist. Lets fix that"
	java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
	-I ${wd}/bam_dorthe_fin/${samp}_finalmap.bam \
	-O ${wd}/bam_dorthe_fin/${samp}_finalmap_RG.bam \
	-LB "library" \
	-PL "ILLumina" \
	-PU "platunit" \
	-SM ${samp}
	date
fi

# Index bam
if [[ -f "${wd}/bam_dorthe_fin/${samp}_finalmap_RG.bai" ]]
then
	echo "Index for RG does exist"
	date
else
	echo "Index doesnt exist. lets fix that"
	java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
	-I ${wd}/bam_dorthe_fin/${samp}_finalmap_RG.bam
	date
fi


echo "Haplotype calling -" "Sample:" $SLURM_ARRAY_TASK_ID

# Haplotype Calling
gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" HaplotypeCaller \
-R ${wd}/totalHiCwithallbestgapclosed.fa \
-I ${wd}/bam_dorthe_fin/${samp}_finalmap_RG.bam \
-O ${wd}/Dorthe_gvcf/${chrom}/${samp}.${chrom}.${SLURM_ARRAY_TASK_ID}.g.vcf \
-L ${wd}/data/${chrom}.bed \
--native-pair-hmm-threads $CPU \
-ERC GVCF

# Bgzip
module load gcc/9.2.0 htslib/1.10.2

# Compress and index with Tabix
bgzip ${wd}/Dorthe_gvcf/${chrom}/${samp}.${chrom}.${SLURM_ARRAY_TASK_ID}.g.vcf
tabix -p vcf ${wd}/Dorthe_gvcf/${chrom}/${samp}.${chrom}.${SLURM_ARRAY_TASK_ID}.g.vcf.gz

# Finish
echo "Complete -" "Sample:" $SLURM_ARRAY_TASK_ID
