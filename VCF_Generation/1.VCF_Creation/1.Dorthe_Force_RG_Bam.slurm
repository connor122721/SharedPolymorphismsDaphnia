#!/usr/bin/env bash
#
#SBATCH -J Force_RG # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-03:00:00 # 3 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/Force_RG.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/Force_RG.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load picard

SLURM_ARRAY_TASK_ID=1

# Parameters
parameterFile="/scratch/csm6hg/daphnia_phylo/dorthe.bams"
wd="/scratch/csm6hg/daphnia_phylo"

# Extract sample name
samp=`sed -n ${SLURM_ARRAY_TASK_ID}p $parameterFile`
out=`echo $samp | sed 's/finalmap.bam//'`

echo "Adding read groups -" "Sample:" $SLURM_ARRAY_TASK_ID

# Move to directory
cd $wd/bam_dorthe_fin

# Force add read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I ${samp} \
-O ${out}finalmap_RG.bam \
-LB "library" \
-PL "ILLumina" \
-PU "platunit" \
-SM ${samp} 

# Index Bam files
java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
-I ${out}finalmap_RG.bam
