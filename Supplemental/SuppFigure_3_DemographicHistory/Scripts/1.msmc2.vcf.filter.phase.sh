#!/usr/bin/env bash
#SBATCH -J MSMC-daphnia # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1-00:00 # Running time of 1 day
#SBATCH --mem 50G # Memory request of 100GB
#SBATCH -o /project/berglandlab/connor/msmc/err/hetsep.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/msmc/err/hetsep.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Load Modules
module load anaconda/2020.11-py3.8
source activate msprime_env
module load bcftools
module load samtools

# Working directory
wd="/project/berglandlab/connor"

# Samples
paramFile=${wd}/msmc/mlg.1.species.country.list

# Extract constants from parameter file
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p $paramFile )

# Chromosome file
intervals=${wd}/metadata/goodChrom.txt

# Reference genome
ref=${wd}/totalHiCwithallbestgapclosed.fa

# Cores
threads=1

# Split samples by chromosome
while read ichrom; do
# ichrom="Scaffold_9201_HRSCAF_10758"

# Progress message
echo "Chromosome:" ${ichrom}

# Create output folders
if [[ -d "${wd}/msmc/${ichrom}" ]]
then
echo "Working tmp folder exist"
echo "lets move on"
date
else
echo "Folder doesnt exist. Let us fix that."
mkdir ${wd}/msmc/${ichrom}
date
fi

# Bam location
bam=${wd}/BACKUP_scratch/all_bam/*/*${sample}_finalmap_*.bam

# Install whatshap (run once)
# pip3 install --user whatshap
export PATH=$HOME/.local/bin:$PATH

tabix -p vcf ${wd}/msmc/${ichrom}/${sample}.${ichrom}.split.vcf.gz

# Read backed phasing whatshap
whatshap \
phase \
-r ${ref} \
-o ${wd}/msmc/${ichrom}/${sample}.${ichrom}.phase.vcf \
--chromosome ${ichrom} \
--sample ${sample} \
${wd}/msmc/${ichrom}/${sample}.${ichrom}.split.vcf.gz \
${bam}

# Bgzip and index
bgzip ${wd}/msmc/${ichrom}/${sample}.${ichrom}.phase.vcf
tabix -p vcf ${wd}/msmc/${ichrom}/${sample}.${ichrom}.phase.vcf.gz

# Run python input conversion script
python ${wd}/msmc-tools/generate_multihetsep.py \
--mask ${wd}/msmc/${ichrom}/${sample}.${ichrom}.split.mask.bed.gz \
--negative_mask ${wd}/data/miss10.daphnia.pulex.merged.RMoutHiCGM.final.bed \
${wd}/msmc/${ichrom}/${sample}.${ichrom}.phase.vcf.gz > \
${wd}/msmc/${ichrom}/${sample}.${ichrom}.phase.filt.samp

done < ${intervals}
