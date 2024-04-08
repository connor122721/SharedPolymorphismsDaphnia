#!/usr/bin/env bash
#
#SBATCH -J bedtools # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hrs
#SBATCH --mem 100G
#SBATCH -o /project/berglandlab/connor/err/bed.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/bed.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Load modules
module load gcc/9.2.0 bedtools/2.29.2

# Working & temp directory
wd="/project/berglandlab/connor/BACKUP_scratch/all_bam"

# Output bam coverage
out="/project/berglandlab/connor/mapping_stats/coverage_bams"

# Bed file with header: chrom start stop
file="/project/berglandlab/connor/metadata/interval_paramList_windows.bed"

# List of individuals
list="/project/berglandlab/connor/metadata/coverage.samples.redo"

# Sample for array jobs
f=`sed -n ${SLURM_ARRAY_TASK_ID}p ${list}`

# Memory for bedtools
mem="100G"

# Create temporary folders
if [[ -d ${out} ]]
then
	echo "Working tmp folder exist"
	echo "lets move on"
	date
else
	echo "Folder doesnt exist. Let us fix that."
	mkdir ${out}
	date
fi

# For loop for every bam - run bedtools coverage
echo "Processing:" ${f}

# European dataset
if [ -f ${wd}/Euro_bams/${f}_finalmap_mdup.bam ]; then
	bedtools coverage -a ${file} -b ${wd}/Euro_bams/${f}_finalmap_mdup.bam -iobuf ${mem} > \
  ${out}/${f}.coverage.10kbp.window
fi

# Dorthe dataset
if [ -f ${wd}/bam_dorthe_fin/${f}*RG.bam ]; then
  bedtools coverage -a ${file} -b ${wd}/bam_dorthe_fin/${f}*RG.bam -iobuf ${mem} > \
   ${out}/${f}.coverage.10kbp.window
fi

# SRA dataset
if [ -f ${wd}/final_bam/${f}*_RG.bam ]; then
	bedtools coverage -a ${file} -b ${wd}/final_bam/${f}*_RG.bam -iobuf ${mem} > \
  ${out}/${f}.coverage.10kbp.window
fi

# Finish
echo "Finish:" ${f}

# Failed runs
# cd $out
# find * -type f -size -1M | sed -r 's/.{22}$//' > ../../metadata/failed.samples.coverage

# Passed runs
# find * -type f -size +1M > passed
