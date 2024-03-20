#!/usr/bin/env bash
#
#SBATCH -J bamqc # A single job name for the array
#SBATCH --ntasks-per-node=1 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 5 hours
#SBATCH --mem 2G
#SBATCH -o /project/berglandlab/connor/BACKUP_scratch/all_bam/err/bam.out # Standard output
#SBATCH -e /project/berglandlab/connor/BACKUP_scratch/all_bam/err/bam.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Working directory
cd /project/berglandlab/connor/BACKUP_scratch/all_bam

# BAMs
files=($(ls */*.bam)) 
taski=${files[${SLURM_ARRAY_TASK_ID}]}
name=$( echo ${taski} | cut -f2 -d"/" | cut -f1 -d"." )
echo ${name}

# Run flagstat
samtools flagstat ${taski} > \
bamqcout/${name}.stats
