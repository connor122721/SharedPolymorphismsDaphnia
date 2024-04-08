#!/usr/bin/env bash
#
#SBATCH -J RepeatMasker
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 40G
#SBATCH -o /scratch/csm6hg/err/RepeatMasker.out # Standard output
#SBATCH -e /scratch/csm6hg/err/RepeatMasker.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules to load
module load singularity

# Genome fasta
gene="/scratch/csm6hg/genomes/pulex_nam/GCF_021134715.1_ASM2113471v1_genomic.fna"
cd /scratch/csm6hg/repeats

# Make repeatmodeler DB
singularity run ~/dfam-tetools-latest.sif \
BuildDatabase \
${gene} \
-name nampulex 

# Run repeatmasker
singularity run ~/dfam-tetools-latest.sif \
RepeatModeler \
-database nampulex \
-threads 20 \
-LTRStruct
