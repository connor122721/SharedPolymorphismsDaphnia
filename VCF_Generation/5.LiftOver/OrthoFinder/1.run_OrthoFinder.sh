#!/usr/bin/env bash
#
#SBATCH -J orthofinder # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-03:00 # 3 hours
#SBATCH --mem 30G
#SBATCH -o /scratch/csm6hg/err/ortho.new.out # Standard output
#SBATCH -e /scratch/csm6hg/err/ortho.new.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab_standard

# Load conda env 
module load anaconda/2020.11-py3.8
source activate msprime_env

# Run OrthoFinder on primary protein transcripts
cd /scratch/csm6hg/genomes/proteins_species/primary_transcripts/
~/OrthoFinder_source/orthofinder.py -f . -t 15