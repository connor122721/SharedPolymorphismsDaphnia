#!/usr/bin/env bash

# Simplified pipeline to: 
# 1) concatenate haplotype sequences.
# 2) align sequences.
# 3) trim redundant sequences.
# 4) build a gene tree.
# 10/5/2023
# Connor Murray

# Load mafft
SLURM_ARRAY_TASK_ID=1
module load gcc/9.2.0 openmpi/3.1.6 mafft/7.475

# Working directory
wd="/scratch/csm6hg/blop"
cd $wd

# Bed file
chromFile="/scratch/csm6hg/data/blop.bed"

# Go through each gene through each window
i=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${chromFile} )

# Parameters
chr=$( echo ${i} | tr ":" '\t' | cut -f1 )
start=$( echo ${i} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f1)
stop=$( echo ${i} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f2 )
echo $chr $start $stop

# Add sample name and haplotype to header
f="${chr}.${start}.${stop}"

# Concat every fasta
cat ${wd}/fasta/*${f}.*fa > \
${wd}/${f}.fa

# Align fasta
mafft --auto \
${wd}/${f}.fa > \
${wd}/${f}.aln.fa

# Trim alignment
module load singularity
singularity run ~/sifs/clipkit_latest.sif \
clipkit -m kpi-smart-gap \
${wd}/${f}.aln.fa

# Realign
mafft --auto \
--thread 4 \
${wd}/${f}.aln.fa.clipkit > \
${wd}/${f}.aln.clip.fa 

# Run iqtree2
~/iqtree2 \
-redo \
-B 1000 \
-s ${wd}/${f}.aln.clip.fa \
-T 1 \
--prefix ${wd}/${f}

# Finish tree
echo "Finished tree"
