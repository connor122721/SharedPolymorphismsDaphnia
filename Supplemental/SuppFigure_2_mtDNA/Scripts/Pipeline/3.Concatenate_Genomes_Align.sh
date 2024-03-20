#!/usr/bin/env bash
#SBATCH -J mtdna.Tree
#SBATCH --ntasks-per-node=5 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 day
#SBATCH --mem 50G
#SBATCH -o /scratch/csm6hg/mito/err/mtdna.Tree.out # Standard output
#SBATCH -e /scratch/csm6hg/mito/err/mtdna.Tree.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gcc/9.2.0 openmpi/3.1.6 mafft/7.475
wd="/scratch/csm6hg/mito"

# Creates gene list (run once)
#grep ">" /scratch/csm6hg/Daphnia_mtdna/euro_pulex/mtdnaD8119.faa | cut -f4 -d";" | sed "s/ //g"  > \
#/scratch/csm6hg/mito/mtdna.cds.genes

# Gene list
genes="/scratch/csm6hg/mito/mtdna.cds.genes"

# Go through each gene
while IFS= read -r prefix; do
  #prefix="atp6"
  echo ${prefix}

  cat ${wd}/cdsgenes_fasta/*${prefix}.fa > \
  ${wd}/cdsgenes_fasta/${prefix}.fa

  mafft --auto \
  ${wd}/cdsgenes_fasta/${prefix}.fa > \
  ${wd}/align/${prefix}.aln.fa

  # Outgroup -> Euro D. magna
  out="mtdna_D8_119.Daphnia_magna.${prefix}"

  # Run iqtree2
  ~/iqtree2 -redo \
  -B 1000 \
  -o ${out} \
  -s ${wd}/align/${prefix}.aln.fa \
  -T 1 \
  --prefix ${wd}/tree/${prefix}.aln

done < ${genes}

# Combine all aligned fastas
while IFS= read -r prefix; do
  #prefix="cox1"
  echo ${prefix}

  # Remove gene names from fasta
  words=("$prefix")
  fmt=$(IFS='|'; printf '%q' "s/${words[*]}//g")
  sed -i "$fmt" ${outwd}/align/${prefix}.sub.aln.fa

done < ${genes}

# Combine all gene sequences
~/seqkit concat \
${outwd}/align/*.sub.aln.fa > \
${outwd}/align/mtdna.12cds.sub.fa

### mtDNA Tree ###

# Outgroup -> Euro D. magna
out="mtdna_D8_119.Daphnia_magna."

# Run iqtree2
~/iqtree2 \
-redo \
-B 1000 \
-o ${out} \
-s ${outwd}/align/mtdna.12cds.sub.fa \
-T 1 \
--prefix /scratch/csm6hg/mito/tree/mtdna.12cds.sub.aln
