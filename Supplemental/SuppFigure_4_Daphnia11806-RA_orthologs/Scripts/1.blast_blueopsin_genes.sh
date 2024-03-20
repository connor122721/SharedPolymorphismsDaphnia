#!/usr/bin/env bash

# Modules to load
module load blast/2.13.0
module load samtools

# Start
echo date

# Working directory
cd /project/berglandlab/connor

# TCL2 gene
"Scaffold_9200_HRSCAF_10757:6676234-6686097"

# Blue wavelength opsin
"Scaffold_1931_HRSCAF_2197:6350219-6351796"

# Extract sequence from fa
samtools faidx \
totalHiCwithallbestgapclosed.fa \
"Scaffold_1931_HRSCAF_2197:6350219-6351796" >
Daphnia11806-RA.fa

# Extract translated fasta from genes
~/seqtk/seqtk subseq \
"Daphnia.proteins.aed.0.6.fasta" \
"toll/tcl2" > \
"toll/tcl2.aa.fa"

# Blast sequence
blastn -query Daphnia11806-RA.fa \
-subject totalHiCwithallbestgapclosed.fa > \
blue_opsin_pulex.blast.out

# Finish
echo "Finish"
echo date
