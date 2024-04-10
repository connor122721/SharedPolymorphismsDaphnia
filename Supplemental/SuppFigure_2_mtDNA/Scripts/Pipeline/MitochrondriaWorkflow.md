# Mitochondrial Tree Creation for the *Daphnia pulex* Species Complex

**Author:** Connor Murray   
**Date:** 4/10/2024

## Description:

This repository contains scripts and workflows to create a mitochondrial tree for the *Daphnia pulex* species complex.

## Workflow Overview:

1. **parameterListGenerator.R**  
   Creates a sample list with locations of files etc.

2. **DownloadMapMitochondria.sh**  
   Downloads and maps samples to their respective reference genome. Outputs consensus mtDNA fasta for each individual.

3. **ProteinExonerateGenomes.md**  
   Aligns protein-coding regions across genomes. Makes GFF files for future use.

4. **Concatenate_Genomes_align.sh**  
   Makes gene-level fastas across all samples. Uses IQTREE to make both individual gene trees and mitochondrial genome tree. Outgroup is European *Daphnia magna*.

5. **AnalysisTrees.R**  
   Plots mtDNA tree.

## Finalized Tree for Protein-Coding mtDNA:

![Finalized mtDNA Tree](https://user-images.githubusercontent.com/55203772/220213780-4afc4091-430c-4515-ba87-d880e923db1a.png)

## Additional Scripts:

- **convert_exonerate_gff_to_gff3.py**  
  Conversion between exonerate output and GFF.

- **RenameIndexFastaHeaders.sh**  
  Renames fasta headers.
