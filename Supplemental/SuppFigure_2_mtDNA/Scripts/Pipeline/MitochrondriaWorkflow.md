# This creates the mitochondrial tree for the Daphnia pulex species complex
# Connor Murray (csm6hg@virginia.edu)
# 2/20/2023

## Creates sample list with locations of files etc:
### 0.parameterListGenerator.R

## Downloads and maps samples to their respective reference genome. Outputs consensus mtDNA fasta for each individual:
### 1.DownloadMapMitochondria.slurm

## Aligns protein-coding regions across genomes. Makes GFF files for future use.
### 2.ProteinExonerateGenomes.md

## Makes gene-level fastas across all samples. Uses IQTREE to make both: 1) individual gene trees and 2) mitochondrial genome tree. Outgroup=European Daphnia magna
### 3.Concatenate_Genomes_align.sh

## Plots mtDNA tree.
### 4.AnalysisTrees.R

### Finalized tree for protein-coding mtDNA:
<img src="https://user-images.githubusercontent.com/55203772/220213780-4afc4091-430c-4515-ba87-d880e923db1a.png" width="500" height="500">

## Conversion between exonerate output and GFF
### convert_exonerate_gff_to_gff3.py

## Renames fasta headers.
### RenameIndexFastaHeaders.sh
