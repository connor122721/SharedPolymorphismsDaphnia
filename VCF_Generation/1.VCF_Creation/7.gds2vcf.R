# GDS 2 VCF files
# 11.18.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(doParallel)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2")

# Executable in command line
in.gds <-  c("combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")
out <- c("daphnia.filt.mlg.genome.11.18.22.vcf")

# Metadata
fin <- data.table(read.csv("../metadata/samples.fin.9.8.22.csv"))

# Read in total SNPs
tot_snps <- data.table(fread("../metadata/snps_new"))

# Load GDS
genofile <- seqOpen(in.gds)

# Filter GDS
seqResetFilter(genofile)
seqSetFilter(genofile,
             sample.id=unique(fin$Sample),
             variant.sel=unique(tot_snps$variant.id))

# Convert VCF to GDS
seqGDS2VCF(gdsfile = genofile, vcf.fn = out, verbose = TRUE)
