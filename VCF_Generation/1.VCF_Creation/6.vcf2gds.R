# GDS 2 VCF files
# 9.7.2022
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
out <-  c("combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")
in.vcf <- c("combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz")

# Register cores
doParallel::registerDoParallel(cores = 15)

# Convert VCF to GDS
seqVCF2GDS(in.vcf, out, parallel = 15, verbose = TRUE)
