# Convert a SeqArray GDS 2 Genlight for use in adegenet
# Connor Murray 9.8.2022
# ijob -A berglandlab --mem=50G -p standard -c 10
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(ape)
library(adegenet)
library(radiator)
library(vcfR)
library(parallel)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2")

# Register cores
doParallel::registerDoParallel(cores = 10)

# All VCFs
file <- as.character(list.files(pattern = "daphnia.filtered.chr.busco.vcf.gz$"))

# Forloop write genlight
foreach(i=1:length(file)) %dopar% {
  
  print(paste("Start:", i))

  # VCF
  vcf <- read.vcfR(file[i], verbose = T)

  # Conversion VCF 2 Genlight
  x <- vcfR2genlight(vcf, n.cores = 10)

  # Save genlight object
  saveRDS(x, file = paste(str_remove(file[i], "vcf.gz"), 
                          "genlight.rds", sep=""))
  
  print(paste("Finish:", i))

}