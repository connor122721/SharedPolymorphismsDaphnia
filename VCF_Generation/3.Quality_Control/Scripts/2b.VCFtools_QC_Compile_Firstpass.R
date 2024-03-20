# Compile and output vcftools QC data 
# This script will concatenate all QC output files from chromosomes per VCF file
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)

# vcftools QC output
setwd('/scratch/csm6hg/daphnia_phylo/vcf/QC/')

# Extract all vcftools files
imiss <- list.files(pattern = "*.imiss$")
idep <- list.files(pattern = "*.idepth$")
dep <- list.files(pattern = "*.ldepth.mean$")
lmiss <- list.files(pattern = "*.lmiss$")
lqual <- list.files(pattern = "*.lqual$")

# QC files
files <- list(imiss, idep, dep, lmiss, lqual)
nam <- c("imiss", "idep", "dep", "lmiss", "lqual")

# Nested loop to output files
foreach(p=1:length(files)) %do% {

  # Extract grouping
  file <- files[[p]]
  print(paste("Percentage complete:", 
              round(p/length(files)*100, 
                    digits = 2), "%"))

  # Read in files as table
  out <- foreach(i=1:length(file), .combine="rbind", .errorhandling = "pass") %do% {

    # Read file
    X <- fread(file=file[i], 
               sep="\t")

    # Final Output
    fin <- data.table(X, chrom=tstrsplit(file[[i]], "\\.")[[1]], 
                      file=file[i])

  }

# Output compiled file
write.csv(x = out, file = paste("/scratch/csm6hg/daphnia_phylo/vcf/raw_vcf/", 
                                nam[p], 
                                ".vcftools.csv", sep = ""),
          row.names = F, quote = F)

}

