# Connor Murray
# 11.4.2024
# Run eSMC on Daphnia dataset

# Install
# git clone https://github.com/TPPSellinger/eSMC
# path="/home/csm6hg/eSMC/eSMC_2.0.5.tar.gz" # Path to the dowloaded eSMC package
#devtools::install_local(path)

# Library
library(eSMC)
library(tidyverse)
library(data.table)

# Working directory
setwd("/project/berglandlab/connor/backup_project/")
#doParallel::registerDoParallel(20)

# Metadata
fin <- data.table(read.csv("/project/berglandlab/connor/backup_project/metadata/samples.fin.9.8.22.csv"))

# Chromosomes
chroms <- fread("metadata/goodChrom.txt", header = F)

# HETSEP files ~14k
files <- data.table(list.files(path = "chapter1/msmc", 
                    pattern = "*.phase.filt.samp", 
                    recursive = T, 
                    full.names = T))

# Arrange metadata
files[,sample:=as.character(tstrsplit(V1, "/")[[4]])]
files[,chrom:=as.character(tstrsplit(V1, "/")[[3]])]
meta <- files %>% 
  mutate(Sample=str_remove_all(str_remove_all(sample, chrom), 
                                "..phase.filt.samp")) %>% 
  left_join(fin %>% select(Sample, Species, Continent)) %>% 
  filter(!is.na(Species)) %>% 
  mutate(id=1:nrow(.))

meta

# Output
write.table(meta[1:9000], file="../eSMC.1.list", 
            quote = F, row.names = F, col.names = F, sep = "\t")

write.table(meta[!1:9000], file="../eSMC.2.list", 
            quote = F, row.names = F, col.names = F, sep = "\t")

