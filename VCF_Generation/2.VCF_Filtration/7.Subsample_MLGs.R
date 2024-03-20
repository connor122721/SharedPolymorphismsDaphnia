# Subsampling by multilocus genotypes
# Connor Murray 9.8.2022

# Run using these commands:
# ijob -A berglandlab_standard --mem=10G -p standard -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2")

# Metadata
fin.dt <- data.table(read.csv(file = "../metadata/samples.9.8.22.csv", header = T))

# Set seed
set.seed(100)

# Reduce redundancy
fin <- data.table(fin.dt %>% 
                    group_by(mlg.country) %>%
                    sample_n(1))

table(table(fin$mlg.country))

# Write new metadata
#write.csv(fin, file = "../metadata/samples.fin.9.8.22.csv", quote = F, row.names = F)

# Samples/nMLGs
nSamp <- fin.dt %>% 
  group_by(cont) %>% 
  summarize(nSamp=n())

nMLG <- fin.dt %>% 
  group_by(cont) %>% 
  summarize(nMLG=length(unique(mlg.country)))

merge(nSamp,nMLG) %>% 
  group_by(cont) %>% 
  summarize(prop=nMLG/nSamp)
