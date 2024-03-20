# Representative individual per country and species
# Connor Murray 12.13.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(tidyverse)

# Working directory
setwd('/project/berglandlab/connor/msmc/')

# Metadata
fin <- data.table(fread("../metadata/samples.fin.3.9.22.csv"))

# Random individual per country and species
set.seed(100)
p <- data.table(fin[cont %like% "Daphnia.pulex."] %>%
                  group_by(mlg.country) %>% 
                  summarize(Sample = sample(Sample, 1)))

# Write genomeList
write_delim(p %>% select(Sample), delim = "\t",col_names = F,
            file ="/project/berglandlab/connor/msmc/mlg.1.species.country.list")

