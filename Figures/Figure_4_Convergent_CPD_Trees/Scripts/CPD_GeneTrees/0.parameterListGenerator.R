# Representative sample gene TSP tree
# Connor Murray 2.22.23
# ijob -A berglandlab --mem=20G -p standard -c 10
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)

# Consistent seed
set.seed(100)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Subsample - 30 per continent
ccc <- data.table(fin %>% 
                    group_by(cont) %>%
                    sample_n(30) %>% 
                    select(Sample, Continent))

# Include mainland european samples
ccc <- rbind(ccc,
             data.table(samps[country %in% c("Daphnia.pulex.Europe.LTU",
                                             "Daphnia.pulex.Europe.CZE")] %>% 
                          select(cont, Sample, Continent)))

# Output sample list
write.table(ccc, 
            file="/project/berglandlab/connor/candgene/samples_phase_var1_meta", 
            quote=F, row.names=F, col.names=F, sep="\t")

write.table(ccc %>% select(Sample, Continent), 
            file="/project/berglandlab/connor/candgene/samples_phase_var1", 
            quote=F, row.names=F, col.names=F, sep="\t")

