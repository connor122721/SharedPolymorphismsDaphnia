# Representative for FST and Pi
# Connor Murray 10.6.2022
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

# Output parameter file
write.table(fin[cont=="Daphnia.pulex.Europe"][pondID=='D8'] %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.euro.D8.pop.fst.txt")

# Output parameter file

write.table(fin[cont=="Daphnia.pulex.NorthAmerica"][Sample.Name=="TEX population data"] %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.nam.TEX.pop.fst.txt")

write.table(fin[cont=="Daphnia.pulex.NorthAmerica"][isolate.name=="KAP"] %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.nam.kap.pop.fst.txt")

write.table(fin[cont=="Daphnia.pulex.NorthAmerica"][isolate.name=="KAP"] %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.nam.kap.pop.fst.txt")

# Output parameter file
write.table(fin[Species %in% "Daphnia pulex"] %>% 
            select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
          file="new_vcf2/species.pop.fst.txt")

write.table(fin[cont %in% "Daphnia.pulex.Europe"] %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.pulex.euro.pop.fst.txt")

write.table(fin[cont %in% "Daphnia.pulex.NorthAmerica"] %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.pulex.nam.pop.fst.txt")

# Random pops
write.table(fin[cont %in% "Daphnia.pulex.NorthAmerica"][country=="Daphnia.pulex.NorthAmerica.CAN"] %>% 
              sample_n(150) %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.pulex.nam.150samps_can.pop.fst.txt")

# Random pops
write.table(fin[cont %in% "Daphnia.pulex.NorthAmerica"][country=="Daphnia.pulex.NorthAmerica.USA"]  %>% 
              sample_n(150) %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.pulex.nam.150samps_usa.pop.fst.txt")

# Random pops
write.table(fin[cont %in% "Daphnia.pulex.Europe"][country=="Daphnia.pulex.Europe.GBR"][pondID %in% c("D8", "DBunk")] %>% 
              sample_n(150) %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.pulex.euro.150samps_dbunk_d8.pop.fst.txt")

# Random pops
write.table(fin[cont %in% "Daphnia.pulex.Europe"][country=="Daphnia.pulex.Europe.GBR"][!pondID %in% c("D8", "DBunk")]  %>% 
              select(Sample, cont), quote=F, col.names = F, row.names=F, sep = "\t",
            file="new_vcf2/species.pulex.euro.85samps_gbr.pop.fst.txt")