# Species for ADMIXTURE
# Connor Murray 6.20.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(tidyverse)
library(foreach)

# Working directory
setwd('/project/berglandlab/connor/admix/new')

# Metadata
fin <- data.table(fread("../../metadata/samples.fin.9.8.22.csv"))

# Exclude outgroup (D.obtusa)
#fin <- data.table(fin[!Species=="Daphnia obtusa"])

# Downsample D.pulex 
set.seed(100)
pulexNam <- sample(fin[cont %in% c("Daphnia.pulex.NorthAmerica")]$Sample, 50)
pulexEuro <- sample(fin[cont %in% c("Daphnia.pulex.Europe")]$Sample, 50)

# Downsample
fin.dt <- data.table(rbind(fin[Sample %in% c(pulexNam, pulexEuro)], 
                           fin[!Species %in% c("Daphnia pulex")]))

# Unique subspecies list with poppr clones
write.table(fin.dt %>% 
            select(Sample, cont),
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            file ="/scratch/csm6hg/daphnia_phylo/admixture/popfile.cont.sub.clust")

# Write genome comparisions list
write.table(fin.dt %>% 
            select(Sample),
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            file ="/scratch/csm6hg/daphnia_phylo/admixture/pop.clust.sub")

### Create independent data for Pulex

# Working directory
setwd('/project/berglandlab/connor/admix/independent')

# Metadata
fin <- data.table(fread("../../metadata/samples.fin.3.9.22.csv"))

# Unique subspecies list
dt <- foreach(i=1:length(unique(fin$cont)), .combine = "rbind") %do% {
  
  # Print progress 
  print(unique(fin$cont)[i])
  
  # Write genome comparisions list
  write.table(fin[cont==unique(fin$cont)[i]] %>% 
                select(Sample),
              sep = "\t",
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE,
              file =paste("/scratch/csm6hg/daphnia_phylo/admixture/popfile", unique(fin$cont)[i], sep="."))
  
}
