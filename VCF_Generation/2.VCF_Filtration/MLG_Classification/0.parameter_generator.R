# Create parameter file for poppr
# Connor Murray 3.3.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Library
library(data.table)
library(tidyverse)

# Working directory
setwd("/project/berglandlab/connor/poppr")

# Metadata
fin.dt <- data.table(read.csv(file = "../metadata/samples.1.31.22.csv", header = T))

# Expanding parameters
dt <- as.data.table(unique(fin.dt$cont))
setnames(dt, names(dt))

# Create variables
dt <- data.table(dt %>% 
                   mutate(slurmID=c(1:c(dim(dt)[1]))))

# Reorder
setcolorder(dt, c("slurmID"))

# Output parameter file
write.csv(dt, quote=F, row.names=F,
         file="model_paramList2")
