# For treemix
# Connor Murray 1.31.2021
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(tidyverse)
library(foreach)
library(slider)
library(GenomicRanges)

# Working directory
setwd('/project/berglandlab/connor/dsuite4/')

# Metadata
fin <- data.table(fread("../metadata/samples.fin.3.9.22.csv"))
fin2 <- data.table(fread("../new_vcf2/daphnia.filt.qual.miss.rep.ann.samps", header=F))
fin2 <- data.table(fin2 %>% select(Sample=V1) %>% filter(!Sample %in% fin$Sample), cont1="xxx")

out <- data.table(fin %>% 
           mutate(cont1=case_when(cont=="Daphnia.obtusa.Europe" ~ "Outgroup",
                                  cont=="Daphnia.pulexcaria.NorthAmerica" ~ "xxx",
                                  cont=="Daphnia.pulicaria.Europe" ~ "xxx",
                                  cont=="Daphnia.obtusa.NorthAmerica" ~ "xxx",
                                  TRUE ~ as.character(cont))) %>% 
           select(Sample, cont1)) 

out2 <- data.table(rbind(out, fin2))

# Write cluster file
write.table(out2,
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            file = "pop.dsuite.country.europul.clust")

# Read in chromosome metadata
chr <- data.table(fread("../metadata/total_snps"))
chr <- data.table(chr %>% summarize(length=n()))

# Parameters
step.bp <- 500
window.bp <- 1000

# Windows
ww <- data.table(start=seq(from=1, to=chr$length, by=step.bp),
                 stop=seq(from=1, to=chr$length, by=step.bp) + window.bp) 
ww[stop > chr$length]$stop <- chr$length 
ww <- data.table(ww %>% mutate(len=stop-start))

# Split sliding window
windows <- data.table(id=c(1:dim(ww)[1]), ww)

# Write output
write.csv(windows, quote=F, row.names=F, file='interval_dsuite_paramlist_1k')
