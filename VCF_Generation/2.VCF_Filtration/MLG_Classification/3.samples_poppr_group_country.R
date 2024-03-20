# Classifying clonal lineages
# Connor Murray 9.8.2022

# Run using these commands:
# ijob -A berglandlab_standard --mem=50G -p standard -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(ggrepel)
library(poppr)
library(vcfR)
library(vegan)
library(adegenet)
library(viridis)
library(cowplot)
library(ggrepel)
library(parallel)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2")

# Executable in command line
arg <- commandArgs(TRUE)

# Execute in parallel
group <- as.character(arg[1])
cores <- as.character(arg[2])
#group="Daphnia.pulex.Europe";cores=10

# Metadata
fin.dt <- data.table(read.csv(file = "../metadata/samples.1.31.22.csv", header = T))

# Read in genlight object
xsnp <- readRDS("daphnia.filtered.chr.busco.genlight.snpclone.rds")

# Individuals within group and species
samps.i <- fin.dt[cont %in% group]$Sample

# Progress message
print(paste("Starting group:", group, sep=" "))

# Restrict to samples in metadata
xsnpi <- xsnp[c(which(xsnp$ind.names %in% samps.i == TRUE))]

# Make image
png(paste("filter_stats", 
          group, "png", sep="."))

# Gather thresholds to test
thresholds <- filter_stats(xsnpi,
                           plot = TRUE,  
                           distance = bitwise.dist, 
                           stats = "thresholds", 
                           threads = cores,
                           threshold = 1e+06 + .Machine$double.eps^0.5)

# Output graph
dev.off()

# Save output
saveRDS(object = thresholds, 
        file = paste("thresholds.busco.chr.MLG", 
                     group, "rds", sep="."))

# We can use these thresholds to find an appropriate cutoff
thresh <- cutoff_predictor(thresholds$farthest)
print(thresh <- cutoff_predictor(thresholds$farthest))

# Reset mlg filter
mll(xsnpi) <- "original" 

# Apply filter - for best threshold
mlg.filter(xsnpi, distance = bitwise.dist) <- thresh

# Save output
saveRDS(object = xsnpi, 
        file = paste("thresholds.busco.chr.MLG.bestthresh.snpclone", 
                     group, "rds", sep="."))

# Extract individuals and MLGs
snp.dt <- data.table(data.table(Sample=indNames(xsnpi),
                                mlg=mll(xsnpi)) %>% 
                       left_join(fin.dt %>% select(Sample, Species, Continent, SC), 
                                 by=c("Sample")))

# Add populations to snpclone object
pop(xsnpi) <- paste(str_replace(snp.dt$Species,pattern = " ", replacement = "."), 
                    snp.dt$Continent, sep=".")

# Exract mlgs
dt <- data.table(n=t(mlg.table(xsnpi, plot=F)),
                 mlg=colnames(mlg.table(xsnpi, plot=F))) 

# Wide to long reformat
dt.fin <- data.table(dt %>% 
                       pivot_longer(cols = colnames(dt)[1], 
                                    names_prefix = "n.", 
                                    names_to="pop", 
                                    values_to="n"),
                     thresh.best=thresh)

# Save output
saveRDS(dt.fin, 
        file=paste("poppr.busco.chr", 
                   group, "rds", sep="."))


# Save thresholds images
png(paste("poppr.busco.chr.thresholds", group, "png", sep="."))

# Test all thresholds
filter_stats(xsnpi, 
             distance = bitwise.dist,
             stats = "all",
             plot = TRUE,
             missing = "ignore")

dev.off()

# Distance matrix
xsnpi.dist <- bitwise.dist(xsnp)
set.seed(120)

# Graph it.
snpi <- poppr.msn(xsnpi, xsnpi.dist, gadj = 15, vertex.label = NA, showplot = F)
png(paste("poppr.busco.chr.thresholds.msn", group, "png", sep="."))
plot_poppr_msn(xsnpi, poppr_msn = snpi, inds = "none")
dev.off()