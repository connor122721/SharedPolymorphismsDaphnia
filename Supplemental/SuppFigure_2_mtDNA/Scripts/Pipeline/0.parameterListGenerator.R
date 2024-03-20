# Samples for mtDNA phylogeny
# Connor Murray 12.12.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)

# Consistent seed
set.seed(100)

# Working directory
setwd("/scratch/csm6hg/")

# Metadata
fin <- data.table(read.csv("data/samples.fin.9.8.22.csv"))

# Samples in snapp tree
snapp <- data.table(fread("/project/berglandlab/connor/snapp3/individuals.2inds.continent.withObtusa.txt", header = F))

# Merge with dorthe fastq meta
h1 <- list.files("all_fqs/dorthe", 
                pattern = "dorthe", full.names = T)
dt.fast1 <- data.table(do.call(cbind, lapply(h1, fread, header=F))) 
colnames(dt.fast1) <- c("for.fasta", "rev.fasta", "Sample")

# Join with metadata
meta <- data.table(left_join(fin, dt.fast1, by=c("Sample")))

# Merge with europe fastq meta
h <- list.files("all_fqs/Euro_fqs/metadata", full.names = T)
dt.fast <- data.table(do.call(rbind, lapply(h, fread, header=F)) %>% 
                        mutate(index=tstrsplit(V1, "SL")[[1]],
                               Sample=tstrsplit(V1, "SL")[[2]],
                               SL=tstrsplit(str_extract(V1, "SL.*"), "_")[[1]]) %>% 
                        mutate(Sample=substring(Sample, 8),
                               index=paste(index, SL, sep=""))) 
colnames(dt.fast)[1] <- "fastq"

# Join with metadata
meta <- data.table(left_join(meta, dt.fast, by=c("Sample")))

# Include set of samples from Snapp phylogeny
snapple <- meta[Sample %in% snapp$V2]

# Combine samples
fin.dt <- data.table(rbind(meta, snapple))
fin.dt <- fin.dt[!duplicated(Sample)]

# Extract sample names and add position in scratch
fin.out <- data.table(fin.dt %>% 
               summarize(Sample=Sample, 
                         Species=Species,
                         Continent=Continent,
                         Origin=Origin,
                         cont=cont,
                         bam = case_when(Origin == "Europe" ~ "/scratch/csm6hg/all_bam/Euro_bams",
                                         Origin == "SRA" ~ "/scratch/csm6hg/all_bam/final_bam",
                                         Origin == "Dorthe" ~ "/scratch/csm6hg/all_bam/bam_dorthe_fin"),
                         fast.wildcard = index,
                         fast.name = fastq,
                         for.fasta, 
                         rev.fasta))

# Add Outgroup - Ref Daphnia magna
out.group <- data.table(Origin="SRA", 
                         bam="", 
                         Species="Daphnia magna", 
                         Continent="Europe",
                         Sample="Daphnia_magna")
fin.out <- data.table(fin.out %>% full_join(out.group))

# Fix slurmid column
dt <- data.table(fin.out %>% 
          mutate(slurmID=c(1:c(dim(fin.out)[1]))))

# Reorder
setcolorder(dt, c("slurmID", "Sample", "Species", 
                  "Continent", "bam", "fast.name", 
                  "fast.wildcard"))

# Output parameter file
write.csv(dt %>% select(-c(cont)), quote=F, row.names=F,
          file="/scratch/csm6hg/mito/samples_mito")
