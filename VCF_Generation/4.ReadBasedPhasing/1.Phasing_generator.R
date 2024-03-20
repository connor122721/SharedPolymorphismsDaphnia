# Create parameter files for whatshap
# Connor Murray
# 11.7.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Library
library(data.table)
library(slider)
library(GenomicRanges)
library(SeqArray)
library(tidyverse)
library(geosphere)
library(foreach)

setwd("/project/berglandlab/connor")

# Output file
out.name = "new_phase/phasing_paramList_phylo"

# All sample neames
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
out <- c("new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samp <- seqGetData(genofile, var.name = "sample.id")
samps <- fin[Sample %in% samp]$Sample

# Chromosomes
chrom = as.character(fread("metadata/goodChrom.txt", header = F)$V1)
dt <- as.data.table(expand.grid(samps, chrom))
setnames(dt, names(dt), c("samp", "chrom"))

# Failed jobs - extract ids to rerun
#fail <- as.character(str_remove(str_remove_all(fread("Scaffold_9199_HRSCAF_10755.list", header = F)$V1, 
#              pattern=c("/scratch/csm6hg/daphnia_phylo/vcf/phasing/Scaffold_9199_HRSCAF_10755/")), 
#              pattern=".Scaffold_9199_HRSCAF_10755.phase.vcf.gz"))

# Write parameter list for fin
write.csv(dt %>% 
            mutate(id=1:nrow(dt[samp %in% samp]), test="phase") %>% 
            select(id, samp, chrom, test), 
          quote = F, row.names = F, file = out.name)

# Write parameter list for fin
write.csv(dt[!1:9999] %>% 
            mutate(id=1:nrow(dt[samp %in% samp][!1:9999]), test="phase") %>% 
            select(id, samp, chrom, test), 
          quote = F, row.names = F, file = "new_phase/phasing_paramList_phylo2")

# Write parameter list for rerunning all Dorthe's samples
write.csv(dt[samp %like% "-"] %>% 
            mutate(id=1:nrow(dt[samp %like% "-"]), test="phase") %>% 
            select(id, samp, chrom, test), 
          quote = F, row.names = F, file = "new_phase/phasing_paramList_phylo3")

# Chromosome
foreach(k=1:length(chrom)) %do% {

  # Progress message
  print(paste("Chromosome:", chrom[k]))
  
  # Write list of samples
  foreach(i=1:length(unique(fin$cont))) %do% {
 
    # Progress message
    print(paste("Continent", unique(fin$cont)[i]))
    
    # Expand samples and chromosome
    dt1 <- as.data.table(expand.grid(fin[cont==unique(fin$cont)[i]]$Sample, chrom[k]))
    setnames(dt1, names(dt1), c("samp", "chrom"))
    
    # Add vcf names
    dt2 <- data.table(dt1 %>% 
      mutate(samp.vcf = paste("/scratch/csm6hg/daphnia_phylo/vcf/phasing2/",
                              "/", chrom, "/", samp, ".", chrom, ".phase.vcf.gz", sep="")) %>% 
      mutate(samp.vcf = str_replace_all(samp.vcf, pattern = "-", replacement = "_")))
    
    # Continent & Species list for phasing
    write.table(dt2 %>% 
       select(samp.vcf),
       col.names = F, 
       quote = F, 
       row.names = F, 
       file = paste("phasing2/", unique(fin$cont)[i], ".",chrom[k],".phasingList", sep=""))
    }
}

# Generate sample list for haplotype extraction 
set.seed(100)
p.dt <- data.table(fin[cont %in% c("Daphnia.pulex.Europe",
                                   "Daphnia.pulex.NorthAmerica",
                                   "Daphnia.obtusa.NorthAmerica")][Sample %in% samps] %>%
                     select(Sample, cont))

# Write genome comparisions list
write.table(p.dt %>% select(Sample, cont),
            sep = " ",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            file ="new_phase/inds.phased.species.out.list")
