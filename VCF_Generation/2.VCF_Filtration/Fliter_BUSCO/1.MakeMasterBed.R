# Compile list of sites to keep
# 9.8.2022
# ijob -c 1 --mem=40G -p largemem -A berglandlab 
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(SeqArray)
library(tidyverse)

# Master bed file
bed <- data.table(fread("/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed"))
colnames(bed) <- c("chrom", "start", "stop")

# Good chromosomes
chrom1 <- fread("/project/berglandlab/connor/metadata/goodChrom.txt", header = F)

# Chromosomes of interest
bed <- bed[chrom %in% unique(chrom1$V1)]

# Go through each chromosome in bed
bed.po <- foreach(i=1:length(unique(bed$chrom)), .combine = "rbind", .errorhandling = "remove") %do% {
  
  # Start
  print(i)
  
  # Chromosome of interest
  chrom.i <- unique(bed$chrom)[i]
  tmp <- bed[chrom==chrom.i]
  
  # Extend each region
  tmp.df <- foreach(j=1:dim(tmp)[1], .combine = "rbind") %do% {
    
    h.h <- data.table(chrom=chrom.i,
                      position=unique(tmp[j]$start:tmp[j]$stop)) %>% 
      mutate(ch=paste(chrom, position, sep = "_"))
    
    return(h.h)
  }
 
   return(tmp.df)
}

# Genofile
out <- c("/project/berglandlab/connor/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Open GDS
genofile <- seqOpen(out)

# Create gds metadata
snp.dt.fin <- data.table(chrom=seqGetData(genofile, var.name = "chromosome"),
                         variant.id=seqGetData(genofile, var.name = "variant.id"),
                         position=seqGetData(genofile, var.name = "position")) %>% 
  mutate(ch=paste(chrom, position, sep = "_"))

# Good positions
snp.dt.fin1 <- snp.dt.fin[!ch %in% unique(bed.po$ch)]

# Output filtered sites
saveRDS(snp.dt.fin1, 
        file = "/project/berglandlab/connor/data/filtered.dep.miss.rep.chr.strict.good.sites")

# Output long format bed
saveRDS(bed.po %>% 
          mutate(ch=paste(chrom,position,sep = "_")), 
        file = "/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.strict.final.bad")
