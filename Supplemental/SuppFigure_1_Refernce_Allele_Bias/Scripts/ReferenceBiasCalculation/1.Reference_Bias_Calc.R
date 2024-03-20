# Reference allele bias at heterozygous sites
# Connor Murray 6.23.2022

# Run using these commands:
# ijob -A berglandlab_standard --mem=50G -p standard -c 15
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(slider)
library(parallel)
library(doParallel)

# Cores
registerDoParallel(cores=15)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2/")

# Metadata
fin <- data.table(read.csv("../metadata/samples.fin.3.9.22.csv"))

# Genofile
out <-  c("daphnia.filt.qual.miss.rep.ann.gds")

# Open GDS
genofile <- seqOpen(out)

# BUSCO SNPs
busco <- data.table(fread("../metadata/busco_snps"))

# Subsample 10 per group
fin <- data.table(rbind(fin[Species %in% c("Daphnia pulex")] %>% 
                    group_by(Species, Continent) %>% 
                    sample_n(100),
                  fin[!Species %in% c("Daphnia pulex")]))

# Individual foreach
p.out <- foreach(j=1:length(fin$Sample), .combine = "rbind") %dopar% { 

    # Progress message
    print(paste("Iteration", j, sep=" "))
  
    # Reset any previous filters
    seqResetFilter(genofile)
    
    # Select all samples from group
    seqSetFilter(genofile, 
                 sample.id=as.vector(fin[j]$Sample),
                 variant.sel = unique(busco$variant.id))
    
    # Find heterozygous sites
    snp.dt <- data.table(variant.id = seqGetData(genofile, "variant.id"),
                         af = seqAlleleFreq(genofile),
                         num = seqGetData(genofile, "$num_allele"))
    
    # All heterozgous biallelic snps
    het.snps <- snp.dt[num==2][af==0.5]
    fix.snps <- snp.dt[num==2][af==1]
    not.snps <- snp.dt[num==2][af==0]
    
    # Reset any previous filters
    seqResetFilter(genofile)
    
    # Select all samples from group
    seqSetFilter(genofile, 
                 sample.id=as.vector(fin[j]$Sample), 
                 variant.sel = sample(c(het.snps$variant.id,
                                        fix.snps$variant.id,
                                        not.snps$variant.id), 
                                      size = 1000, 
                                      replace = F))

    # Get depth information
    ad <- seqGetData(genofile, "annotation/format/AD")
    dp <- seqGetData(genofile, "annotation/format/DP")
    
    # Alelle frequency at heterozygous sites
    snp.dt1 <- data.table(variant.id = seqGetData(genofile, "variant.id"),
                         length = seqGetData(genofile, "$num_allele"),
                         dos.alt = data.table(data.table(t(ad$data)) %>%
                                   filter(row_number() %% 2 == 1))$V1,
                         dos.ref = as.numeric(data.table(t(dp))$V1),
                         Species = fin[j]$Species,
                         Continent = fin[j]$Continent,
                         Country = fin[j]$country,
                         Sample = fin[j]$Sample,
                         iteration = j) %>% 
                    left_join(snp.dt, by="variant.id") %>% 
                    mutate(af.dos = dos.alt/dos.ref)
    
    # Summarize data
    snp.out <- data.table(snp.dt1 %>%
                            group_by(Species, Continent, 
                                     Country, Sample, 
                                     iteration, af, num) %>% 
                            summarize(avg.dos.ref = mean(dos.ref, na.rm = T)/2,
                                      avg.dos.alt = mean(dos.alt, na.rm = T)/2) %>% 
                            mutate(avg.af = avg.dos.alt/avg.dos.ref))
    
    # Finish individual
    return(snp.dt1)
}

# Save output
write.csv(p.out, file="../data/refallelebias_subsamp_busco1000snps.csv", quote = F, row.names = F)
