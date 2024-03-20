# Explore Candidate Genes & Mutations & Exons
# 3.23.2023
# ijob -c 15 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(readxl)
library(SeqArray)
library(foreach)

# Working directory
setwd("/project/berglandlab/connor/")

# Register cores
doParallel::registerDoParallel(cores = 15)

# SNP metadata 
tot <- data.table(readRDS(file="/project/berglandlab/connor/candgene/classified_snps_filt_exon.rds"))
tot.tsp <- data.table(rbind(tot[c(Daphnia.pulex.NorthAmerica>0.1 & Daphnia.pulex.Europe>0.1)],
                            tot[Daphnia.pulex.NorthAmerica>0.1 | Daphnia.pulex.Europe>0.1][!classified=="shared_poly"]))

# Add number of snps around focal SNP
bp = 250
snp.density <- foreach(i=1:dim(tot.tsp)[1], .combine = "rbind") %dopar% {

  # Start  
  print(i)
  
  # Extract data
  focal <- tot.tsp[i]
  tot.i <- tot[chrom==focal$chrom][position %in% c((focal$position-bp):(focal$position+bp))]
  
  # Summarize
  pp <- data.table(density_snps=length(unique(tot.i$variant.id))-1,
                   density_tsp=length(unique(tot.i[classified=="shared_poly"]$variant.id))-1,
                   i=i)
  
  # Finish
  return(pp)
}

# Add density info
tot.tsp <- cbind(tot.tsp, snp.density)
#saveRDS(tot.tsp, file = "data/density.tsps.focalsnps_250bp_all.rds")

# Filter to snps without other TSPs in sequence
tot.tsp1 <- data.table(tot.tsp %>% 
                filter(snp.density>20 & snp.density<200 & density_tsp==0))

# Unique genes in both conditions
tsp.genes <- unique(tot.tsp1[classified=="shared_poly"]$gene)
control.genes <- unique(tot.tsp1[!classified=="shared_poly"]$gene)

# Independent genes
tot.tsp2 <- data.table(tot.tsp1 %>% 
          filter(!c(!classified=="shared_poly" & gene %in% tsp.genes)) %>% 
          filter(!c(classified=="shared_poly" & gene %in% control.genes)))

# Subset by 1 tree per gene
tot.tsp3 <- data.table(tot.tsp2 %>% 
              group_by(classified, gene) %>% 
              slice(1))

# Select TSPs - add bps around sequence
gene.gtf3 <- data.table(tot.tsp3 %>% 
  select(chrom, gene, position, simpleAnnot, classified) %>% 
  mutate(start = position - bp,
         end = position + bp) %>% 
  mutate(ch = paste(chrom, ":", start, "-", end, sep="")))

# Take 1000 random non-tsp
#gene.gtf3 <- gene.gtf3[!classified=="shared_poly"] %>% sample_n(1000)

# Write output
write.table(data.table(gene.gtf3$ch,
                       gene.gtf3$gene,
                       gene.gtf3$classified,
                       gene.gtf3$simpleAnnot), 
            file="candgene/exons.genome.list.tspset6", 
            quote=F, row.names=F, col.names=F, sep = ",")
