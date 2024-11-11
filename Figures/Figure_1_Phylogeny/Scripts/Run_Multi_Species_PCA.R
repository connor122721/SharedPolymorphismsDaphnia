# Make Principal component analysis of WGS data
# 6.24.2023
# ijob -c 10 --mem=40G -p standard -A berglandlab_standard
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(doParallel)
library(SNPRelate)
library(ggforce)

# Working directory
setwd("/project/berglandlab/connor/metadata")
doParallel::registerDoParallel(cores = 10)

# SNP GDS
out <-  c("../new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Open GDS
genofile <- seqOpen(out)

# Samples full
samps <- seqGetData(genofile, var.name = "sample.id")

# Sample metadata
fin <- data.table(read.csv("../metadata/samples.fin.9.8.22.csv"))

# SNP metadata
snps <- data.table(fread("../metadata/snps_new"))

# Run PCA
ccm_pca <- snpgdsPCA(genofile, 
                     autosome.only = F, 
                     sample.id = as.character(fin$Sample),
                     num.thread = 10, 
                     snp.id = unique(snps$variant.id),
                     maf=0.01)

# Output PCA
saveRDS(ccm_pca, file = "../new_vcf2/pca.daphnia.genome.filt.mlgsub.outgroup.rds")

# Skree plot
plot(ccm_pca$varprop[1:20]*100)

# Merge by Sample
pca <- data.table(sample=ccm_pca$sample.id, PC1=ccm_pca$eigenvect[,1],
                  PC2=ccm_pca$eigenvect[,2], PC3=ccm_pca$eigenvect[,3],
                  PC4=ccm_pca$eigenvect[,4], PC5=ccm_pca$eigenvect[,5])

# Remove missing samples
pca <- data.table(pca[sample %in% fin$Sample])

# Merge PCA and metadata
pca <- data.table(merge(pca, fin,
                        by.x="sample", by.y="Sample"))

# PC 1/2 Plot
pc <- {
  
  pca %>%
  ggplot(., aes(x=PC1, 
                y=PC2, 
                fill=Continent,
                shape=Species)) +
  #facet_wrap(~Species+Continent) +
  geom_point(size=6, alpha=0.6) +
  geom_mark_ellipse() +
  theme_minimal() + 
  labs(x=paste("PC1 (", round(ccm_pca$varprop[[1]], digits=3)*100, " %)", sep=""),
       y=paste("PC2 (", round(ccm_pca$varprop[[2]], digits=3)*100, " %)", sep="")) +
  scale_fill_brewer(name = "Species complex", palette = "Set1") +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.text = element_text(size=16, face="bold.italic"),
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))
}

# Output plot
ggsave(pc, "../figures/pca12.filt.qual.miss.rep.dep.chr.ann.busco.nooutgroup.pdf")
