# Compile and output vcftools QC data 
# This script will concatenate all QC output files from chromosomes per VCF file
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(SeqArray)

# Working directory
setwd("/oldscratch/csm6hg/daphnia_phylo/vcf/raw_vcf")

# Data
dep <- data.table(read.csv("idep.vcftools.csv"))
miss <- data.table(read.csv("imiss.vcftools.csv"))

# Join
fin <- right_join(dep, miss, 
                  by = c("INDV"))

# Summarize stats
fin.map <- data.table(fin %>% 
                        group_by(INDV) %>% 
                        summarise(mean.dep = mean(MEAN_DEPTH),
                                  mean.miss = mean(F_MISS)))

# Individuals
total <- unique(fin$INDV)

# Metadata
meta <- data.table(read.csv("../../data/raw.samples.2.6.22.csv"))
meta <- meta[Sample %in% total]

# Fix species id
meta <- data.table(meta %>% 
            mutate(Species=case_when(
              is.na(Species) ~ "Confirm",
              Species == "soil metagenome" ~ "Daphnia pulex - Metagenome",
              Species == "UNK" ~ "Confirm",
              TRUE ~ as.character(Species))))

# Merge with metadata
fin <- data.table(fin.map %>% 
                left_join(meta, 
                          by = c("INDV"="Sample")) %>% 
                 pivot_longer(cols = c('mean.dep', 'mean.miss')))

#write.csv(fin[INDV %like% "index"], "/project/berglandlab/connor/metadata/dorthe.samples.wdep.csv", quote = F, row.names = F)

pdf("../raw_vcf_mean_miss_depth_jitter.pdf", width = 12, height = 8)

# Individual based missing and depth
fin %>% 
  ggplot(.,
         aes(x=value,
             y=Species)) + 
  geom_jitter(size =1.2, 
             shape = 21, 
             fill ="grey") + 
  facet_wrap(~name, scales = "free_x") +
  geom_vline(data=fin[name=="mean.miss"], 
             aes(xintercept = c(0.2)), linetype=2, size=1.3) +
  geom_vline(data=fin[name=="mean.dep"], 
             aes(xintercept = c(10)), linetype=2, size=1.3) +
  theme_bw() + 
  labs(x="Average missingness or depth",
       y="Species",
       title="Raw VCF Quality Control") +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold.italic", size = 14),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        title = element_text(face = "bold", size=18), 
        strip.text =  element_text(face = "bold", size = 20))

dev.off()

pdf("../raw_vcf_mean_miss_depth_boxplot.pdf", width = 12, height = 8)

# Individual based missing and depth
fin %>% 
  ggplot(.,
         aes(x=value,
             y=Species)) + 
  geom_boxplot(size =1.2, 
              shape = 21, 
              fill ="grey") + 
  facet_wrap(~name, scales = "free_x") +
  geom_vline(data=fin[name=="mean.miss"], 
             aes(xintercept = c(0.2)), linetype=2, size=1.3) +
  geom_vline(data=fin[name=="mean.dep"], 
             aes(xintercept = c(10)), linetype=2, size=1.3) +
  theme_bw() + 
  labs(x="Average missingness or depth",
       y="Species",
       title="Raw VCF Quality Control") +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold.italic", size = 14),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        title = element_text(face = "bold", size=18), 
        strip.text =  element_text(face = "bold", size = 20))

dev.off()

# Find D.pulex individuals over than 20% missingness
drop.miss <- unique(fin[name=="mean.miss"][value >= 0.2][Species=="Daphnia pulex"]$INDV)

# Ensure dropped samples are not Obtusa or Pulicaria
meta2 <- data.table(read.csv("/project/berglandlab/connor/metadata/samples.1.31.22.csv"))
table(drop.miss %in% meta2$Sample)
samps_to_drop <- c(drop.miss[which(drop.miss %in% meta2$Sample==F)])

# Write table
write.table(samps_to_drop, 
            file = "/project/berglandlab/connor/high_missingness_samples_raw", 
            quote = F, 
            row.names = F, 
            col.names = F)

# Samples in metadata
samps2 <- fread("/project/berglandlab/connor/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.recode.samples", header = F)
table(meta2$Sample %in% samps2$V1)

var_depth <- fread("/scratch/csm6hg/daphnia_phylo/data/depth.site.daphnia.test.ldepth.mean")
colnames(var_depth) = c("chr", "pos", "mean_depth", "var_depth")

a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill="steelblue", colour = "black", alpha = 0.3)
a + theme_light() +xlim(c(0, 40))

pdf(file = "/scratch/csm6hg/daphnia_phylo/figures/meandepthpersite_vcf.pdf")
a
dev.off()

# Upper and lower quantiles of depth
quantile(var_depth$mean_depth, probs = 0.05)
mean(var_depth$mean_depth)
quantile(var_depth$mean_depth, probs = 0.95)

# Generate bed of low and high depth regions
bed.coverage <- var_depth[mean_depth<=10 | mean_depth>=20]

# Output coverage bed
write.table(bed.coverage %>% select(chr, pos) %>% mutate(pos2=pos), 
            file = "/project/berglandlab/connor/data/meanDepth.10low.20high.bed", 
            quote = F, sep = "\t", row.names = F, col.names = F)


