# Shared polymorphism for VCF files
# 11.2.2022
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(viridis)
library(doParallel)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# output files fst
files <- system("ls -f -R /project/berglandlab/connor/fixed/data/*.rds", intern = TRUE)
list <- lapply(files, readRDS)
setattr(list, 'names', system("ls -f -R /project/berglandlab/connor/fixed/data/*.rds", intern = TRUE))

# Bind list 
dt <- data.table(rbindlist(list, use.names = T, idcol = T))

# Read in data
o.fin <- data.table(readRDS(file="/project/berglandlab/connor/fixed/fixed.50reps.4samps.rds"))

# Classify 
o.fin2 <- data.table(o.fin %>%
  mutate(cont.x=tstrsplit(Spp1_comp, ".", fixed=T)[[3]],
         cont.y=tstrsplit(Spp2_comp, ".", fixed=T)[[3]]) %>% 
  mutate(Spp=paste(Species,Continent, sep="."),
         Species.x=tstrsplit(Spp1_comp, ".", fixed=T)[[2]],
         Species.y=tstrsplit(Spp2_comp, ".", fixed=T)[[2]]) %>% 
  mutate(spp_comp=paste(Species.x, Species.y, sep="_"),
         spp_comp.tot=paste(Spp1_comp,Spp2_comp, sep="_"),
         cont_comp=paste(cont.x, cont.y, sep="_")) %>% 
  mutate(spp_comp = case_when(spp_comp == 'pulex_pulicaria'~"pulicaria_pulex",
                              spp_comp == 'pulex_pulexcaria' ~ "pulexcaria_pulex",
                              TRUE ~ as.character(spp_comp)),
         cont_comp = case_when(cont_comp == 'Europe_NorthAmerica' ~ 'NorthAmerica_Europe',
                               TRUE ~ as.character(cont_comp))) %>% 
  mutate(cont_comp = case_when(cont_comp == 'NorthAmerica_Europe' ~ "Between Continents",
                               cont_comp == 'Europe_Europe' ~ "Europe",
                               cont_comp == 'NorthAmerica_NorthAmerica' ~ "NorthAmerica"),
         spp_comp = case_when(spp_comp == "pulicaria_pulex" ~ "D.pulicaria x D.pulex",
                              spp_comp == "pulexcaria_pulex" ~ "D.pulexcaria x D.pulex",
                              spp_comp == "pulex_pulex" ~ "D.pulex",
                              spp_comp == "pulexcaria_pulicaria" ~ "D.pulexcaria x D.pulicaria",
                              spp_comp == "pulexcaria_pulexcaria" ~ "D.pulexcaria",
                              spp_comp == "pulicaria_pulicaria" ~ "D.pulicaria"),
         Spp=case_when(Spp=="Daphnia pulex.Europe" ~ "Euro D. pulex",
                       Spp=="Daphnia pulicaria.Europe" ~ "Euro D. pulicaria",
                       Spp=="Daphnia obtusa.Europe" ~ "Euro D. obtusa",
                       Spp=="Daphnia pulex.NorthAmerica" ~ "NAm. D. pulex",
                       Spp=="Daphnia pulicaria.NorthAmerica" ~ "NAm. D. pulicaria",
                       Spp=="Daphnia obtusa.NorthAmerica" ~ "NAm. D. obtusa",
                       Spp=="Daphnia pulexcaria.NorthAmerica" ~ "NAm. D. pulex x D. pulicaria")))

# saveRDS(o.fin2, "/project/berglandlab/connor/data/format.fixed.50reps.4samps.rds")
o.fin2 <- data.table(readRDS("/project/berglandlab/connor/fixed/format.fixed.50reps.4samps.rds"))

# Order for plotting
o.fin2$Spp1 <- factor(o.fin2$Spp,  
      levels = c("Euro D. obtusa", "NAm. D. obtusa", 
                 "Euro D. pulex", "Euro D. pulicaria", 
                 "NAm. D. pulex x D. pulicaria"))

# Shared heterozygous alleles - prop nam pulic pulex
shared.prop <- {o.fin2[spp_comp.tot=="Daphnia.pulicaria.NorthAmerica_Daphnia.pulex.NorthAmerica"][!Spp %in% 
         c("NAm. D. pulex", "NAm. D. pulicaria")][Dosage==1] %>% 
  ggplot(., aes(x=N/mean.total_fixed_snps,
           y=Spp,
           group=as.factor(Dosage))) +
  geom_boxplot(fill="cyan4") +
  scale_x_log10() +
  annotation_logticks(scaled = TRUE, sides = "b") +
  facet_wrap(~Spp1, ncol = 1, scales = "free_y") +
  theme_minimal() +
  labs(x="Proportion of heterozygous SNPs", 
       y="",
       title="Fixed differences between: NAm. D.pulicaria x NAm D.pulex") +
  theme_bw() +
  theme(legend.text = element_text(face="bold", size=14),
        title = element_text(face="bold.italic", size = 18),
        legend.title = element_text(face="bold", size=16),
        legend.background = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.text.y = element_text(face="bold.italic", size=14),
        axis.text.x = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))}
