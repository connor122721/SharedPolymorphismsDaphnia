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
# saveRDS(dt, "/project/berglandlab/connor/fixed/fixed.50reps.4samps.rds")
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
                              spp_comp == "pulicaria_pulicaria" ~ "D.pulicaria")))

# Fixed and polymorphic SNPs
dt2 <- data.table(o.fin2 %>% 
  group_by(cont_comp, spp_comp, Species, Dosage) %>% 
  summarize(mean.N=mean(N, na.rm = T),
            mean.total_shared_snps=mean(mean.total_shared_snps, na.rm = T)) %>% 
  mutate(cont_comp=factor(cont_comp,
         levels = c("Europe","Between Continents","NorthAmerica"))))

# Shared heterozygous alleles - num nam pulic pulex
shared.n <- {ggplot(o.fin2[spp_comp.tot=="Daphnia.pulicaria.NorthAmerica_Daphnia.pulex.NorthAmerica"][Spp %in% 
                      c("Daphnia pulex.Europe", "Daphnia pulicaria.Europe", 
                        "Daphnia obtusa.Europe", "Daphnia obtusa.NorthAmerica",
                        "Daphnia pulexcaria.NorthAmerica")][Dosage==1], 
       aes(x=N,
           y=Spp,
           group=as.factor(Dosage))) +
  geom_boxplot() +
  scale_x_log10() +
  annotation_logticks(scaled = TRUE, sides = "b") +
  facet_wrap(~Spp, ncol = 1, scales = "free_y") +
  theme_minimal() +
  labs(x="Mean number of heterozygous SNPs", 
       y="",
       title="Cross: NAm D.pulicaria x NAm D.pulex",
       fill="Dosage of Shared Polymorphism") +
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

# Shared heterozygous alleles - prop nam pulic pulex
shared.prop <- {ggplot(o.fin2[spp_comp.tot=="Daphnia.pulicaria.NorthAmerica_Daphnia.pulex.NorthAmerica"][Spp %in% 
                              c("Daphnia pulex.Europe", "Daphnia pulicaria.Europe", 
                                "Daphnia obtusa.Europe", "Daphnia obtusa.NorthAmerica",
                                "Daphnia pulexcaria.NorthAmerica")][Dosage==1], 
       aes(x=N/mean.total_fixed_snps,
           y=Spp,
           group=as.factor(Dosage))) +
  geom_boxplot() +
  scale_x_log10() +
  annotation_logticks(scaled = TRUE, sides = "b") +
  facet_wrap(~Spp, ncol = 1, scales = "free_y") +
  theme_minimal() +
  labs(x="Proportion of heterozygous SNPs", 
       y="",
       title="Cross: NAm D.pulicaria x NAm D.pulex",
       fill="Dosage of Shared Polymorphism") +
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

# Shared heterozygous alleles - num euro nam pulex
num.euro <- {ggplot(o.fin2[spp_comp.tot=="Daphnia.pulex.NorthAmerica_Daphnia.pulex.Europe"][Spp %in% 
                                         c("Daphnia pulicaria.Europe", "Daphnia pulicaria.NorthAmerica",
                                           "Daphnia obtusa.Europe", "Daphnia obtusa.NorthAmerica",
                                           "Daphnia pulexcaria.NorthAmerica")][Dosage==1], 
                     aes(x=N,
                         y=Spp,
                         group=as.factor(Dosage))) +
    geom_boxplot() +
    scale_x_log10() +
    annotation_logticks(scaled = TRUE, sides = "b") +
    facet_wrap(~Spp, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(x="Mean number of heterozygous SNPs", 
         y="",
         title="Cross: NAm D.pulex x Euro D.pulex",
         fill="Dosage of Shared Polymorphism") +
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

# Shared heterozygous alleles - prop euro nam pulex
prop.euro <- {ggplot(o.fin2[spp_comp.tot=="Daphnia.pulex.NorthAmerica_Daphnia.pulex.Europe"][Spp %in% 
                                        c("Daphnia pulicaria.Europe", "Daphnia pulicaria.NorthAmerica",
                                          "Daphnia obtusa.Europe", "Daphnia obtusa.NorthAmerica",
                                          "Daphnia pulexcaria.NorthAmerica")][Dosage==1], 
                     aes(x=N/mean.total_fixed_snps,
                         y=Spp,
                         group=as.factor(Dosage))) +
    geom_boxplot() +
    scale_x_log10() +
    annotation_logticks(scaled = TRUE, sides = "b") +
    facet_wrap(~Spp, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(x="Proportion of heterozygous SNPs", 
         y="",
         title="Cross: NAm D.pulex x Euro D.pulex",
         fill="Dosage of Shared Polymorphism") +
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

# Output figures
pdf("figures/nam_pulex_pulic.pdf", width = 12, height = 6)
shared.n
dev.off()

pdf("figures/nam_pulex_pulic_prop.pdf", width = 12, height = 6)
shared.prop
dev.off()

pdf("figures/nam_pulex_euro_pulex.pdf", width = 12, height = 6)
num.euro
dev.off()

pdf("figures/nam_pulex_euro_pulex_prop.pdf", width = 12, height = 6)
prop.euro
dev.off()
