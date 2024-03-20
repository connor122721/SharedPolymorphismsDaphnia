# Make Figure 2
# Connor Murray 12.9.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(colortools)
library(forcats)
library(patchwork)

# Working directory
setwd('/project/berglandlab/connor/')

# Metadata
meta <- data.table(fread("metadata/samples.fin.9.8.22.csv"))

# Admixture data
fin.dt <- data.table(read.csv(file = "data/admix.sub.data.k2_23.csv") %>% 
    mutate(Species=case_when(Species=="Daphnia pulex" & Continent=="Europe" ~ "Euro D. pulex",
                   Species=="Daphnia pulicaria" & Continent=="Europe" ~ "Euro D. pulicaria",
                   Species=="Daphnia obtusa" & Continent=="Europe" ~ "Euro D. obtusa",
                   Species=="Daphnia pulex" & Continent=="NorthAmerica" ~ "NAm. D. pulex",
                   Species=="Daphnia pulicaria" & Continent=="NorthAmerica" ~ "NAm. D. pulicaria",
                   Species=="Daphnia obtusa" & Continent=="NorthAmerica" ~ "NAm. D. obtusa",
                   Species=="Daphnia pulexcaria" ~ "NAm. D. pulex x D. pulicaria")))

# Cross validation (CV) scores
err <- data.table(fread("data/daphnia.sub.cv.error"))
colnames(err) <- c("K", "CV_error")

# K=9 has the minimal CV
err[CV_error == min(err$CV_error)]

# Create labels
lab <- data.table(lab.cont = unique(fin.dt$Species))

# Order for plotting
fin.dt$Species1 <- factor(fin.dt$Species,  
                         levels = c("Euro D. pulex", "Euro D. pulicaria", "Euro D. obtusa", 
                                    "NAm. D. pulex", "NAm. D. pulex x D. pulicaria", "NAm. D. pulicaria",
                                    "NAm. D. obtusa"))
# Plot K of various admixtures
admix <- {fin.dt[k %in% c(9)] %>% 
  mutate(kk = paste("K=", k, sep="")) %>% 
  mutate(kk=factor(kk, 
                   levels = c("K=9"))) %>% 
  ggplot(., aes(x = factor(id), y = value, fill = factor(pop))) +
  geom_col(size = 0.01) +
  facet_grid(kk~Species1, 
             switch = "x", 
             scales = "free") +
  theme_minimal() + 
  labs(x = "", 
       y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_fill_brewer(palette ="Paired") +
  theme(panel.spacing.x = unit(0.1, "lines"),
        panel.grid = element_blank(),
        strip.text = element_text(face="bold.italic", 
                                  size=12, 
                                  angle=30),
        legend.text = element_text(size=16),
        legend.position = "none",
        legend.title = element_text(face="bold", size=22),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size=22),
        axis.title.y = element_text(face="bold", size=22),
        axis.title = element_text(face="bold", size=23))}

# Read in fixed differences data
o.fin2 <- data.table(readRDS("data/format.fixed.50reps.4samps.rds"))

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

# Compile figures
mega.fig <- ( admix / shared.prop )

# Output
pdf("figures/Figure2.part.pdf", width = 12, height = 10)
mega.fig
dev.off()
