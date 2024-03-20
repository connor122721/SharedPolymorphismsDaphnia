# Trans species polymorphism - FIS
# 3.12.2024
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(ggplot2)
library(tidyverse)
library(patchwork)

# Gene FIS averages
gene.ag <- data.table(readRDS("/scratch/csm6hg/data/Crosses_geneclass_conservative_100boot_raw.rds"))

# Blue opsin wavelength gene FIS
gene.ag[perm==0][gene %in% "Daphnia11806-RA"]

# Wild by Lab FIS
gene.ag <- data.table(gene.ag[N>2] %>% 
  mutate(exp1=case_when(exp%like%"Lab"~"Lab",
                        TRUE~"Wild"),
         sim=case_when(perm==0~"Empirical",
                       TRUE~"Simulation w/RD HWE")))

# Reorder
gene.ag$sim <- factor(gene.ag$sim, levels = c("Simulation w/RD HWE", "Empirical"))  

# Plot FIS distribution
fis.plot <- {
  gene.ag %>% 
    ggplot() +
    geom_density(aes(x = FIS, fill = sim), alpha=0.7) +
    geom_vline(data=gene.ag[sim=="Empirical"][gene %in% "Daphnia11806-RA"], 
               aes(xintercept = FIS), size= 0.8,
               color="steelblue2") +
    scale_fill_manual(values = c("Empirical"="blue",
                                 "Simulation w/RD HWE"="rosybrown")) +
    facet_wrap(~exp, ncol = 1) +
    theme_bw() +
    labs(title = "", 
         y="Density of Genes",
         fill = "",
         x = expression(italic(italic("F")["IS"]))) +
    theme(legend.position = "bottom",
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=18), 
          strip.text = element_text(face="bold", size=18), 
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))
}

# Save output
ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_Dist_Sim.pdf", plot = fis.plot, width=6, height=10, dpi=300)

### F1 frequency for genotypes ###

# Data for blue opsin gene
dt1 <- data.table(readRDS("/scratch/csm6hg/data/AxC_lab_tot_het_withmet.rds"), exp="Lab AxC F1s")
dt2 <- data.table(readRDS("/scratch/csm6hg/data/AxC_Nolab_subsamp_tot_het_withmet.rds"), exp="Wild AxC F1s Subsampled")
dt3 <- data.table(readRDS("/scratch/csm6hg/data/AxC_Nolab_tot_het_withmet.rds"), exp="Wild AxC F1s")
dt4 <- data.table(readRDS("/scratch/csm6hg/data/CxC_tot_het_withmet.rds"), exp="Lab Selfed C F1s")
dat <- data.table(rbind(dt1, dt2, dt3, dt4))
dat[all==2,geno:="hom_alt"]
dat[all==0,geno:="hom_ref"]
dat[all==1,geno:="het"]
dat[af.dos>.01 & af.dos<.99, geno:="het"]
dat[,rd:=(dos.ref + dos.alt)]
dat[,posBin:=round(position, -4)]

# Denominator for each exp
denom <- data.table(dat[all %in% c(0,1,2)][!is.na(variant.id.y)] %>% 
            group_by(exp) %>%       
            summarize(prop.d=length(unique((Sample)))))

# Distribution of F1s for genome
toti <- data.table(dat[all %in% c(0,1,2)][!is.na(variant.id.y)] %>% 
           left_join(denom) %>% 
           group_by(all, exp, gene, variant.id.x) %>% 
           summarize(n=(n()/prop.d)*100) %>% 
           group_by(all, exp, gene) %>% 
           summarize(n_mean=mean(n),
                     n_med=median(n),
                     n_uci=quantile(n, probs = 0.95),
                     n_lci=quantile(n, probs = 0.05),
                     snpset=unique(exp)))

# Data for simulated blue opsin gene genotypes based on RD
blop <- data.table(readRDS("/scratch/csm6hg/data/Crosses.daphnia11806-ra.sims.rds"))
blop.sum <- blop %>% group_by(snpset, exp) %>% summarize(num=n())  
blop <- data.table(blop %>% 
          group_by(sim.geno, snpset, exp) %>% 
          summarize(y=n()) %>% 
          left_join(blop.sum) %>% 
  mutate(y.prop=y/num*100,
         snpset=case_when(snpset=="TSP"~"SPs BLOP",
                          snpset=="Not TSP"~"Sim. BLOP"),
         variant.id=paste("test12", snpset),
         all=case_when(sim.geno=="het"~1,
                       sim.geno=="hom ref"~0,
                       sim.geno=="hom alt"~2)))

blop.tot <- blop[snpset=="SPs BLOP"] %>% 
  group_by(sim.geno, exp) %>% 
  summarize(y.prop=(sum(y)/sum(num))*100,
            snpset="Sim. BLOP",
            variant.id=paste("test12", snpset),
            all=case_when(sim.geno=="het" ~ 1,
                          sim.geno=="hom ref" ~ 0,
                          sim.geno=="hom alt" ~ 2))

# HWE expectation 
hwe <- data.table(all=c(0,1,2), y=c(25,50,25), snpset="HWE", variant.id="test")

# Distribution of genotype frequencies in F1s for blue opsin gene
f1_daphnia <- {
   toti %>%
    ggplot(., 
           aes(x=as.factor(all), 
               y=n_med, 
               color=snpset,
               group=gene), size=1.5) +
    geom_line(alpha=0.05) + 
    facet_wrap(~exp, ncol = 1) +
    geom_line(data=hwe, 
              aes(x=as.factor(all), y=y, color=snpset, group=1), size=1.5) +
    geom_line(data=blop.tot, 
              aes(x=as.factor(all), y=y.prop, color=snpset, group=1), size=1.5) +
    ylim(c(0,100)) +
    theme_bw() +
    labs(x="Allele dosage",
         y="F1 genotype frequency (%)",
         color="") +
    theme(legend.position = "bottom", 
          strip.text = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))
}

### Output
ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_BLOP_new.pdf", f1_daphnia, width = 8, height=6, dpi=300)

#### Table of frequencies over time ###

# Denominator for each exp
denom1 <- data.table(dat[all %in% c(0,1,2)][!is.na(variant.id.y)][gene=="Daphnia11806-RA"] %>% 
             group_by(exp) %>%       
             summarize(prop.d=n()))

denomFin <- data.table(dat[all %in% c(0,1,2)][!is.na(variant.id.y)][gene=="Daphnia11806-RA"] %>% 
             left_join(denom1) %>% 
             group_by(exp, geno) %>%       
             summarize(prop=n()/unique(prop.d)*100))

# Reorder
denomFin$geno <- factor(denomFin$geno, levels = c("hom_ref", "het", "hom_alt"))  
