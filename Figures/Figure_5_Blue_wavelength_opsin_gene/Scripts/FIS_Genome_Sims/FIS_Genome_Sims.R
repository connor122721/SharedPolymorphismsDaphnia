# Trans species polymorphism - FIS
# 3.22.2024
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
q_blop=ecdf(gene.ag[exp=="Wild AxC F1s"][perm==0]$FIS)
q_blop(-0.5437215)

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

fis.plot2 <- {
  gene.ag[!exp=="Wild AxC F1s"] %>% 
    ggplot() +
    geom_density(aes(x = FIS, fill = sim), alpha=0.7) +
    geom_vline(data=gene.ag[sim=="Empirical"][!exp=="Wild AxC F1s"][gene %in% "Daphnia11806-RA"], 
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
ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_dist_crossesAll.pdf", fis.plot2, width = 6, height=12, dpi=300)

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
toti <- data.table(dat[all %in% c(0,1,2)][!is.na(variant.id.y)][!gene=="Daphnia11806-RA"] %>% 
           left_join(denom) %>% 
           group_by(all, exp, gene, variant.id.x) %>% 
           summarize(n=(n()/prop.d)*100) %>% 
           group_by(all, exp) %>% 
           summarize(n_mean=mean(n),
                     n_med=median(n),
                     n_uci=quantile(n, probs = 0.95),
                     n_lci=quantile(n, probs = 0.05),
                     snpset="Genome-wide"))

# Empirical data for BLOP
toti.blop <- data.table(dat[all %in% c(0,1,2)][!is.na(variant.id.y)][gene=="Daphnia11806-RA"] %>% 
                     left_join(denom) %>% 
                     group_by(all, exp, gene, variant.id.x, classified) %>% 
                     summarize(n=(n()/prop.d)*100) %>% 
                     group_by(all, exp, gene, variant.id.x, classified) %>% 
                     summarize(n_mean=mean(n),
                               n_med=median(n),
                               n_uci=quantile(n, probs = 0.95),
                               n_lci=quantile(n, probs = 0.05),
                               snpset="BLOP"))

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

blop.tot <- data.table(blop[snpset=="SPs BLOP"] %>% 
  group_by(sim.geno, exp) %>% 
  summarize(y.prop=(sum(y)/sum(num))*100,
            snpset="Sim. BLOP",
            variant.id=paste("test12", snpset),
            all=case_when(sim.geno=="het" ~ 1,
                          sim.geno=="hom ref" ~ 0,
                          sim.geno=="hom alt" ~ 2)))

# HWE expectation 
hwe <- data.table(all=c(0,1,2), y=c(25,50,25), snpset="HWE", variant.id="test")

# Distribution of genotype frequencies in F1s for blue opsin gene
f1_daphnia <- {
   toti[exp=="Wild AxC F1s"] %>%
    ggplot(., 
           aes(x=as.factor(all), 
               y=n_med, 
               color=snpset,
               group=1)) +
    geom_line(size=1.5) + 
    geom_line(data=hwe, 
              aes(x=as.factor(all), y=y, color=snpset, group=1), size=1.5) +
    geom_line(data=blop.tot[exp=="Wild AxC F1s"], 
              aes(x=as.factor(all), y=y.prop, color=snpset, group=1), size=1.5) +
    geom_line(data=toti.blop[exp=="Wild AxC F1s"][classified=="shared_poly"], 
              aes(x=as.factor(all), y=n_med, color=snpset, group=variant.id.x), size=1.5, alpha=0.7) +
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
ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_BLOP_new_crossAxCWild.pdf", f1_daphnia, width = 8, height=4, dpi=300)

# Distribution of mutations for BLOP across crosses
f1_daphnia_all <- {
  toti[!exp=="Wild AxC F1s"] %>%
    ggplot(., 
           aes(x=as.factor(all), 
               y=n_med, 
               color=snpset,
               group=1)) +
    facet_wrap(~exp, ncol = 1) +
    geom_line(size=1.5) + 
    geom_line(data=hwe, 
              aes(x=as.factor(all), y=y, color=snpset, group=1), size=1.5) +
    geom_line(data=blop.tot[!exp=="Wild AxC F1s"], 
              aes(x=as.factor(all), y=y.prop, color=snpset, group=1), size=1.5) +
    geom_line(data=toti.blop[!exp=="Wild AxC F1s"][classified=="shared_poly"], 
              aes(x=as.factor(all), y=n_med, color=snpset, group=variant.id.x), size=1.5, alpha=0.7) +
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

ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_BLOP_new_crossAll.pdf", f1_daphnia_all, width = 6, height=12, dpi=300)

#### Table of frequency in F1s ###

# Metadata
fin <- data.table(read.csv("/project/berglandlab/connor/metadata/samples.9.8.22.csv"))
samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
sc <- samps[Nonindependent==0]
sc[,SC.uniq:=SC]
sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
sc.ag <- sc[LabGenerated==F & Species=="pulex", list(clone=clone[which.max(medrd)][1]), 
            list(SC.uniq, Species, AxCF1Hybrid)]
samps[LabGenerated==1][SC=="selfedC"]  # CxC lab
samps[LabGenerated==1][AxCF1Hybrid==1] # AxC lab

# Metadata of blop genotype
blop.meta <- data.table(fread("/scratch/csm6hg/data/euroPulex.haplo.clustk3.bluopsin.txt") %>% 
            left_join(samps, by=c("Sample"="clone")))

dat.blop <- data.table(dat[exp=="Wild AxC F1s"] %>% 
                  distinct(Sample, exp) %>% 
                  left_join(blop.meta, by="Sample") %>% 
                  filter(Genotype %in% c("1|1", "1|2", "2|2") &
                         WildSequenced == "1" & year=="2018") %>% 
                  mutate(sc_n=case_when(SC=="OO"~"OO",
                                        TRUE ~ "SC")))

dat.blop %>% 
  group_by(sc_n, Genotype) %>% 
  summarize(n())

# How many were wild-samples collected in 2018?
samps[year=="2018"][population %in% c("D8","DBunk")][WildSequenced==1]

