# Trans species polymorphism - make CPD tree figure
# 3.12.2024
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(patchwork)

# Working directory
setwd("/scratch/csm6hg/data")

# Read in data
ds.tree2 <- data.table(readRDS("CPD_tsp_tree_figure.rds") %>% 
              mutate(position=as.numeric(tstrsplit(tstrsplit(id, "-")[[1]], ":")[[2]]) + 500,
                     chrom=tstrsplit(id,":")[[1]]) %>% 
              mutate(ch=paste(chrom,position,sep = "_")))

# Restrict to conservative SNP set
con <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
dat <- ds.tree2[ch %in% unique(con$ch)][!simpleAnnot=="Inter"][dataset=="500bps"]

# Plot median CPD distance output - 500 bps flanking focal SNP
dist.plot.med_a <- {
  
  dat %>% 
    mutate(classified1=case_when(classified=="shared_poly"~"TSP",
                                 TRUE~"Not TSP")) %>% 
    ggplot(., 
           aes(x=med_cpd,
               fill=classified1)) +
    geom_histogram(alpha=0.8, position="identity", color="black") +
    scale_fill_manual(values = c("TSP"="darkcyan",
                                 "Not TSP"='orange')) +
    theme_bw() +
    labs(x = "Median CPD (Within-Between)", 
         fill = "Focal SNP",
         y = "Number of Allele trees") +
    theme(legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=18), 
          strip.text = element_text(face="bold", size=18), 
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))
  
}

cpd.dist.tree_b <- {
  
  dat %>% 
    ggplot(., 
           aes(x=reorder(id, med_cpd), 
               y=med_cpd,
               color=classified1,
               group=1)) +
    geom_line() +
    geom_point(size=2) +
    facet_wrap(~classified1, nrow=1, scales = "free_x") +
    scale_color_manual(values = c("TSP"="darkcyan",
                                  "Not TSP"="gold")) +
    geom_hline(yintercept = 0, linetype=2) +
    theme_classic() +
    labs(x = "Allele tree", 
         color = "Focal SNP",
         y = "Median CPD (Within-Between)") +
    theme(legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=18), 
          strip.text = element_text(face="bold", size=18), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))
  
}

# Concatenate figures
p1 <- (dist.plot.med_a / cpd.dist.tree_b + plot_layout(guides = 'collect'))

# Save figure
ggsave("Figure_4_CPD.pdf", p1, width=8, height=8)

# Stats
t.test(data = dat, med_cpd~classified1)

# Normality tests
qqline(dat[classified1=="TSP"]$med_cpd)
qqline(dat[classified1=="Not TSP"]$med_cpd)

ks.test(dat[classified1=="TSP"]$med_cpd, y = "pnorm")
ks.test(dat[classified1=="Not TSP"]$med_cpd, y="pnorm")
