# Read in betascan output and visualize
# Connor Murray 10.26.2023
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Packages
library(data.table)
library(tidyverse)
library(foreach)
library(colortools)
library(forcats)
library(ggridges)
library(cowplot)

# Working directory
setwd('/scratch/csm6hg/betascan/')

# betascan output files
file1 <- list.files(pattern = "euro.betafin",
                    recursive = TRUE)

# betascan output files
file2 <- list.files(pattern = "nam.betafin",
                    recursive = TRUE)

len <- unique(length(file1), length(file2))

# Go through each file
fin.dt <- foreach(i=1:len, .combine = "rbind", .errorhandling = "pass") %do% {

    # Progress message
    print(i)
    
    # Load in betascan output europe
    tbl1 <- data.table(fread(file1[i]),
                      chrom=tstrsplit(file1[i], ".", fixed=T)[[1]],
                      snpset=tstrsplit(file1[i], ".", fixed=T)[[3]],
                      cont=tstrsplit(file1[i], ".", fixed=T)[[4]],
                      file=file1[i])
    colnames(tbl1)[2] <- "Beta"
    
    # Load in betascan output north america
    tbl2 <- data.table(fread(file2[i]),
                       chrom=tstrsplit(file2[i], ".", fixed=T)[[1]],
                       snpset=tstrsplit(file2[i], ".", fixed=T)[[3]],
                       cont=tstrsplit(file2[i], ".", fixed=T)[[4]],
                       file=file2[i])
    colnames(tbl2)[2] <- "Beta"
    
    
    fin <- rbind(tbl1, tbl2)
    
    # Finish
    return(fin)
    
}

#saveRDS(fin.dt, "/scratch/csm6hg/data/betascan.daphnia.phased.rds")
fin.dt <- readRDS("/scratch/csm6hg/data/betascan.daphnia.phased.rds")

# SNP metadata - TSPs
tot <- readRDS(file="/scratch/csm6hg/data/classified_snps_filt.rds")
fin.dt2 <- data.table(na.omit(fin.dt %>% left_join(tot, by=c("chrom", "Position"="position"))) %>% 
                mutate(snp=case_when(classified=="shared_poly"~"TSP",
                                     TRUE ~ "Not TSP"),
                       cont=case_when(cont=="nam"~ "NAm. D. pulex",
                                      cont=="euro"~ "Euro. D. pulex")))

# Restrict to conservative TSP set
con.tsp <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
fin.dt2 <- fin.dt2[ch %in% con.tsp$ch]
table(fin.dt2$cont)

# Gene average
fin.dt3 <- data.table(na.omit(fin.dt2 %>% 
              group_by(gene, snp, cont, simpleAnnot) %>% 
              summarize(Beta1=mean(Beta, na.rm = T),
                        sd=sd(Beta, na.rm = T),
                        n=n())) %>% 
              mutate(se=sd/sqrt(n)) %>% 
              group_by(snp, cont, simpleAnnot) %>% 
              summarize(Beta2=mean(Beta1, na.rm = T),
                        sd=sd(Beta1, na.rm = T),
                        n=n()) %>% 
              mutate(se=sd/sqrt(n)))

# Genome average
g.avg <- data.table(fin.dt2 %>% 
            group_by(cont) %>%
            summarize(mean=mean(Beta, na.rm = T)))

# Candidate gene avg.
cand.avg <- data.table(fin.dt2[gene=="Daphnia11806-RA"] %>% 
                         group_by(cont) %>% 
                         summarize(mean=mean(Beta, na.rm = T)))

# difference in Betastat across pop
beta.plot <- {
  fin.dt3[!simpleAnnot%in% c("Start", "Stop", "Splice")] %>% 
  ggplot(., 
       aes(x = simpleAnnot, 
           y = Beta2,
           color = as.character(snp),
           group = cont)) +
  geom_hline(data=g.avg[cont %like% "NAm."],
             aes(yintercept = mean), linetype=2, size=1) +
  geom_hline(data=g.avg[cont %like% "Euro."],
             aes(yintercept = mean), linetype=2, size=1) +
  geom_pointrange(aes(x=simpleAnnot, 
                      y=Beta2, 
                      ymin=Beta2-2*se, 
                      ymax=Beta2+2*se), 
               position = position_dodge2(0.6)) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("TSP"="turquoise",
                                "Not TSP"="deeppink4")) +
  facet_wrap(~cont,  nrow = 2) +
  labs(x = "SNP Classification",
       color="",
       y = "Beta1") +
  theme(strip.text = element_text(face="bold.italic", size=18),
        legend.text = element_text(face="bold", size=18),
        legend.position="bottom",
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20)) 
}

# Stats
t.test(fin.dt2[simpleAnnot%in%c("NS")][cont=="NAm. D. pulex"][snp=="TSP"]$Beta)
t.test(fin.dt2[simpleAnnot%in%c("NS")][cont=="Euro. D. pulex"][snp=="TSP"]$Beta)

t.test(fin.dt2[simpleAnnot%in%c("Syn")][cont=="NAm. D. pulex"][snp=="TSP"]$Beta)
t.test(fin.dt2[simpleAnnot%in%c("Syn")][cont=="Euro. D. pulex"][snp=="TSP"]$Beta)

# Output
ggsave(plot = beta.plot, 
       filename = "/scratch/csm6hg/figs/betascan.species.genome.gene.se.conservative.pdf", 
       width = 7, 
       height = 7)
