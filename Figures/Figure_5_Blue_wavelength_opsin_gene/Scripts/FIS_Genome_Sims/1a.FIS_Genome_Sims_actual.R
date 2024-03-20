# Trans species polymorphism - FIS
# 1.30.2024
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)

# Read data
o.ag <- readRDS("/scratch/csm6hg/data/AxC_geneclass_conservative_100boot_Wild.rds")

# Gene averages
gene.ag <- data.table(readRDS("/scratch/csm6hg/data/AxC_geneclass_conservative_100boot_Wild_raw.rds"))
gene.ag <- gene.ag %>% 
  mutate(snpset=case_when(snpset %like% ";" ~ "Contains Both",
                          snpset =="Genome-wide" ~ "Not TSP",
                          TRUE~snpset))

# Blue opsin wavelength gene
gene.ag[perm==0][gene %in% "Daphnia11806-RA"]

# Add metadata
dt <- data.table(o.ag%>% 
              mutate(pp=case_when(perm==0~"Empirical",
                                  !perm==0~"Simulation"),
                     snpset=case_when(snpset %like% ";" ~ "Contains Both",
                                      TRUE~snpset)))

# Reorder
dt$snpset <- factor(dt$snpset, levels = c("Not TSP", "Contains Both", "TSP"))
gene.ag$snpset <- factor(gene.ag$snpset, levels = c("Not TSP", "Contains Both", "TSP"))  

# Plot
fis.plot <- { 
  
  gene.ag[FIS>-0.7 & FIS <0.7][N>=2][!snpset=="Contains Both"] %>% 
    ggplot() +
    geom_density(aes(x=FIS, fill=snpset), alpha=0.7) +
    geom_vline(data=gene.ag[perm==0][gene %in% "Daphnia11806-RA"], 
               aes(xintercept = FIS, size=1.1),
               color="steelblue2") +
    scale_fill_manual(values = c("Not TSP"="darkgreen",
                                 "TSP"="rosybrown")) +
    theme_bw() +
    labs(title = "", 
         y="Density",
         fill = "",
         x = expression(italic(italic("F")["IS"]))) +
    theme(legend.position = "none",
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=18), 
          strip.text = element_text(face="bold", size=18), 
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))
  }

# Enrichment
library(clusterProfiler)
library(enrichplot)
genelist <- paste("pulex.euro_",unique(gene.ag[FIS< -0.3][perm==0][N>=2][!snpset=="Not TSP"]$gene), sep="")
background <- paste("pulex.euro_",unique(gene.ag[FIS> -0.3][perm==0][!snpset=="Not TSP"]$gene), sep="")

spp.gmap = readRDS(list.files(path = "/scratch/csm6hg/genomes/proteins_species/hmmer_sep8", 
                              pattern = paste("gmap.", "pulex.euro",".rds", sep=""), 
                              full.names = T))

# Test enrichment
go.test <- enricher(gene = genelist,
                    pvalueCutoff = 0.01,
                    pAdjustMethod = "fdr",
                    minGSSize = 10,
                    maxGSSize = 500,
                    qvalueCutoff = 0.01,
                    universe = background, 
                    TERM2GENE = spp.gmap)

plot.go <- {pp %>% 
    filter(p.adjust<0.01) %>% 
    arrange(desc(p.adjust)) %>% 
    ggplot(aes(x = log2(odd),
               y = reorder(term, log2(odd)),
               size = Count,
               color = log10(p.adjust))) +
    geom_point() +
    labs(x="log2(Enrichment)", 
         y="",
         color="FDR P-value",
         size="Candidate gene number") +
    theme_bw() +
    scale_color_viridis(option="magma") +
    #paletteer::scale_color_paletteer_c("gameofthrones::targaryen") +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=14),
          legend.title = element_text(face="bold", size=16),
          legend.background = element_blank(),
          axis.text.x = element_text(face="bold", size=15),
          axis.text.y = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18))}

# Data for blue opsin gene
dat <- readRDS("/project/berglandlab/connor/chapter1/paralogs/AxC_tot_het_withmeta_Daphnia11806-RA_Daphnia00056-RA.rds")
tot <- readRDS("/project/berglandlab/connor/chapter1/paralogs/AxC_tot_het_withmeta.rds")

# Restrict to Wildcaught F1s
fin <- data.table(read.csv("/project/berglandlab/connor/metadata/samples.9.8.22.csv"))
dat <- data.table(dat %>% left_join(fin %>% select(c(id,WildSequenced,year,pondID)), by=c("Sample"="id")))
tot <- data.table(tot %>% left_join(fin %>% select(c(id,WildSequenced,year,pondID)), by=c("Sample"="id")))

# Restrict to Wild-caught F1s
dat <- data.table(dat[WildSequenced==1])
tot <- data.table(tot[WildSequenced==1])

# Data for simulated blue opsin gene genotypes based on RD
blop <- data.table(readRDS("/scratch/csm6hg/data/daphnia11806-ra.sims.rds"))
blop <- data.table(blop %>% left_join(fin %>% select(c(id,WildSequenced,year,pondID)), by=c("Sample"="id")))
blop <- data.table(blop[WildSequenced==1])

blop.sum <- blop %>% group_by(snpset) %>% summarize(num=n())  
blop <- data.table(blop %>% 
          group_by(sim.geno, snpset) %>% 
          summarize(y=n()) %>% 
          left_join(blop.sum) %>% 
  mutate(y.prop=y/num*100,
         snpset=case_when(snpset=="TSP"~"TSP - Sim. BLOP",
                          snpset=="Not TSP"~"Not TSP - Sim. BLOP"),
         variant.id=paste("test12", snpset),
         all=case_when(sim.geno=="het"~1,
                       sim.geno=="hom ref"~0,
                       sim.geno=="hom alt"~2)))

blop.tot <- blop %>% 
  group_by(sim.geno) %>% 
  summarize(y.prop=sum(y)/744*100,
            snpset="Sim. BLOP",
            variant.id=paste("test12",snpset),
            all=case_when(sim.geno=="het"~1,
                          sim.geno=="hom ref"~0,
                          sim.geno=="hom alt"~2))

# Restrict to conservative SNPs
con <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
dat <- dat[ch %in% unique(con$ch)]
dat.tot <- tot[ch %in% unique(con$ch)]

# Distribution of F1s for genome
toti <- data.table(dat.tot[all %in% c(0,1,2)] %>% 
                     group_by(all, snpset, gene, variant.id) %>% 
                     summarize(n=n()) %>% 
                     group_by(all, snpset) %>% 
                     summarize(mean=(mean(n)/48)*100) %>% 
                     mutate(snpset=paste("Genome", snpset),
                            variant.id=paste("test1", snpset)))

# HWE expectation 
hwe <- data.table(all=c(0,1,2), y=c(25,50,25), snpset="HWE", variant.id="test")

# Distribution of genotype frequencies in F1s for blue opsin gene
f1_daphnia <- {
  dat[all %in% c(0,1,2)][!snpset=="Not TSP"]%>% 
    group_by(all, variant.id, snpset, col, gene) %>% 
    summarize(mean=n(),
              snpset="TSP - Blue opsin gene") %>% 
    ggplot(., 
           aes(x=as.factor(all), 
               y=(mean/48)*100, 
               color=snpset, 
               group=variant.id)) +
    geom_line(data=hwe, aes(x=as.factor(all), y=y, color=snpset), size=1.5) +
    geom_line(data=toti, aes(x=as.factor(all), y=mean, color=snpset), size=1.5) +
    geom_line(data=blop.tot, aes(x=as.factor(all), y=y.prop, color=snpset), size=1.5) +
    geom_line() +
    geom_point() + 
    scale_color_manual(values = c("HWE"="black", 
                                  "TSP - Blue opsin gene"="steelblue2",
                                  "Sim. BLOP"="darkgreen",
                                  "Genome TSP"="rosybrown", 
                                  "Genome Not TSP"= "blue")) +
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
ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_Dist.pdf", plot = fis.plot, width=8, height=6, dpi=300)
ggsave("/scratch/csm6hg/figs/SuppFigure_FIS_BLOP_new.pdf", f1_daphnia, width = 8, height=6, dpi=300)
