# Reference allele bias at heterozygous sites - Supplemental Figure 1
# Connor Murray 10.17.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(patchwork)
library(ggalluvial)

# Working directory
setwd("/scratch/csm6hg/mapping/vcf/concordance/")

# Read in BUSCO gene SNP data
p.out <- data.table(read.csv("/scratch/csm6hg/data/refallelebias_subsamp_busco1000snps.csv"))

# Plot Reference allele frequency 
plot <- {
  p.out[af==0.5] %>% 
    ggplot(aes(x = paste(Species, Continent),
               y = dos.alt/dos.ref,
               fill = paste(Species, Continent)
    )) +
    geom_violin(draw_quantiles = c(0.025,0.5,0.975)) +
    coord_flip() +
    geom_hline(yintercept = 0.5, linetype=2, size=1.2) +
    labs(title="BUSCO gene SNPs",
         x="", 
         y="Allele Frequency",
         fill="Species") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(title = element_text(face="bold", size=15),
          legend.position = "none",
          legend.text = element_text(face="bold.italic", size=10),
          legend.title = element_text(face="bold", size=20),
          legend.background = element_blank(), 
          strip.text = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=16),
          axis.text.y = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18))}

# Data for alluvial plot
tot.new <- data.table(readRDS("/scratch/csm6hg/data/classified_snps_filt_remapped_ann.rds"))
tot.old <- data.table(readRDS("/scratch/csm6hg/data/classified_snps_filt.rds"))

tot.new1 <- tot.new[!classified %in% c("invariant", "filtered")]
tot.old1 <- tot.old[!classified %in% c("invariant", "filtered")]

# New mapped could be upwardly biased by not filtering out the same regions of the original mapping
bed.tot <- data.table(fread("/scratch/csm6hg/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed"))
colnames(bed.tot) <- c("chrom", "start", "stop")
bed.tot[,bed:="filt"]

### we need to do a little back-end work to align our snp object with the bed
tot.new1[,start:=position]
tot.new1[,stop:=position]

### we need to specify which columns we are merging these two files on
setkey(bed.tot, chrom, start, stop)
setkey(tot.new1, chrom, start, stop)

# Overlap with gene regions
sync.m <- foverlaps(tot.new1, bed.tot, nomatch=NA)

# Make non-overlapping regions = False
sync.m[is.na(bed), bed:=F]
tot.new2 <- sync.m[!bed=="filt"]

# Restrict to single copy orthologs identified between D. pulex species
sco <- fread("/scratch/csm6hg/genomes/proteins_species/primary_transcripts/OrthoFinder/Results_Aug31/euro.sco.genes.txt", header=F)
tot.new2 <- tot.new2[gene %in% sco$V1]
tot.old2 <- tot.old1[gene %in% sco$V1]

# Convert data
p1 <- data.table(rbind(table(tot.new2$classified),
                       table(tot.old2$classified)), 
                 ref=c("KAP4", "D84A")) %>% 
  pivot_longer(cols=c("fixed","polymorphic_1",
                      "polymorphic_2","shared_poly"))

# Change names
p2 <- data.table(p1 %>% 
  mutate(name=case_when(name=="filtered" ~ "Filtered",
                        name=="fixed" ~ "Fixed",
                        name=="invariant" ~ "Invariant",
                        name=="polymorphic_1" ~ "Private NAm.",
                        name=="polymorphic_2" ~ "Private Euro.",
                        name=="shared_poly" ~ "TSP")))

# Alluvial plot
tp1 <- {
  p2[!name=="Private Euro."] %>% 
    ggplot(.,
         aes(x = ref, 
             stratum = name, 
             alluvium = name,
             y = value,
             fill=name,
             label = name)) +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", size = 6) +
    theme_bw() +
    labs(x="Reference Assembly", 
         y="Number of SNPs",
         title="SNP Classifications for 2DSFS",
         fill="") +
    theme(title = element_text(face="bold", size=20),
          legend.text = element_text(face="bold", size=26),
          legend.title = element_text(face="bold", size=24),
          legend.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(face="bold", size=20),
          strip.text = element_text(face = "bold", size=16),
          axis.text.y = element_text(face="bold", size=20),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))
}

table(unique(tot.new2$ch) %in% unique(tot.old2$ch))
tot.new2[!which(unique(tot.new2$ch) %in% unique(tot.old2$ch))]
table(tot.new2[!which(unique(tot.new2$ch) %in% unique(tot.old2$ch))]$classified)

# Identify how SNP classification change according to reference assembly
m.new <- tot.new2 %>% select(chrom,position,classified,ch,gene,simpleAnnot)
m.old <- tot.old2 %>% select(chrom,position,classified,ch,gene,simpleAnnot)

m.tot <- data.table(m.new %>% 
  full_join(m.old, 
            by=c("chrom", "position", 
                 "gene", "ch", "simpleAnnot")) %>% 
    mutate(Var = case_when(
           is.na(classified.x) & is.na(classified.y) ~ "Missing in Both",
           is.na(classified.x) & classified.y %like% "polymorphic_2" ~ "Unique Private in D84A",
           is.na(classified.x) & classified.y %like% "shared" ~ "Unique TSPs in D84A",
           is.na(classified.x) & classified.y %like% "fixed" ~ "Unique Fixed in D84A",
           is.na(classified.y) & classified.x %like% "polymorphic_1" ~ "Unique Private in KAP4",
           is.na(classified.y) & classified.x %like% "shared" ~ "Unique TSPs in KAP4",
           is.na(classified.y) & classified.x %like% "fixed" ~ "Unique Fixed in KAP4",
           classified.x == classified.y ~ "Unchanged",
           classified.x == "polymorphic_1" & classified.y == "polymorphic_2" ~ "Fliped Poly.",
           classified.x == "polymorphic_2" & classified.y == "polymorphic_1" ~ "Fliped Poly.",
           classified.x %in% c("polymorphic_1", "polymorphic_2") & classified.y=="shared_poly" ~ "TSP to Private",
           classified.y %in% c("polymorphic_1", "polymorphic_2") & classified.x=="shared_poly" ~ "Private to TSP",
           classified.x %in% c("polymorphic_1", "polymorphic_2", "shared_poly") & classified.y=="fixed" ~ "Poly. to Fixed",
           classified.y %in% c("polymorphic_1", "polymorphic_2", "shared_poly") & classified.x=="fixed" ~ "Poly. to Fixed",
           TRUE ~ "UNK")))

# Pie chart of Whole dataset
pie <- {
  m.tot[!Var %like% "Unique"][!Var %like% "UNK"] %>%
    group_by(Var) %>% 
    summarize(n=length(Var)) %>% 
    ggplot(., 
           aes(x=1,
               y=n,
               fill=Var
           )) + 
    geom_bar(stat="identity", width=1, color="white") +
    geom_text(aes(label=paste(round((n/dim(m.tot[!Var %like% "Unique"][!Var %like% "UNK"])[1])*100, 1),
                              "%", sep="")),
              position = position_stack(vjust = 0.5), size=8) +
    coord_polar("y", start=0) + 
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    labs(x="",
         fill="Loci Category",
         title="Merged SNP dataset and classification",
         y="") +
    theme(text = element_text(face="bold", size=18),
          strip.text = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

m.eu <- data.table(m.tot[Var %like% "D84A"] %>%
  group_by(name=classified.y, Var) %>% 
  summarize(n=length(Var)) %>% 
  mutate(name=case_when(name=="filtered" ~ "Filtered",
                        name=="fixed" ~ "Fixed",
                        name=="invariant" ~ "Invariant",
                        name=="polymorphic_1" ~ "Private NAm.",
                        name=="polymorphic_2" ~ "Private Euro.",
                        name=="shared_poly" ~ "TSP"),
         brief=str_remove(Var, " in D84A")),
  ref="D84A")

m.na <- data.table(m.tot[Var %like% "KAP4"] %>% 
  group_by(name=classified.x, Var) %>% 
  summarize(n=length(Var))%>% 
  mutate(name=case_when(name=="filtered" ~ "Filtered",
                        name=="fixed" ~ "Fixed",
                        name=="invariant" ~ "Invariant",
                        name=="polymorphic_1" ~ "Private NAm.",
                        name=="polymorphic_2" ~ "Private Euro.",
                        name=="shared_poly" ~ "TSP"),
         brief=str_remove(Var, " in KAP4")),
  ref="KAP4")

# Bar chart of missing SNPs
bar <- {
    ggplot() + 
    geom_bar(data=m.eu, 
             aes(x=ref, y=n, fill=brief), 
             stat="identity", width=1, color="white") +
    geom_bar(data=m.na, 
             aes(x=ref, y=n, fill=brief), 
             stat="identity", width=1, color="white") +
    theme_bw() +
    scale_fill_manual(values = c("Unique Fixed"="purple",
                                  "Unique Private"="orange",
                                  "Unique TSPs"="steelblue")) +
    #facet_wrap(~ref) +
    labs(x="Reference Assembly",
         title="Unique SNPs, post-mapping",
         y="Number of SNPs",
         fill="") +
    theme(text = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18), 
          legend.position = "bottom",
          axis.text.x = element_text(face="bold", size=20),
          axis.text.y = element_text(face="bold", size=20),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))
}

### Compile large figure ###
ggsave(plot, "../figures/refallelebias_boxplot_subsamp_busco1000snps.pdf", width=12, height=8)
ggsave(plot=bar, "/scratch/csm6hg/figs/supp1_barplot_missing_SNPs.pdf", width=8, height=8)
ggsave(plot=tp1, "/scratch/csm6hg/figs/supp1_alluvial_SNPs.pdf", width=8, height=8)
ggsave(plot=pie, "/scratch/csm6hg/figs/supp1_pie_SNPs.pdf", width=8, height=8)
