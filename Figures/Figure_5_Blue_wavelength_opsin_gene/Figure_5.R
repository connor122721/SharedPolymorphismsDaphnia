# Trans species polymorphism - gene trees
# 1.30.24
# ijob -c 10 --mem=100G -p standard -A berglandlab_standard
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(ape)
library(stringr)
library(stringi)
library(magrittr)
library(tidyverse)
library(ggtree)
library(foreach)
library(colortools)
library(cowplot)
library(viridis)
library(readxl)
library(patchwork)
library(adegenet)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Read in gtf file
gene.gtf <- data.table(fread("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

# Get gene info
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))

# SNP metadata
tot <- readRDS(file="/project/berglandlab/connor/data/tot_snpsclassified_0.01thresh.RDS")

# Trees
files.tree <- system("ls -f -R /project/berglandlab/connor/candgene/tree/*treefile", intern = T)

  # Read tree
  tre.n <- read.tree(files.tree[i])
  
  # Only include samples in metadata
  tre.n <- keep.tip(phy = tre.n, tip = tre.n$tip.label[!tre.n$tip.label %like% "D.magna"])
  
  # Color wheel
  color <- wheel(color = "tomato", num = 15)
  
  # Open tree data
  tree <- data.table(File=as.character(tre.n$tip.label))
  
  # Extract naming
  tree <- data.table(tree %>%
                  mutate(Continent=tstrsplit(File,"_", fixed=T)[[1]],
                    Sample=as.character(tstrsplit(File, ".", fixed=T)[[1]]),
                    haplotype=stri_reverse(tstrsplit(stri_reverse(File), ".", fixed=T)[[4]]),
                    Region=paste(stri_reverse(tstrsplit(stri_reverse(File), ".", fixed=T)[[3]]),
                                 stri_reverse(tstrsplit(stri_reverse(File), ".", fixed=T)[[2]]),
                                 stri_reverse(tstrsplit(stri_reverse(File), ".", fixed=T)[[1]]), sep="."),
                    Start=stri_reverse(tstrsplit(stri_reverse(File), ".", fixed=T)[[2]]),
                    Stop=stri_reverse(tstrsplit(stri_reverse(File), ".", fixed=T)[[1]])) %>% 
                  mutate(Sample=case_when(Sample %like% "Europe_" ~ 
                                            str_remove_all(Sample, pattern = "Europe_"),
                                          Sample %like% "NorthAmerica_" ~ 
                                            str_remove_all(Sample, pattern = "NorthAmerica_"),
                                          TRUE ~ Sample)))
  
  # Any missing comparisons?
  missing <- tree$Sample[which(unique(tree$Sample) %in% unique(fin$Sample)==F)]
  
  # Merge with metadata
  tree <- data.table(rbind(na.omit(tree %>%
                               left_join(fin %>%
                                           select(Sample, Continent, cont, Species, Origin) %>%
                                           mutate(Sample=as.character(Sample)), by=c("Sample", "Continent"))) %>%
                       mutate(sample.n=paste(cont, Sample, sep=".")), 
            tree[File=="D.magna"], fill=T))
  
  # Color tips by Origin & Continent & Species
  tips <- data.table(tree %>%
                       mutate(color.con = case_when(Continent == "Europe" ~ color[1],
                                                    Continent == "NorthAmerica" ~ color[5],
                                                    File == "D.magna" ~ "Black"),
                              color.spp = case_when(Species == "Daphnia pulex" ~ color[1],
                                                    Species == "Daphnia pulicaria" ~ color[2],
                                                    Species == "Daphnia obtusa" ~ color[3],
                                                    File == "D.magna" ~ "Black"),
                              color.sppcont = case_when(cont == "Daphnia.pulex.NorthAmerica" ~ color[1],
                                                    cont == "Daphnia.pulexcaria.NorthAmerica" ~ color[2],
                                                    cont == "Daphnia.pulex.Europe" ~ color[5],
                                                    File == "D.magna" ~ "Black")))
  
  # Extract name
  name.fasta <- unique(tree$Region)[1]
  gene.name <- unique(gene.gtf[chrom == tstrsplit(unique(tree$Region)[1], ".", 
              fixed=T)[[1]]][start==unique(tree$Start)[1]][stop==unique(tree$Stop)[1]]$gene)
  
  # Read in aligned haplotypes
  msa <- fasta2genlight(paste("candgene/tree/", name.fasta, ".aln.fa", sep=""))
  msa.dt <- data.frame(tab(msa), sampleid = rownames(tab(msa)))
  #saveRDS(msa.dt, file = "data/msadt_Daphnia00056-RA.rds")
  
  # Wide to long format
  msa.dt2 <- data.table(msa.dt %>%  
        separate(sampleid, remove = F, into = c("sp"), sep = "_") %>% 
        melt(id = c("sampleid", "sp")) %>% 
        mutate(Continent=tstrsplit(sampleid,"_", fixed=T)[[1]],
            Sample=as.character(tstrsplit(sampleid, ".", fixed=T)[[1]]),
            haplotype=stri_reverse(tstrsplit(stri_reverse(sampleid), ".", fixed=T)[[4]]),
            Region=paste(stri_reverse(tstrsplit(stri_reverse(sampleid), ".", fixed=T)[[3]]),
                         stri_reverse(tstrsplit(stri_reverse(sampleid), ".", fixed=T)[[2]]),
                         stri_reverse(tstrsplit(stri_reverse(sampleid), ".", fixed=T)[[1]]), sep="."),
            Start=as.numeric(stri_reverse(tstrsplit(stri_reverse(sampleid), ".", fixed=T)[[2]])),
            Stop=as.numeric(stri_reverse(tstrsplit(stri_reverse(sampleid), ".", fixed=T)[[1]]))) %>% 
        mutate(Sample=case_when(Sample %like% "Daphnia.pulex.Europe_" ~ 
                                  str_remove_all(Sample, pattern = "Daphnia.pulex.Europe_"),
                                Sample %like% "Daphnia.pulex.NorthAmerica_" ~ 
                                  str_remove_all(Sample, pattern = "Daphnia.pulex.NorthAmerica_"),
                                TRUE ~ Sample),
               chrom=tstrsplit(Region, ".", fixed=T)[[1]]) %>% 
        mutate(variant=as.numeric(str_remove_all(tstrsplit(variable, ".", fixed=T)[[1]], pattern = "X")),
               ref=tstrsplit(variable, ".", fixed=T)[[2]],
               alt=tstrsplit(variable, ".", fixed=T)[[3]]) %>% 
        mutate(position=(variant+Start-1)) %>% 
        left_join(tot %>% select(-c(variant.id, col)), by=c("position", "chrom")))
  
  ### Gene tree and MSA ###
  
  # Candidate gene msa plot
  msa.plot <- {data.frame(msa.dt[-c(dim(msa.dt)[2])], sampleid = rownames(tab(msa))) %>% 
    separate(sampleid, remove = F, into = c("sp"), sep = "_") %>% 
    melt(id = c("sampleid", "sp") ) %>% 
    ggplot(.,
           aes(
      x=variable,
      y=sampleid,
      fill = as.factor(value))) +
    geom_tile() +
    labs(x="", 
         y="", 
         fill="") +
    ggtitle(name.fasta) +
    theme_classic() +
    facet_wrap(sp~., scales = "free", ncol = 1) +
    scale_fill_manual(values = c("skyblue", "brown")) +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold.italic", size=16),
          legend.background = element_blank(),
          text = element_text(face="bold.italic", size=16),
          axis.text.x=element_text(face="bold", size=14, angle=-40),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())}
  
  # Tree plot
  p <- {
    ggtree(tre.n,
              size=1.7, 
              root.position = 1) %<+% 
    tips +
    geom_tippoint(aes(color = cont), 
                  size=5) +
    geom_nodelab(aes(x=branch, label=label), 
                 vjust=-.5, size=4, fontface=2) +
    scale_color_manual(values = c("skyblue", "brown", "black")) +
    labs(color="", title=paste(gene.name, name.fasta)) +
    geom_treescale(x = 1.02, linesize = 1.5, fontsize = 6) +
    theme_tree(title = element_text(face="bold", size=20),
               legend.position='none',
               legend.text = element_text(face="bold.italic", size=16),
               legend.background = element_blank())
  }
  
  # Tree + MSA plot
  pp <- gheatmap(p, 
                 msa.dt[-c(dim(msa.dt)[2])],
                 colnames = F, 
                 width=1,
                 high="purple",
                 low="orange")
  
  ### Gene structure ###
  
  # Gene structure plot
  #saveRDS(gene.gtf[gene==gene.name][sec=="exon"], file="data/Daphnia00056-RA.exon.list.rds")
  gene.struc <- {gene.gtf[gene==gene.name][sec=="exon"] %>% 
    ggplot(aes(start/1000, y=paste("Exon", start, sep="_"))) +
    geom_linerange(aes(xmin=start/1000, xmax=stop/1000), size=4, color='blueviolet') +
    geom_vline(xintercept = 6350433/1000, linetype=2, size=1.5, color="red") +
      geom_vline(xintercept = 6351231/1000, linetype=2, size=1.5, color="red") +
    annotate('rect', xmin=6350433/1000, xmax=6351231/1000, 
             ymin = 0, ymax=9, alpha=0.2, fill='red') +
    theme_classic() +
    labs(x="Position (Kbp)",
         y="") +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.position = "none",
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}
  
  ### LD Matrix ###
  
  # Read in LD matrix
  ld.file <- system(paste("ls -f -R /scratch/csm6hg/matrix/pulex_euro_*", 
                          gene.name, "*.ld*", sep=""), intern = T)
  
  # Read in matrix
  ld <- data.table(fread(ld.file))
  ld <- ld %>% 
    left_join(tot[chrom==unique(ld$CHR_A)] %>% select(classified_A=classified, simpleAnnot_A=simpleAnnot, BP_A=position)) %>% 
    left_join(tot[chrom==unique(ld$CHR_A)] %>% select(classified_B=classified, simpleAnnot_B=simpleAnnot, BP_B=position))
  
  #saveRDS(ld, file = "data/LD_Daphnia00056-RA.euro.rds")
  
  # Plot LD Matrix
  ld.plot <- {ld[SNP_B %in% SNP_A] %>%
      ggplot(.) + 
      geom_tile(aes(x=SNP_A, y=SNP_B, fill=R2)) +
      scale_fill_viridis(option = "plasma") + 
      theme_classic() +
      coord_fixed(ratio = 1) +
      labs(x = "", 
           fill = expression(bold("r"^2)),
           y = "") +
      theme(legend.text = element_text(face="bold", size=16),
            legend.title = element_text(face="bold", size=20),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_text(face="bold", size=18),
            axis.title.y = element_text(face="bold", size=18),
            axis.title = element_text(face="bold", size=20))}
  
  ld.axis.x <- {ld[SNP_B %in% SNP_A] %>% 
    group_by(BP_A) %>% 
    distinct(classified_A) %>%
    mutate(classified_A=case_when(classified_A=="shared_poly" ~ "TSP",
                                  classified_A %like% "poly" ~ "Poly.",
                                  TRUE ~ "Other")) %>% 
    ggplot(.) +
    geom_tile(aes(x=as.factor(BP_A), y=1, fill=classified_A)) +
    scale_fill_manual(values = c("TSP"="darkcyan",
                                 "Poly."="beige")) + 
    theme_classic() +
    labs(x = "SNP A", 
         fill = "SNP classification",
         y = "") +
    theme(legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=18), 
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold", size=18, angle = -40),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}
  
  ld.axis.y <- {ld[SNP_B %in% SNP_A] %>% 
      group_by(BP_B) %>% 
      distinct(classified_B) %>%
      mutate(classified_B=case_when(classified_B=="shared_poly" ~ "TSP",
                                    classified_B %like% "poly" ~ "Poly.",
                                    TRUE ~ "Other")) %>% 
      ggplot(.) +
      geom_tile(aes(x=as.factor(BP_B), y=1, fill=classified_B)) +
      coord_flip() +
      scale_fill_manual(values = c("TSP"="darkcyan",
                                   "Poly."="beige")) + 
      theme_classic() +
      labs(x = "SNP B", 
           fill = "",
           y = "") +
      theme(strip.text = element_text(face="bold.italic", size=12),
            legend.text = element_text(face="bold", size=16),
            legend.position="none",
            legend.title = element_text(face="bold", size=20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=18),
            axis.title.y = element_text(face="bold", size=18),
            axis.title = element_text(face="bold", size=20))}
  
  # Master LD matrix
  ld.plot.tot <- { (ld.axis.y + 
                    ld.plot + 
                    plot_spacer() +
                    ld.axis.x + 
                    plot_layout(guides="collect", 
                                heights = c(1, 0.1),
                                widths = c(0.1, 1))) }
  
  ### Output figures ###
  comp.fig <- { ( gene.struc / ld.plot.tot | pp + 
                    plot_layout(guides="collect") ) }
  
  ggsave(paste("candgene/figs/Daphnia11806-RA.new.pdf", sep=""), comp.fig, width=25, height=16)
