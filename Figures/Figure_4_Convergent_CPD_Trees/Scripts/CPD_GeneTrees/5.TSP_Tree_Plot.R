# Trans species polymorphism - allele trees calculate distance
# 3.23.2023
# ijob -c 15 --mem=100G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

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
library(readxl)
library(adegenet)
library(patchwork)

# Register cores
doParallel::registerDoParallel(cores = 15)

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
tot <- readRDS(file="/project/berglandlab/connor/candgene/classified_snps_filt_exon.rds")

# Trees
files.tree <- system("ls -f -R /project/berglandlab/connor/candgene/tree_exon_tsp/*treefile", intern = T)

# Exon metadata
exon <- data.table(rbind(data.table(read.csv("candgene/exons.genome.list.tspset2", header = F))))
colnames(exon) <- c("ch", "gene", "classified")

# Remove Not TSPs
tsp.genes <- unique(exon[classified=="shared_poly"]$gene)
#exon <- exon[c(!classified=="shared_poly" & gene %in% tsp.genes)]

# Color wheel
color <- wheel(color = "tomato", num = 15)
 
# Go through each gene - make a composite figure
ds <- foreach(i=1:length(files.tree), .combine = "rbind", .errorhandling = "remove") %dopar% {
  
  # Start gene; i=133
  print(paste("Gene tree:", i, sep=" "))
    
  # Read tree
  tre.n <- read.tree(files.tree[i])
   
  # Open tree data
  tree <- data.table(File=as.character(tre.n$tip.label))
  
  # Extract naming
  tree1 <- data.table(tree %>%
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
  missing <- tree1$Sample[which(unique(tree1$Sample) %in% unique(fin$Sample)==F)]
  
  # Merge with metadata
  tree2 <- data.table(rbind(na.omit(tree1 %>%
                left_join(fin %>%
                select(Sample, Continent, cont, Species, Origin) %>%
                mutate(Sample=as.character(Sample)), by=c("Sample", "Continent"))) %>%
                mutate(sample.n=paste(cont, Sample, sep="."))))
    
  # Color tips by Origin & Continent & Species
  tips <- data.table(tree2 %>%
                         mutate(color.con = case_when(Continent == "Europe" ~ color[1],
                                                      Continent == "NorthAmerica" ~ color[5]),
                                color.spp = case_when(Species == "Daphnia pulex" ~ color[1]),
                                color.sppcont = case_when(cont == "Daphnia.pulex.NorthAmerica" ~ color[1],
                                                      cont == "Daphnia.pulex.Europe" ~ color[5])))
  
  # Drop missing tips
  tre.n1 <- keep.tip(phy = tre.n, tip = unique(tips[!Sample %in% missing]$File))
    
  # Extract name
  name.fasta <- unique(tree2$Region)[1]
  gene.name <- unique(gene.gtf[chrom == tstrsplit(unique(tree2$Region)[1], ".", 
                fixed=T)[[1]]][start == unique(tree2$Start)[1]][stop==unique(tree2$Stop)[1]]$gene)
  
  # Distance formula
  co.dist <- cophenetic.phylo(tre.n1)
  co.dist1 <- data.table(reshape2::melt(co.dist))
  colnames(co.dist1) <- c("samp1", "samp2", "x")
  
  # Add metadata & rm within sample comparison 
  co.dist2 <- data.table(co.dist1[!samp1==samp2] %>% 
                  left_join(tips, by=c("samp1"="File")) %>% 
                  left_join(tips, by=c("samp2"="File")) %>% 
              mutate(chrom.x=tstrsplit(Region.x, ".", fixed=T)[[1]],
                     chrom.y=tstrsplit(Region.y, ".", fixed=T)[[1]]) %>% 
              mutate(id_x1=paste(chrom.x, paste(Start.x, Stop.x, sep="-"), sep=":"),
                     id_y1=paste(chrom.y, paste(Start.y, Stop.y, sep="-"), sep=":")) %>% 
              left_join(exon, by=c("id_x1"="ch")) %>% 
              left_join(exon, by=c("id_y1"="ch")))
  
  # Distance per species & haplotype & class
  co.dist3 <- data.table(co.dist2 %>% 
    mutate(id=unique(id_x1), 
           classified=unique(classified.x)) %>% 
    mutate(cont_comp = case_when(cont.x=="Daphnia.pulex.Europe" & cont.y=="Daphnia.pulex.NorthAmerica" ~ "Between",
                                 cont.x=="Daphnia.pulex.NorthAmerica" & cont.y=="Daphnia.pulex.Europe" ~ "Between",
                                 cont.x %like% "Euro" & cont.y %like% "Euro" ~ "Within Euro",
                                 cont.x %like% "North" & cont.y %like% "North" ~ "Within NorthAmerica",
                                 TRUE ~ "Other")) %>%
    mutate(haplo_comp = case_when(cont_comp == "Within Euro" & !haplotype.x==haplotype.y & Sample.x==Sample.y 
                                  ~ "Within_euro",
                                  cont_comp == "Within Euro" & !Sample.x==Sample.y ~ "Between_euro",
                                  cont_comp == "Between" ~ "Between_species",
                                  cont_comp == "Within NorthAmerica" & !haplotype.x==haplotype.y & Sample.x==Sample.y
                                  ~ "Within_nam",
                                  cont_comp == "Within NorthAmerica" & !Sample.x==Sample.y ~ "Between_nam")) %>% 
    group_by(id, classified, cont_comp) %>% 
    summarize(mean.x=mean(x),
              median.x=median(x)), 
    i=i)

  # Finish
  return(co.dist3)
}
#ds <- data.table(readRDS("/project/berglandlab/connor/candgene/phylo.dist.828tsp.nohaplo.med.rds"))

# Calculate CPD(w-b)
ds1 <- data.table(na.omit(ds) %>% 
                    mutate(cont_comp=str_replace(cont_comp, " ", "_")) %>%
                    group_by(classified) %>% 
                    pivot_wider(values_from = c(mean.x, median.x), 
                                names_from = c(cont_comp)) %>% 
                    mutate(cpd_wb_euro_med = median.x_Within_Euro - median.x_Between,
                           cpd_wb_nam_med = median.x_Within_NorthAmerica - median.x_Between,
                           cpd_wb_euro_avg = mean.x_Within_Euro - mean.x_Between,
                           cpd_wb_nam_avg = mean.x_Within_NorthAmerica - mean.x_Between) %>% 
                    mutate(chrom=tstrsplit(id, ":")[[1]],
                           position.str=tstrsplit(id, ":")[[2]]) %>% 
                    mutate(position=as.numeric(tstrsplit(position.str,"-")[[1]])+500) %>% 
                    left_join(tot, by=c("chrom", "position", "classified")))

# Average
ds.tree2 <- data.table(ds1 %>% 
                         group_by(id, classified, simpleAnnot, gene) %>%
                         summarize(med_cpd=median(cpd_wb_euro_med, cpd_wb_nam_med),
                                   avg_cpd=mean(cpd_wb_euro_avg, cpd_wb_nam_avg)) %>% 
                         mutate(classified1=case_when(classified=="shared_poly"~"TSP",
                                                      TRUE~"Not TSP")))

ds.tree2 %>% group_by(classified1, simpleAnnot) %>% summarize(mean(avg_cpd), median(med_cpd))

# Save output
saveRDS(ds.tree2, file = "/project/berglandlab/connor/candgene/phylo.dist.500bptsp.all.nohaplo.med.rds")
#ds.tree2 <- data.table(readRDS("/project/berglandlab/connor/candgene/phylo.dist.552tsp.nohaplo.rds"))

# Go through each gene - make a composite figure
foreach(i=1:length(files.tree), .errorhandling = "remove") %do% {
  #i=133
  
  # Start gene
  print(paste("Gene tree:", i, sep=" "))
  
  # Read tree
  tre.n <- read.tree(files.tree[i])
  
  # Open tree data
  tree <- data.table(File=as.character(tre.n$tip.label))
  
  # Extract naming
  tree1 <- data.table(tree %>%
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
  missing <- tree1$Sample[which(unique(tree1$Sample) %in% unique(fin$Sample)==F)]
  
  # Merge with metadata
  tree2 <- data.table(rbind(na.omit(tree1 %>%
                                      left_join(fin %>%
                                                  select(Sample, Continent, cont, Species, Origin) %>%
                                                  mutate(Sample=as.character(Sample)), by=c("Sample", "Continent"))) %>%
                              mutate(sample.n=paste(cont, Sample, sep="."))))
  
  # Color tips by Origin & Continent & Species
  tips <- data.table(tree2 %>%
                       mutate(color.con = case_when(Continent == "Europe" ~ color[1],
                                                    Continent == "NorthAmerica" ~ color[5]),
                              color.spp = case_when(Species == "Daphnia pulex" ~ color[1]),
                              color.sppcont = case_when(cont == "Daphnia.pulex.NorthAmerica" ~ color[1],
                                                        cont == "Daphnia.pulex.Europe" ~ color[5])))
  
  # Drop missing tips
  tre.n1 <- keep.tip(phy = tre.n, tip = unique(tips[!Sample %in% missing]$File))
  
  # Extract name
  name.fasta <- unique(ds.tree2[bed==unique(tree2$Region)]$id)
  gene.name <- unique(gene.gtf[chrom == tstrsplit(unique(tree2$Region)[1], ".", 
                  fixed=T)[[1]]][start %in% unique(tree2$Start)[1]:unique(tree2$Stop)[1]][
                  stop %in% unique(tree2$Start)[1]:unique(tree2$Stop)[1]]$gene)
  Mean_cpd <- round(ds.tree2[bed==unique(tree2$Region)]$avg_cpd, digits = 3)
  Median_cpd <- round(ds.tree2[bed==unique(tree2$Region)]$med_cpd, digits = 3)
  snp <- ds.tree2[bed==unique(tree2$Region)]$classified
  bed.name <- unique(ds.tree2[bed==unique(tree2$Region)]$bed)
    
  ### Gene tree ###
  
  # Tree plot
  p <- {
    ggtree(tre.n1,
               size=1.7, 
               root.position = 1) %<+% 
      tips +
      geom_tippoint(aes(color = cont), 
                    size=4) +
      scale_x_ggtree() +
      geom_treescale(x=0, y=45, fontsize=4, linesize=1) +
      scale_color_manual(values = c("Daphnia.pulex.Europe"="skyblue", 
                                    "Daphnia.pulex.NorthAmerica"="brown", 
                                    "D.pulicaria"="black")) +
      #labs(color=paste("Med. CPD:", Median_cpd, "Avg. CPD:", Mean_cpd), 
      #     title=paste("Gene:", gene.name, name.fasta)) +
      theme_tree(title = element_text(face="bold", size=20),
                 legend.text = element_text(face="bold.italic", size=16),
                 legend.background = element_blank())
  }
  
  # Output
  ggsave(paste("candgene/TSP_Run_3_23_250/", snp, ".", bed.name, ".", i, ".pdf", sep=""), p, height = 10, width = 10)
}

# Plot avg distance output
dist.plot.avg <- {
  
  ds.tree2[!simpleAnnot =="Inter"] %>% 
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

cpd.dist.tree <- {
  
  ds.tree2 %>% 
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

p1 <- (dist.plot.avg / cpd.dist.tree + plot_layout(guides = 'collect'))

ggsave("candgene/cpd.nohaplo.500bps.nohaplo.median.pdf", p1, width=12, height=12)

ds.tree2[avg_cpd>0]

# Read in gene annotations
panth <- data.table(read_excel("../daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

panth[qseqid == "Daphnia11806-RA"]

