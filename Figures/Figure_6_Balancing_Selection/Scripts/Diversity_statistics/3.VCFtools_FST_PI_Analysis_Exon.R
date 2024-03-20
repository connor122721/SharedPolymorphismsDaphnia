# Read in betascan output and visualize
# Connor Murray 10.26.2023
# ijob -c 10 --mem=10G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(tidyverse)
library(foreach)
library(colortools)
library(forcats)
library(readxl)
library(cowplot)
require(scales)

# Working directory
setwd('/project/berglandlab/connor/')

# Sample metadata
meta <- data.table(fread("metadata/samples.9.8.22.csv"))

# SNP metadata 
tot <- data.table(readRDS(file="/project/berglandlab/connor/data/classified_snps_filt.rds"))

# Restrict to conservative TSP set
con.tsp <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
tot <- tot[ch %in% con.tsp$ch]

# Read in gtf file
gene.gtf <- data.table(fread("../daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))

# Read in gene annotations
panth <- data.table(read_excel("../daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Write list of exons
gene.exon <- data.table(gene.gtf[sec=="exon"][chrom %in% unique(tot$chrom)][gene %in% unique(tot$gene)] %>% 
  select(chrom, start, stop, sec, gene) %>% 
  group_by(gene) %>% 
  mutate(exon.id=paste(gene, sec, row_number(), sep="_"),
         exon.size=stop-start,
         exon.id2=paste(chrom,start,stop, sep="_")) %>% 
  group_by(exon.id2) %>% 
  mutate(num.exon=n()))

# Remove duplicates
gene.exon <- gene.exon[!duplicated(exon.id2)]

#write.table(gene.exon, file = "data/exon.list", quote = F, sep = "\t", row.names = F, col.names = F)

# Overlap with exons
gene.gtfi <- data.table(gene.exon %>% select(start, stop, chrom, gene, exon.id))
toti <- tot %>% select(start=position, stop=position, chrom)
setkey(gene.gtfi, chrom, start, stop)

# Calculate exon length
laps <- na.omit(foverlaps(toti, gene.gtfi, type="within") %>% 
                  mutate(exon.size=abs(stop-start)) %>% 
                  select(-c(i.stop)))
colnames(laps)[c(2:3,6)] <- c("exon.start", "exon.stop", "position")

# Merge w/metadata
tottest <- data.table(tot %>% 
            left_join(laps[!duplicated(laps)], 
                      by=c("chrom", "position", "gene")))

# Number of tsps / exon
tot.exon.tsp <- data.table(tottest %>% 
            mutate(classified=case_when(classified == "shared_poly" ~ "TSP",
                                        classified %like% "poly" ~ "Not TSP",
                                        TRUE ~ "Other")) %>% 
            group_by(gene, exon.id, exon.size, exon.start, exon.stop, classified) %>% 
            summarize(num.tsp=length(variant.id)))

# Exons w/TSPs
exons.tsp.ex <- unique(tot.exon.tsp[classified=="TSP"]$exon.id) 

# Some exons have both TSP and Not TSP - remove that
tot.exon.tsp <- data.table(tot.exon.tsp[!c(classified=="Not TSP" & exon.id %in% exons.tsp.ex)][!classified=="Other"])

# Number of tsps / gene
tot.gene.tsp <- data.table(tot[classified=="shared_poly"] %>% 
        group_by(gene) %>% 
        summarize(gene.tsp=mean(length(unique(variant.id)))) %>% 
        left_join(gene.gtf %>% 
              filter(sec=="transcript") %>% 
              select(gene, gene.start=start, gene.stop=stop), 
              by=c("gene")) %>% 
        mutate(tsp.stan=(gene.tsp)/(gene.stop-gene.start)))

# Bind gene & exon
tot.tsp <- data.table(na.omit(tot.exon.tsp %>% 
          left_join(tot.gene.tsp, by=c("gene"))))

# output files fst
files.fst <- system("ls -f -R /project/berglandlab/connor/new_vcf2/vcftools/*_fst_nomiss_*", intern = TRUE)
list.fst <- lapply(files.fst, fread)
setattr(list.fst, 'names', system("ls -f -R /project/berglandlab/connor/new_vcf2/vcftools/*_fst_nomiss_*", intern = TRUE))

# Bind list 
fst.dt <- data.table(rbindlist(list.fst, use.names = T, idcol = T) %>% 
                       mutate(cont=case_when(.id %like% "euro" ~ "Europe D. pulex",
                                             .id %like% "nam" ~ "North America D. pulex",
                                             .id %like% "spp" ~ "Across species")) %>% 
                       select(-c(".id")))
colnames(fst.dt)[1:9] <- c("chrom", "start", "stop", "num.sites", "weighted_fst",
                           "fst", "exon.start", "exon.stop", "gene")

# Remove duplicated exons
fst.dt1 <- data.table(fst.dt %>% mutate(exon.id2=paste(chrom,exon.start,exon.stop,sep="_")))
fst.dt1 <- fst.dt1[!duplicated(exon.id2)]

# Summarize fst by gene
fst.dt <- na.omit(fst.dt1 %>% 
  left_join(tot.tsp, 
            by=c("gene", "exon.start", "exon.stop")))

fst.tot.cat <- data.table(na.omit(fst.dt %>% 
  mutate(exon.size.cat=case_when(exon.size<=100 ~ "1.Small (< 0.1 kb)",
                                 exon.size>100 & exon.size<=1000 ~ "2.Intermediate (0.1-1 kb)",
                                 exon.size>1000 ~ "3.Large (> 1 kb)")) %>% 
    group_by(gene, cont, exon.id, exon.size.cat, classified) %>% 
    summarize(fst=mean(fst))))

#saveRDS(fst.tot.cat, file = "data/fst_comp.rds")

############# PI #############
files <- system("ls -f -R /project/berglandlab/connor/new_vcf2/vcftools/*_pi_nomiss_*", intern = TRUE)
listi <- lapply(files, fread)
setattr(listi, 'names', list.files(path = "/project/berglandlab/connor/new_vcf2/vcftools", pattern = "*_pi_nomiss_*"))

# Bind list 
pi.dt <- data.table(rbindlist(listi, use.names = T, idcol = T) %>% 
                       mutate(cont=case_when(.id %like% "euro" ~ "Europe D. pulex",
                                             .id %like% "nam" ~ "North America D. pulex",
                                             .id %like% "spp" ~ "Across Species")) %>% 
                       select(-c(".id")))
colnames(pi.dt)[1:8] <- c("chrom", "start", "stop", "num.sites", "pi", 
                           "exon.start", "exon.stop", "gene")

# Summarize pi by gene
pi.dt <- pi.dt %>% 
  left_join(tot.tsp, 
            by=c("gene", "exon.start", "exon.stop"))

pi.tot <- data.table(pi.dt %>% 
  group_by(classified, exon.size, cont) %>% 
  summarize(pi=mean(pi, na.rm = T)))

pi.tot.cat <- na.omit(pi.dt %>% 
  mutate(exon.size.cat=case_when(exon.size<=100 ~ "1.Small (<0.1 kb)",
                                 exon.size>100 & exon.size<=1000 ~ "2.Intermediate (0.1-1 kb)",
                                 exon.size>1000 ~ "3.Large (>1 kb)")) %>% 
  group_by(gene, cont, exon.id, exon.size.cat, classified) %>% 
  summarize(pi=mean(pi)))

#saveRDS(pi.tot.cat, file = "data/pi_cat_exon.rds")

################ Tajima's D ###############
files.taj <- system("ls -f -R /project/berglandlab/connor/new_vcf2/vcftools/*_tajimaD_nomiss_*", intern = TRUE)
list.taj <- lapply(files.taj, fread)
setattr(list.taj, 'names', system("ls -f -R /project/berglandlab/connor/new_vcf2/vcftools/*_tajimaD_nomiss_*", intern = TRUE))

# Bind list 
taj.dt <- data.table(rbindlist(list.taj, use.names = T, idcol = T) %>% 
                       mutate(cont=case_when(.id %like% "euro" ~ "Europe D. pulex",
                                             .id %like% "nam" ~ "North America D. pulex")) %>% 
                       select(-c(".id")))
colnames(taj.dt)[1:7] <- c("chrom", "start", "num.sites", "taj",
                           "exon.start", "exon.stop", "gene")

# Summarize taj by gene
taj.dt <- taj.dt %>% 
  left_join(tot.tsp, 
            by=c("gene", "exon.start", "exon.stop"))

taj.tot <- data.table(taj.dt %>% 
                        group_by(classified, exon.size, cont) %>% 
                        summarize(taj=mean(taj, na.rm = T),
                                  taj.stan=mean(taj, na.rm = T)/num.sites))

taj.tot.cat <- na.omit(taj.dt %>% 
  mutate(exon.size.cat=case_when(exon.size <= 100 ~ "1.Small (<0.1 kb)",
                                 exon.size > 100 & exon.size <= 1000 ~ "2.Intermediate (0.1-1 kb)",
                                 exon.size > 1000 ~ "3.Large (>1 kb)")) %>% 
  group_by(gene, cont, exon.id, exon.size.cat, classified) %>% 
  summarize(taj=mean(taj/num.sites)))

#saveRDS(taj.tot.cat, file = "data/tajD_cat_exon.rds")

# Combine
exon.stats <- data.table(fst.tot.cat %>% 
  full_join(pi.tot.cat, by=c("exon.id", "gene", "cont", "exon.size.cat", "classified")) %>% 
  full_join(taj.tot.cat, by=c("exon.id", "gene", "cont", "exon.size.cat", "classified")))

#saveRDS(exon.stats, "/scratch/csm6hg/data/exon.diversity.stats.conservative.rds")
