# Trans species polymorphism - LD decay
# 10.26.2023
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(viridis)
library(doParallel)
library(readxl)
library(cowplot)
library(adegenet)
library(seqinr)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
out <- c("new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Read in total SNPs
tot_snps <- data.table(fread("metadata/snps_new"))

# Register cores
#doParallel::registerDoParallel(cores = 10)

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samps <- seqGetData(genofile, var.name = "sample.id")
fin <- fin[Sample %in% samps]

# Read in gene annotations
panth <- data.table(read_excel("../daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Gene analyses
pro <- read.fasta("/project/berglandlab/Karen/genomefiles/Daphnia.proteins.aed.0.6.fasta", seqtype="AA")
pro <- data.table(gene=getName(pro), AA.length=getLength(pro))
pro <- pro[,splice:=tstrsplit(gene, "-")[[2]]]

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

# SNP metadata - TSPs
tot <- data.table(readRDS(file="/project/berglandlab/connor/data/classified_snps_filt.rds"))

# Restrict to conservative TSP set
con.tsp <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
tot <- tot[ch %in% con.tsp$ch]

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(fin[cont %in% as.character(compi[1])]$Sample)
samp2 <- as.vector(fin[cont %in% as.character(compi[2])]$Sample)

# Read in SNP LD info
files <- system("ls -f -R /project/berglandlab/connor/chapter1/linkage/out/*.ld", intern = TRUE)
listi <- lapply(files, fread)
setattr(listi, 'names', list.files(path = "/project/berglandlab/connor/chapter1/linkage/out", pattern = "ld"))

# Bind list 
dt <- rbindlist(listi, use.names = T, idcol = T) %>% 
  left_join(tot %>% select(-c(ch, col)), 
            by=c("BP_A"="position", 
                 "CHR_A"="chrom")) %>% 
  left_join(tot %>% select(-c(ch, col)), 
            by=c("BP_B"="position",
                 "CHR_B"="chrom")) %>% 
  mutate(dist = as.numeric(abs(BP_A - BP_B))) %>% 
  mutate(cont = tstrsplit(.id, "_")[[2]])

# Save output
# saveRDS(dt, "/scratch/csm6hg/data/linkage.r2genes.conservative.rds")
dt <- data.table(readRDS("/scratch/csm6hg/data/linkage.r2genes.conservative.rds"))

# Add intervals
upp=10
dt2 <- data.table(dt[gene.x==gene.y] %>% 
  #group_by(.id, SNP_A) %>% 
  mutate(distc = as.factor(cut(dist, 
                breaks = seq(from = 0,
                             to = max(dist) + upp,
                             by = upp), right = F))))

# Add Gene x Gene comps
dt3 <- data.table(dt2 %>% 
      mutate(gene_comp = case_when(classified.x=="shared_poly" & 
                                   classified.y=="shared_poly" ~ "TSP",
                                   TRUE ~ "Not TSP")))

# Summarize
dt.fin <- data.table(dt3 %>% 
  group_by(dist, cont, gene_comp) %>% 
  summarise(mean = mean(R2, na.rm = T),
            uci.r2 = quantile(R2, probs = 0.95, na.rm = T),
            lci.r2 = quantile(R2, probs = 0.05, na.rm = T),
            median = median(R2, na.rm = T)))

# saveRDS(dt.fin, "/scratch/csm6hg/data/linkage.r2genes.filt.conservative.rds")
dt <- data.table(readRDS("/project/berglandlab/connor/data/linkage.r2genes.filt.rds"))

# Plot LD decay
require(scales)
ld.cont <- {dt.fin[dist>0][gene_comp=="TSP"] %>% 
  ggplot(., 
     aes(x=dist, 
         y=mean, 
         color=gene_comp)) + 
  facet_wrap(~cont, nrow = 2) +
  geom_line(size=2, alpha=0.85) +
  geom_line(data=dt.fin[dist>0][gene_comp=="Not TSP"], 
            aes(x=dist, y=mean, color=gene_comp), size=2, alpha=0.85) +
  theme_bw() +
  annotation_logticks(scaled = TRUE) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(1, 2000)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.005, 1)) +
  labs(x = "Distance (bp)", 
       color = "",
       y = "Average R2") +
  scale_color_manual(values = c("TSP"="darkgreen", "Not TSP"="orange")) +
  theme(strip.text = element_text(face="bold.italic", size=20),
        legend.text = element_text(size=20), 
        legend.position = "bottom",
        legend.title = element_text(face="bold", size=20),
        axis.text.x = element_text(face="bold", size=20),
        axis.text.y = element_text(face="bold", size=20),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))}

# Plot all LD matrices
pdf("figures/LD_candgene_decay_neg_below2kb_logline_new.pdf")
ld.cont
dev.off()

# LM - Europe
aov1 <- lm(mean ~ dist, data=dt.fin[dist>0][cont=="euro"][gene_comp=="TSP"])
aov2 <- lm(mean ~ dist, data=dt.fin[dist>0][cont=="euro"][gene_comp=="Not TSP"])

summary(aov1)
summary(aov2)

# NAm.
aov3 <- lm(mean ~ dist, data=dt.fin[dist>0][cont=="nam"][gene_comp=="TSP"])
aov4 <- lm(mean ~ dist, data=dt.fin[dist>0][cont=="nam"][gene_comp=="Not TSP"])

summary(aov3)
summary(aov4)


### Boxplot ###

# Summarize
dt.fin3 <- data.table(dt3 %>% 
    group_by(distc, cont, gene_comp, .id) %>% 
    summarise(mean = mean(R2, na.rm = T),
              median = median(R2, na.rm = T)))

# Add limits by binning
dt.fin4 <- data.table(dt.fin3 %>% 
        mutate(start = as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
               end = as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
               mid = start+((end-start)/2)))

# Plot LD decay
ld.cont <- {dt.fin4[mid<1000] %>% 
    ggplot(., 
           aes(x=as.factor(end/1000), 
               y=mean,
               group=interaction(as.factor(end/1000), gene_comp),
               fill=gene_comp)) + 
    facet_wrap(~cont, nrow = 2) +
    geom_boxplot(outlier.alpha = 0) +
    theme_bw() +
    labs(x = "Distance (kbp)", 
         fill = "",
         linetype = "Species",
         y = "Average Gene R2") +
    scale_fill_manual(values = c("TSP"="darkgreen", "Not TSP"="orange")) +
    theme(strip.text = element_text(face="bold.italic", size=20),
          legend.text = element_text(size=20), 
          legend.position = "bottom",
          legend.title = element_text(face="bold", size=20),
          axis.text.x = element_text(face="bold", size=20),
          axis.text.y = element_text(face="bold", size=20),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

### Wilcox test ###
ddd <- unique(dt.fin4$end)

# Wilcox test 
wily <- foreach(i=1:length(ddd), .combine = "rbind", .errorhandling = "remove") %do% {
  
  print(i)
  e <- wilcox.test(mean ~ gene_comp, 
              data = dt.fin4[end==ddd[i]],
              subset = cont=="euro") 
  n <- wilcox.test(mean ~ gene_comp, 
                  data = dt.fin4[end==ddd[i]],
                  subset = cont=="nam") 
  end <- data.table(rbind(data.table(stat=e$statistic, pval=e$p.value, cont="euro"),
                          data.table(stat=n$statistic, pval=n$p.value, cont="nam")),
                    i=i, end=ddd)
  
  return(end)
}

hist(wily[cont=="euro"][end<1000]$pval)
