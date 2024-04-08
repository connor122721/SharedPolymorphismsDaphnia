# Compile and output vcftools QC data 
# This script will concatenate all QC output files from chromosomes per VCF file
# 2.11.2022
# Connor Murray

# Ran with this in command line:
# ijob -A berglandlab_standard -p standard --mem=50G -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)

# vcftools QC output
setwd('/project/berglandlab/connor/mapping_stats/coverage_bams/')

# Metadata
meta <- data.table(read.csv("../../metadata/samples.1.31.22.csv"))

# Extract all vcftools files
miss <- fread("passed", header = F)

# Files not run
miss1 <- data.table(miss %>% mutate(Sample=tstrsplit(V1, ".coverage.10kbp.window")[[1]]))
meta[!Sample%in%miss1$Sample]

# Output missing individuals 
# write.table(meta[!Sample%in%miss1$Sample]$Sample, file = "../../metadata/coverage.samples.redo", quote = F, row.names = F, col.names = F)

# QC files
files <- list(miss$V1)
nam <- c("coverage")

# Nested loop to output files
foreach(p=1:length(files)) %do% {
  
  # Extract grouping
  file <- files[[p]]
  
  # Read in files as table
  out <- foreach(i=1:length(file), .combine="rbind", .errorhandling = "remove") %do% {
    
    # Read file
    fin <- data.table(fread(file=file[[i]]),
                      Sample=as.character(gsub(as.character(file[[i]]), 
                              pattern = ".coverage.10kbp.window", 
                              replacement = "")),
                      iteration=i)
    
      # Message  
      print(i)
      
      # Finish
      return(fin)
    
  }
  
  # Merge with metadata
  fin <- data.table(out %>% 
          left_join(meta %>% 
             select(Sample, Species, Continent, Origin), 
                        by = "Sample"))
  
  # Aggregate of Continent and Origin
  fin.dt <- data.table(fin[Species %in% c("Daphnia obtusa","Daphnia pulicaria","Daphnia pulex")][Sample %in% unique(meta$Sample)] %>% 
                         group_by(V1, V2, V3, Continent, Species) %>% 
                         select(c(V4, V7)) %>%  
                         summarise_all(list(mean = ~ mean(.),
                                            median = ~ median(.),
                                            upp.ci = ~ quantile(., probs = 0.975),
                                            low.ci = ~ quantile(., probs = 0.025))))
  
  # Output information
  write.table(fin.dt, file = paste("/project/berglandlab/connor/data/", 
                                   nam[p], ".final.csv", sep = ""), 
              sep = ",", col.names = F, row.names = F, quote = F, append = T)
  
}

# Coverage by bam
out <- data.table(read.csv("../../data/coverage.final.csv", header = F))

# Fix column names
colnames(out) <- c("Chrom", "start", "stop", 
                   "Continent", "Species", 
                   "mean.n.snps", "mean.miss", 
                   "med.n.snps", "med.miss", "uci.n.snps",
                   "uci.miss", "lci.n.snps", "lci.miss")

pdf("../../figures/coverage.10kpb.mean.new.species.pdf", width = 12, height = 8)

# Bam missingness per window
ggplot(out) +
  facet_wrap(~Species, 
             scales = "free_x", 
             nrow = 3) +
  geom_histogram(aes(x=1-mean.miss, 
                     fill=Continent), 
                     color="black",
                     position = "dodge", 
                     bins = 50) +
  geom_vline(xintercept = 0.1, 
             linetype=2, 
             size=1) +
  theme_bw() +
  labs(x="Average Missingness", 
       y="Counts") +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))

dev.off()

# Create bam files for missingness based off 10kbp windows - 10% cut off within Daphnia pulex
bed <- data.table(out[1-mean.miss>=0.10][Species=="Daphnia pulex"][,region:=paste(Chrom, start, stop, sep=".")])

# Correct missingness
bed <- data.table(bed %>% 
                    select(region, mean.miss, Chrom, start, stop) %>% 
                    mutate(mean.miss=1-mean.miss))

# Make distinct
bed <- data.table(distinct(bed %>% select(Chrom, start, stop) %>% 
                             mutate(index=as.character(paste(Chrom, ":", start, "-", stop, sep="")),
                                    Chrom=as.character(Chrom))))

# Write bed file
write.table(bed %>% select(Chrom, start, stop), 
            file="../../data/miss10.daphnia.pulex.final.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")

# Write index files
dt <- data.table(fread("../../data/miss10.daphnia.pulex.merged.RMoutHiCGM.final.bed"))

# Make distinct
dt <- data.table(distinct(dt %>% 
                            mutate(index=as.character(paste(V1, ":", V2, "-", V3, sep="")))))

# Write output
write.table(dt %>% select(index), 
            file="../../data/miss10.daphnia.pulex.merged.RMoutHiCGM.final.index.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
