# Classifying clonal lineages
# Connor Murray 3.9.2022

# Run using these commands:
# ijob -A berglandlab_standard --mem=10G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(ggrepel)
library(poppr)
library(vcfR)
library(vegan)
library(adegenet)
library(viridis)
library(cowplot)
library(ggrepel)
library(parallel)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2")

# Metadata
fin.dt <- data.table(read.csv(file = "../metadata/samples.1.31.22.csv", header = T))

# Read in poppr output
files <- list.files(pattern = ".rds")

# With threshold applied
snpclone <- files[which(files %like% "bestthresh.snpclone"==TRUE)]

# Analyze output for each group
out <- foreach(i=1:length(snpclone), .combine = "rbind") %do% {
  
  # Progress message
  print(i)
  
  # Read in snpclone object with best threshold applied
  xsnp <- readRDS(snpclone[i])
  
  # Extract individuals and MLGs
  snp.dt <- data.table(data.table(Sample=indNames(xsnp),
                                  mlg=mll(xsnp)) %>% 
                         left_join(fin.dt %>% select(Sample, Species, Continent, country, SC), 
                                   by=c("Sample")))
  
  # Create unique MLG names
  snp.dt <- data.table(snp.dt %>% 
              group_by(country) %>% 
              mutate(mlg.country=paste(country, mlg, sep=".")))
  
  # Add populations to snpclone object
  pop(xsnp) <- paste(str_replace(snp.dt$Species,pattern = " ", replacement = "."), 
                      snp.dt$Continent, sep=".")
  
  # Exract mlgs
  dt <- data.table(n=t(mlg.table(xsnp, plot=F)),
                   mlg=colnames(mlg.table(xsnp, plot=F))) 
  
  # Wide to long reformat
  dt.fin <- data.table(dt %>% 
                         pivot_longer(cols = colnames(dt)[1], 
                                      names_prefix = "n.", 
                                      names_to="pop", 
                                      values_to="n"),
                       thresh.best=as.numeric(xsnp@mlg@cutoff[2]))
  
  # Finish
  return(snp.dt)
  
}

# Unique preidentified superclones from Europe
euro.sc <- unique(out[Continent=="Europe"][!SC==""]$SC)

# Compare to european clonal classifications
euro.valif <- foreach(i=1:length(euro.sc), .combine = "rbind") %do% {
  
  print(i)
  p <- out[Continent=="Europe"][SC==euro.sc[i]]
  p.list <- data.table(unique.mlg=length(unique(p$mlg)),
                       n.min.mlg=min(table(p$mlg)),
                       n.max.mlg=max(table(p$mlg)),
                       n.mean.mlg=mean(table(p$mlg)),
                       prop_correct=(length(unique(p$mlg))/length(p$SC))*100,
                       SC=euro.sc[i],
                       SC.size=length(p$SC),
                       iteration=i)
  
  return(p.list)
  
}

a <- euro.valif %>% 
  group_by(SC) %>% 
  summarise(prop.right=n.max.mlg/SC.size,
            unique.mlg,
            SC.size) %>% 
  filter(SC.size>1, unique.mlg>1) %>% 
  ggplot(aes(x=reorder(SC, -SC.size), y=prop.right*100)) +
  geom_col(aes(fill=log10(unique.mlg))) +
  geom_label(aes(y=0, label=SC.size)) +
  scale_fill_viridis(option = "inferno") + 
  labs(x="Super clone", 
       y="Proportion of samples within major MLG (%)",
       fill="log10(Number of MLGs)") +
  theme_classic() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.background = element_blank(),
        axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

b <- euro.valif %>% 
  group_by(SC) %>% 
  summarise(prop.right=n.max.mlg/SC.size,
            unique.mlg,
            SC.size) %>% 
  filter(SC.size>1, unique.mlg>1) %>% 
  ggplot(aes(x="", y=prop.right*100, label=SC)) +
  geom_boxplot() +
  geom_label_repel(aes(fill=log10(unique.mlg)), color="white", max.overlaps=20) +
  scale_fill_viridis(option = "inferno") + 
  ylim(0,100) +
  labs(x="Super clone", 
       y="",
       fill="log10(Number of MLGs)") +
  theme_classic() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=15),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

png("poppr.SC.wholegenome.png")

# Combine
plot_grid(a, b, rel_widths = c(2,1))

dev.off()

# Merge metadata and MLG assignment
out.fin <- data.table(out %>% 
              select(Sample, mlg.country) %>% 
              left_join(fin.dt %>% select(-mlg.country)))

# Write MLG assignment file
write.csv(out.fin, 
          file = "../metadata/samples.9.8.22.csv",
          quote = FALSE,
          row.names = FALSE)

# Analyze clonal assignment
fin.dt <- data.table(read.csv("../metadata/samples.9.8.22.csv"))

require(scales)

pdf("../figures/numbermlgcontinent_new_filt.pdf", width=12, height=8)

fin.dt[!Species %like% "Daphnia.pulexcaria"] %>% 
  group_by(cont, Species, Continent) %>% 
  summarise(n=length(unique(mlg.country))) %>%
  ggplot(., aes(x=n, 
                y=reorder(str_replace_all(string = cont, pattern = "[.]", replacement = " "), n),
                #fill=Species
                )) +
  geom_col() +
  labs(x="N. multi-locus genotyoes", 
       y="") +
  theme_classic() +
  annotation_logticks(scaled = TRUE, sides = "b") +
  scale_x_log10() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold.italic", size=13),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=18),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold.italic", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()
