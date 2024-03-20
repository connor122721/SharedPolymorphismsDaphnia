# Read in admixture output and visualize
# Connor Murray 11.2.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(tidyverse)
library(foreach)
library(colortools)
library(forcats)

# Working directory
setwd('/scratch/csm6hg/oldscratch/daphnia_phylo/admixture/data_sub')

# Metadata
meta <- data.table(fread("/project/berglandlab/connor/metadata/samples.fin.9.8.22.csv"))

# MSMC output files
file1 <- list.files(path = "data_sub",
                    pattern = "sub.out$",
                    recursive = TRUE)

# Go through each K
fin.dt <- foreach(i=2:23, .combine = "rbind", .errorhandling = "remove") %do% {
  
  # Progress message
  print(i)
  
  # Load in admixture output
  tbl <- data.table(read.table(paste("daphnia.new.sub.", 
                                     i, 
                                     ".Q", 
                                     sep="")))
  
  # Metadata from plink
  indTable <- data.table(read.table("../daphnia.new.sub.fam",
                         col.names = c("Sample", "rep", "v1", "v2", "v3", "v4")))
  
  # Merge with metadata
  dt <- data.table(cbind(tbl, indTable %>% 
                           select(Sample)) %>% 
                     left_join(meta %>% 
                                 select(Sample, country, cont, Species, Continent), 
                               by="Sample"))
  
  # Order by species
  ordered <- dt[order(dt$cont),]
  
  # Merge
  fin <- data.table(dt %>% 
                      mutate(id = row_number()) %>% 
                      gather(key = 'pop', 'value',V1:paste("V",i, sep="")) %>%
                      group_by(id) %>% 
                      mutate(likely_assignment = pop[which.max(value)],
                             assingment_prob = max(value)) %>% 
                      arrange(likely_assignment, desc(assingment_prob)) %>% 
                      ungroup() %>% 
                      mutate(id = forcats::fct_inorder(factor(id))),
                    k=i)
  
  # Finish
  return(fin)
  
}

# Output
write.csv(fin.dt, file = "/project/berglandlab/connor/data/admix.sub.data.k2_23.csv")

# Create labels
lab <- data.table(lab.cont = unique(paste(fin.dt$Species, fin.dt$Continent, sep=".")))

#pdf("admixture.new.sub.k8.pdf", width = 15, height = 8)

# Plot K of various admixtures
admix <- {fin.dt[k %in% c(9)] %>% 
  mutate(kk = paste("K=", k, sep="")) %>% 
  mutate(kk=factor(kk, 
                   levels = c("K=9"))) %>% 
  ggplot(., 
         aes(x = factor(id), 
             y = value, 
             fill = factor(pop))) +
  geom_col(size = 0.01) +
  facet_grid(kk~fct_inorder(paste(Species, Continent)), 
             switch = "x", 
             scales = "free") +
  theme_minimal() + 
  labs(x = "Samples", 
       y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_fill_brewer(palette ="Paired") +
  theme(panel.spacing.x = unit(0.1, "lines"),
        panel.grid = element_blank(),
        strip.text = element_text(face="bold.italic", 
                                  size=12, 
                                  angle=30),
        legend.text = element_text(size=16),
        legend.position = "none",
        legend.title = element_text(face="bold", size=18),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))}

dev.off()

# Cross validation score
err <- data.table(fread("daphnia.sub.cv.error"))

# Create CV-error plot decay
#pdf("cv.error.new.sub.pdf", width = 8, height = 8)

ggplot(err, 
       aes(x = V1, 
           y = V2, 
           group=1)) +
  geom_point(size=1.2) +
  annotate("rect",fill="steelblue", xmin = 7.7, xmax=8.3, ymin = 0.1, ymax = 0.55,
           alpha = 0.2) +
  geom_line() +
  theme_bw() + 
  labs(x = "k", 
       y = "Cross-validation error") +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face = "bold", size=20),
        axis.text.y = element_text(face = "bold", size=20),
        axis.title.x = element_text(face = "bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))

dev.off()