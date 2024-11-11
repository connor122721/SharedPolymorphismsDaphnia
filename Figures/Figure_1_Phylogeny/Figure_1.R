# Make Figure 1
# Connor Murray 3.29.2024
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(ggforce)
library(patchwork)
library(ggtree)
library(ape)
library(treeio)

# Working directory
setwd("/scratch/csm6hg/data/")

# Metadata
fin <- data.table(read.csv("samples.fin.9.8.22.csv"))

# Open PCA rds
ccm_pca <- readRDS("pca.daphnia.genome.filt.mlgsub.nooutgroup.rds")

# Merge by Sample
pca <- data.table(sample=ccm_pca$sample.id, PC1=ccm_pca$eigenvect[,1],
                  PC2=ccm_pca$eigenvect[,2], PC3=ccm_pca$eigenvect[,3],
                  PC4=ccm_pca$eigenvect[,4], PC5=ccm_pca$eigenvect[,5])

# Remove missing samples
pca <- data.table(pca[sample %in% fin$Sample])

# Merge PCA and metadata
pca <- data.table(merge(pca, fin,
                        by.x="sample", by.y="Sample") %>% 
                mutate(Species=case_when(Species == "Daphnia pulex" ~ "D. pulex",
                                         Species == "Daphnia pulicaria" ~ "D. pulicaria",
                                         Species == "Daphnia pulexcaria" ~ 
                                           "D. pulex x D. pulicaria Hybrids",
                                         Species == "Daphnia obtusa" ~ "D. obtusa"),
                       Continent=case_when(Continent == "NorthAmerica" ~ "North America",
                                           TRUE ~ Continent)))

# PC 1/2 Plot
pca.plot <- {
  
  pca %>%
  ggplot(., aes(x=PC1, 
                y=PC2, 
                fill=Continent, 
                shape=Species)) +
  geom_point(color="black",
             size=6, 
             alpha=0.8) +
  geom_mark_ellipse() +
  theme_minimal() + 
  labs(x=paste("PC1 (", round(ccm_pca$varprop[[1]], digits=3)*100, " %)", sep=""),
       y=paste("PC2 (", round(ccm_pca$varprop[[2]], digits=3)*100, " %)", sep=""),
       title="Genome-wide SNPs") +
  scale_fill_manual(values = c("Europe"="cyan",
                               "North America"="deeppink1"),
                    name = "Species complex") +
  theme(strip.text = element_text(face="bold.italic", size=16),
        title = element_text(size=18, face="bold"),
        legend.text = element_text(size=16, face="bold.italic"),
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))}

ggsave(pca.plot, filename = "figures/pca.new.outgroup.pdf", height = 10, width = 12)

# GLM
library(multcomp) 
pca$cont <- as.factor(pca$cont)
aov1 <- glm(PC1 + PC2 ~ cont, data=pca)
summary(glht(aov1, linfct = mcp(cont = "Tukey")))

### Map ###

# Metadata
fin <- data.table(read.csv(file = "metadata/samples.12.3.21.csv", header = T))

# Counts of samples
fin <- data.table(fin[!is.na(lat)] %>% 
                    group_by(Sample, Species, lat, long) %>% 
                    mutate(count=n()))

# Rename
fin.map <- data.table(fin %>% 
           group_by(Species, Continent, lat, long) %>% 
           summarise(count=n()) %>% 
           mutate(lat=as.numeric(lat),
                  long=as.numeric(long)) %>% 
           mutate(Species=case_when(Species=="Daphnia pulex" & Continent=="Europe" ~ "Euro D. pulex",
                  Species=="Daphnia pulicaria" & Continent=="Europe" ~ "Euro D. pulicaria",
                  Species=="Daphnia obtusa" & Continent=="Europe" ~ "Euro D. obtusa",
                  Species=="Daphnia pulex" & Continent=="NorthAmerica" ~ "NAm. D. pulex",
                  Species=="Daphnia pulicaria" & Continent=="NorthAmerica" ~ "NAm. D. pulicaria",
                  Species=="Daphnia obtusa" & Continent=="NorthAmerica" ~ "NAm. D. obtusa"))) 

# World map object
world <- map_data('world')

# Colors for species
c <- hcl.colors(n=3, "Oslo")
c1 <- hcl.colors(n=3, "Plasma")
colors1 <- c("Euro D. obtusa" = "red", 
             "Euro D. pulex" = c[2], 
             "Euro D. pulicaria" = "darkgreen")

colors2 <- c("NAm. D. obtusa" = c1[1], 
             "NAm. D. pulex" = c1[2], 
             "NAm. D. pulicaria" = c1[3])

# Plotting with scatter pie
euro <- {
  
  ggplot(world, aes(x=long, y=lat, group = group)) +   
    geom_shape(col=NA, fill = "tan") +
    coord_map("bonne", 
              lat0 = 0, 
              xlim = c(-11,3), 
              ylim=c(49.5, 59)) +
    geom_point(data=fin.map[Species=="Euro D. pulex"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    geom_point(data=fin.map[Species=="Euro D. pulicaria"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    geom_point(data=fin.map[Species=="Euro D. obtusa"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    scale_color_manual(values = colors1) +
    labs(color="European clade",
         x="Longitude",
         y="") +
    theme_void() + 
    theme(legend.position = "none",
          legend.background = element_rect(fill = "white", size = 0.5),
          legend.text = element_text(face = "bold.italic", size = 16),
          legend.title = element_text(face = "bold", size = 16),
          legend.box = "vertical",
          panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(colour = "black", fill=NA, size=3),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.ticks = element_blank(),
          panel.grid  = element_blank())

}

# Europe - CZE + LITH
euro2 <- {
  
  ggplot(world, aes(x=long, y=lat, group = group)) +   
    geom_shape(col=NA, fill = "tan") +
    coord_map("bonne", lat0=5, xlim = c(10, 28), ylim = c(40, 60)) + 
    geom_point(data=fin.map[Continent=="Europe"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    scale_color_manual(values = colors1) +
    labs(color="European clade",
         x="Longitude",
         y="") +
    theme_void() + 
    theme(legend.position = "none",
          legend.background = element_rect(fill = "white", size = 0.5),
          legend.text = element_text(face = "bold.italic", size = 16),
          legend.title = element_text(face = "bold", size = 16),
          legend.box = "vertical",
          panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(colour = "black", fill=NA, size=3),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.ticks = element_blank(),
          panel.grid  = element_blank())
  
}

# Plotting with scatter pie everything else
nam <- {
  
  ggplot(world, aes(x=long, y=lat, group = group)) +   
    geom_shape(fill="tan") +
    coord_map(projection = "bonne", 
              lat0 = 20,
              xlim = c(-120,-70), 
              ylim = c(25,60)) +
    geom_point(data=fin.map[Species=="NAm. D. pulex"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    geom_point(data=fin.map[Species=="NAm. D. pulicaria"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    geom_point(data=fin.map[Species=="NAm. D. obtusa"], 
               aes(x=long, y=lat, color=Species, group=1), alpha=0.9, size=10) +
    labs(color="North American clade",
         x="Longitude",
         y="Latitude") +
    theme_void() + 
    scale_color_manual(values = colors2) +
    theme(legend.position = "top",
          legend.background = element_rect(fill = "white", size = 0.5),
          legend.text = element_text(face = "bold.italic", size = 16),
          legend.title = element_text(face = "bold", size = 16),
          legend.box = "vertical", 
          panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(colour = "black", fill=NA, size=3),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.ticks = element_blank(),
          panel.grid  = element_blank())

}

# Compile map
map.plot <- ((nam + euro + euro2) + 
             plot_layout(guides = 'collect') & 
               theme(legend.position = 'top'))

### Phylogenetic tree ###

# Read in tree
myTree <- read.beast("snapp/snapp.tree.edit.nexus")

# Plot tree
tree <- {
  
  myTree %>% 
  ggtree(., size=1.5, 
         root.position = 1) +
  geom_range(range='length_0.95_HPD', color='cyan', size=3) +
  geom_nodelab(aes(label=round(posterior, 2)), nudge_x = 0.1, size=3) +
  geom_tiplab() +
  theme_tree2() +
  theme(text = element_text(size = 20, face = "bold.italic"), 
        legend.title = element_text(size=24, face="bold")) +
  geom_treescale(fontsize=6, linesize=2)
}

# Compile map
mega.plot <- ( map.plot / ( pca.plot + tree + plot_layout(widths = c(1, 1))))

# Output
pdf("figures/Figure1.pdf", width = 20, height = 16, )
mega.plot
dev.off()