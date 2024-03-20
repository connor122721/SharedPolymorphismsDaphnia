# Representative individuals per country and species
# Connor Murray 11.1.2022 
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(tidyverse)
library(foreach)

# Working directory
setwd('/project/berglandlab/connor/snapp5/')

# Metadata
fin <- data.table(fread("../metadata/samples.fin.9.8.22.csv"))

# Set seed
set.seed(100)

# Collect 2 random species per continent
p.dt <- data.table(fin[mean.depth <= 50] %>%
                     group_by(Species, Continent, cont) %>% 
                     top_n(2, mean.depth) %>% 
                     summarize(Sample = sample(Sample, 2)))

# Write genome comparisions list
write.table(p.dt %>% select(Sample),
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            file ="highdep.2inds.species.continent.list")

ccm_pca <- readRDS("/project/berglandlab/connor/metadata/pca.filt.busco.rds")

# Skree plot
plot(ccm_pca$varprop[1:20]*100)

# Merge by Sample
pca <- data.table(sample=ccm_pca$sample.id, PC1=ccm_pca$eigenvect[,1],
                  PC2=ccm_pca$eigenvect[,2], PC3=ccm_pca$eigenvect[,3],
                  PC4=ccm_pca$eigenvect[,4], PC5=ccm_pca$eigenvect[,5])

# Remove missing samples
pca <- data.table(pca[sample %in% fin$Sample])

# Merge PCA and metadata
pca <- data.table(merge(pca, fin,
                        by.x="sample", by.y="Sample"))

# PC 1/2 Plot
pca[!sample %in% p.dt$Sample] %>%
  ggplot(., aes(x=PC1, y=PC2, color=Continent, shape=Species)) +
  geom_point(size=7, alpha=0.5) +
  geom_point(data=pca[sample%in%p.dt$Sample],
             aes(x=PC1, y=PC2, shape=Species),
             size=7, color="black") +
  theme_bw() +
  labs(x=paste("PC1 (", round(ccm_pca$varprop[[1]], digits=3)*100, " %)", sep=""),
       y=paste("PC2 (", round(ccm_pca$varprop[[2]], digits=3)*100, " %)", sep="")) +
  scale_fill_brewer(name = "Sp. complex", palette = "Set2") +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))

# Write genome comparisions list
write.table(p.dt %>% 
              mutate(species=paste(str_replace(Species, 
                                               pattern = " ", 
                                               replacement = "."), 
                                   Continent, sep=".")) %>% 
              select(species, individual=Sample),
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            file ="individuals.2inds.continent.withObtusa.txt")

# Create constraints file
lognormal(0,31.5,1) crown Daphnia.obtusa.Europe,Daphnia.obtusa.NorthAmerica,Daphnia.pulex.Europe,Daphnia.pulicaria.Europe,Daphnia.pulex.NorthAmerica,Daphnia.pulexcaria.NorthAmerica,Daphnia.pulicaria.NorthAmerica
