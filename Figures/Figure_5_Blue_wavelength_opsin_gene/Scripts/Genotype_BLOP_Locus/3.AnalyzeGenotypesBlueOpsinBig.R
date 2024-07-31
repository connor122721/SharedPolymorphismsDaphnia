# Trans species polymorphism - genotypes of blue wavelength opsin gene
# 10.1.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(ggtree)
library(adegenet)
library(cluster)
library(ape)
library(foreach)

# Working directory
setwd("C:/Users/Conno/Desktop/gitmaster/Blue_opsin_gene/Full_Europe/")

# Metadata
fin <- data.table(fread("Euro_samples_meta_7_27_23.csv"))

# Read tree
tre.n <- read.tree("Scaffold_1931_HRSCAF_2197.6350219.6351796.treefile")

# Filter tree to Euro pulex
tips <- data.table(sampleid=tre.n$tip.label) %>% 
  mutate(Sample=str_replace(sampleid, 
                            pattern = "...Scaffold_1931_HRSCAF_2197.6350219.6351796",
                            replacement = "")) %>% 
  mutate(Sample=str_replace(Sample, 
                            pattern = "NorthAmerica_", 
                            replacement = "")) %>% 
  mutate(Sample=str_replace(Sample, 
                            pattern = "Europe_", 
                            replacement = "")) %>% 
  mutate(Sample=gsub("^(([^_]+_){8})(.*)", "\\3", Sample)) %>% 
  left_join(fin, fill=T, by=c("Sample"="Sample_Name")) %>% 
  mutate(haplotype=case_when(sampleid %like% ".1.Scaffold"~1,
                             sampleid %like% ".2.Scaffold"~2))

tips <- tips[Organism=="Daphnia pulex" | is.na(Organism)]

# Read in aligned haplotypes
msa <- fasta2genlight("Scaffold_1931_HRSCAF_2197.6350219.6351796.aln.clip.fa",
                      parallel = F)

# Filter MSA to pulex species
msa.pul <- msa[msa@ind.names %in% tips[Organism == "Daphnia pulex" | is.na(Organism)]$sampleid]

# Extract samples and haplotypes
msa.dt <- data.frame(tab(msa.pul), 
                     sampleid = rownames(tab(msa.pul))) %>% 
  left_join(tips, by=c('sampleid'))

msa.dt2 <- data.table(msa.dt %>%
                    group_by(sampleid) %>% 
                    distinct(haplotype, Sample, Organism))

###### Cluster Haplotypes ######

# Perform K-Means clustering
n_clusters=10

# Look over 1 to n possible clusters
wss <- foreach(i=1:n_clusters, .combine = "rbind") %do% {
  # Fit the model: km.out
  print(i)
  km.out <- kmeans(msa.pul, centers = i)
  # Save the within cluster sum of squares
  dtii <- data.table(wsse=km.out$tot.withinss,
                     clusters=i)
  return(dtii)
}

# Shows decrease in performance of k
scree_plot <- {
  ggplot(data=wss, 
         aes(x = clusters, y = wsse, group = 1)) +
  geom_point(size = 4) +
  geom_line() +
  geom_hline(aes(yintercept = wsse), linetype=2) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
}

ggsave(plot = scree_plot, "skree.blop.europulex.kmeans3.pdf", width=8, height=8)

# From gene tree we can see 3-4 predefined groups
num_clusters <- 3
set.seed(123)
kmean <- kmeans(msa.pul, centers = num_clusters)
kmean.dt <- data.table(clust=kmean$cluster, 
                       sampleid = colnames(t(kmean$cluster))) %>% 
  left_join(msa.dt2, by=c("sampleid"))

# Remove duplicates
samps2rm <- kmean.dt[Sample %in% names(which(table(kmean.dt$Sample)>2))][sampleid %like% "index"]$sampleid
kmean.dt <- kmean.dt[!sampleid %in% samps2rm]

# Number of times a sample has a cluster
hap <- data.table(kmean.dt %>% 
            select(sampleid, Sample, haplotype, clust) %>% 
            mutate(haplotype=paste("hap_", haplotype, sep="")) %>% 
            mutate(Sample.row=paste(Sample,row_number(), sep="_")) %>% 
            pivot_wider(id_cols=Sample, 
                        names_from=haplotype, 
                        values_from=clust) %>% 
            mutate(Genotype = ifelse(hap_1 < hap_2, 
                        paste(hap_1, hap_2, sep="|"), 
                        paste(hap_2, hap_1, sep="|"))) %>% 
            mutate(hap.geno=case_when(hap_1==hap_2 ~ "Homo",
                                      !hap_1==hap_2 ~ "Alt")) %>% 
            left_join(fin %>% 
                        select(Sample_Name, Pond_ID), by=c("Sample"="Sample_Name")))

# Output data
write.table(hap, file = "euroPulex.haplo.clustk3.bluopsin.txt", quote = F, row.names = F)

# Only include samples in metadata
tips2keep = tre.n$tip.label[tre.n$tip.label %in% kmean.dt$sampleid]
tre.n <- keep.tip(phy = tre.n, tip = tips2keep)

### Gene tree ###

# Add cluster info to tips metadata
tips2 <- tips %>% left_join(kmean.dt, by=c("sampleid"))

# Tree plot
g <- {
    ggtree(tre.n, 
           size=1.7, 
           root.position = 1) %<+% 
    tips2 +
    geom_tippoint(aes(color = as.factor(clust)), 
                  size=3.5) +
    #geom_tiplab(aes(label=paste(Sample,haplotype,sep="_")), vjust=-.5, size=4, fontface=2) +
    labs(color="Haplotype cluster (k=3)", title="Blue wavelength opsin (Daphnia11806)") +
    #geom_treescale(x = 1.04, linesize = 1.5, fontsize = 6) +
    theme_tree(title = element_text(face="bold.italic", size=20),
               legend.text = element_text(face="bold.italic", size=16),
               legend.background = element_blank())
}

g.samp <- {
  ggtree(tre.n, 
         size=1.7, 
         root.position = 1) %<+% 
    data.table(tips %>% left_join(hap)) +
    geom_tippoint(aes(color = Genotype), 
                  size=3.5) +
    #geom_tiplab(aes(label=paste(Sample,haplotype,sep="_")), vjust=-.5, size=4, fontface=2) +
    labs(color="Genotype", title="Blue wavelength opsin (Daphnia11806)") +
    #geom_treescale(x = 1.04, linesize = 1.5, fontsize = 6) +
    theme_tree(title = element_text(face="bold.italic", size=20),
               legend.text = element_text(face="bold.italic", size=16),
               legend.background = element_blank())
}

# Output
ggsave("Daphnia11806-RA.tree.haplotype.pdf", dpi=300, g, width = 12, height=7)
ggsave("Daphnia11806-RA.tree.haplotype.sampsgeno.pdf", dpi=300, g.samp, width = 12, height=7)

# Tree + MSA plot
rownames(msa.dt) <- msa.dt$sampleid
pp <- gheatmap(g, 
               msa.dt[-c(131:167)],
               colnames = F, 
               width=1,
               high="purple",
               low="orange")

ggsave("Daphnia11806-RA.tree.haplotype.msa.png", 
       dpi = 300, pp, width = 12, height=7)
