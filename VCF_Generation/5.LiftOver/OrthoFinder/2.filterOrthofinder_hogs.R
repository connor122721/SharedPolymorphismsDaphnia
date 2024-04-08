# Extract SCO genes across species
# Connor Murray 9.3.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

library(data.table)
library(tidyverse)

# WD
setwd("/scratch/csm6hg/genomes/proteins_species/primary_transcripts/OrthoFinder/Results_Aug31/")

# Single copy genes
system("cat Single_Copy_Orthologue_Sequences/* | grep ">" | cut -f2 -d">" | sort | uniq > unique.sco.genes.txt")

# SCO genes
ortho <- data.table(fread("unique.sco.genes.txt", header = F))
colnames(ortho) <- c("gene")

# Euro SCO
euro <- unique(ortho[gene %like% "Daphnia"]$gene)

# NAm SCO
nam <- unique(ortho[!gene %like% "Daphnia"]$gene)

# Exclude HOGs present in only 1 species
fwrite(x = data.table(euro), file = 'euro.sco.genes.txt', sep='\t', col.names = F)
fwrite(x = data.table(nam), file = 'nam.sco.genes.txt', sep='\t', col.names = F)

gff.nam <- data.table(fread("/scratch/csm6hg/genomes/pulex_nam/pulex.nam.genomic.gff", drop = 8))
