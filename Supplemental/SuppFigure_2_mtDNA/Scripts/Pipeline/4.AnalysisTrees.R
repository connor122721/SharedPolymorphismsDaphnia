# Collect and visualize iqtree2 output 
# Connor Murray 12.13.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Packages
library(data.table)
library(ape)
library(foreach)
library(phylobase)
library(tidyverse)
library(ggtree)
library(colortools)
library(phylotools)
library(phytools)
library(thacklr)
library(treespace)
library(TreeTools)

# Working directory
setwd('/scratch/csm6hg/from_old_scratch/daphnia_phylo/mito/tree')

# Metadata
fin.dt <- data.table(read.csv("/project/berglandlab/connor/metadata/samples.fin.3.9.22.csv"))

# Add Daphnia magna outgroup
out.group <- data.table(Origin="SRA",
                         Species="Daphnia magna",
                         Continent="Europe",
                         cont="Daphnia.magna.Europe",
                         Sample="Daphnia_magna.")
fin.dt <- data.table(fin.dt %>% full_join(out.group))

# All samples used
file <- data.table(fread("../samples_mito_all", header = F))

# Load trees
filenames <- list.files(pattern="mtdna.12cds.sub.aln.treefile")

# Color wheel
color <- wheel(color = "tomato", num = 15)

# Forloop through each window genome tree
out <- foreach(i=1:length(filenames), .combine = "rbind") %do% {
  
  # Read tree
  tre.n <- read.tree(filenames[[i]])

  # Progress message
  print(paste("Tree number:", i, sep=" "))

  # Open tree data
  tree <- data.table(File=as.character(tre.n$tip.label))

  # Extract naming
  tree <- data.table(tree %>%
              mutate(Sample=tstrsplit(File,".", fixed=T)[[2]],
                     tree=filenames[[i]]))
  
  # Rename species
  tree[Sample %like% "magna"]$File <- "mtdna_D8_119.Daphnia_magna."
  tree[Sample %like% "Spring_2016_W6_6"]$File <- "mtdna_D8_119.Spring_2016_W6_6.1"
  tree[Sample %like% "Spring_2016_W1_1"]$File <- "mtdna_D8_119.Spring_2016_W1_1.1"
  tree[Sample %like% "Spring_2016_D10_10"]$File <- "mtdna_D8_119.Spring_2016_D10_10.5"
  tree[Sample %like% "magna"]$Sample <- "Daphnia_magna."
  tree[Sample %like% "Spring_2016_W6_6"]$Sample <- "Spring_2016_W6_6.1"
  tree[Sample %like% "Spring_2016_W1_1"]$Sample <- "Spring_2016_W1_1.1"
  tree[Sample %like% "Spring_2016_D10_10"]$Sample <- "Spring_2016_D10_10.5"
  
  # Any missing comparisons?
  missing <- tree$Sample[which(unique(tree$Sample) %in% unique(fin.dt$Sample)==F)]

  # Merge with metadata
  tree <- data.table(na.omit(tree %>%
                      left_join(fin.dt %>%
                         select(Sample, Continent, cont, Species, Origin) %>%
                         mutate(Sample=as.character(Sample)), by=c("Sample"))) %>%
                         mutate(sample.n=paste(cont, Sample, sep=".")))

  # Color tips by Origin & Continent & Species
  tips <- data.table(tree %>%
            mutate(color.con = case_when(Continent == "Europe" ~ color[1],
                   Continent == "NorthAmerica" ~ color[5]),
                   color.spp = case_when(Species == "Daphnia pulex" ~ color[1],
                            Species == "Daphnia pulicaria" ~ color[2],
                            Species == "Daphnia obtusa" ~ color[3],
                            Species == "Simocephalus" ~ "Black"),
                   color.sppcont = case_when(cont == "Daphnia.pulex.NorthAmerica" ~ color[1],
                          cont == "Daphnia.pulexcaria.NorthAmerica" ~ color[2],
                            cont == "Daphnia.pulex.Europe" ~ color[5],
                            cont == "Daphnia.pulicaria.NorthAmerica" ~ color[6],
                            cont == "Daphnia.pulicaria.Europe" ~ color[7],
                            cont == "Daphnia.obtusa.NorthAmerica" ~ color[8],
                            cont == "Daphnia.obtusa.Europe" ~ color[9],
                            cont == "Daphnia.magna.NorthAmerica" ~ color[10],
                            cont == "Daphnia.magna.Europe" ~ "Black",
                            cont == "Daphnia.longispina.NorthAmerica" ~ color[12],
                            cont == "Daphnia.longispina.Europe" ~ color[13],
                            cont == "Curvirostris.Europe" ~ color[14],
                            cont == "Daphnia.longiremus.Europe" ~ color[15],
                            cont == "Simocephalus.Europe" ~ "Black")))
  
  # Fix names
  tre.n$tip.label[tre.n$tip.label %in% tree$File ==F] <- c(
    "mtdna_D8_119.Spring_2016_W6_6.1",
    "mtdna_D8_119.Spring_2016_D10_10.5",
    "mtdna_D8_119.Spring_2016_W1_1.1")
  
  # Only include samples in metadata
  tre.n1 <- keep.tip(phy = tre.n, tip = unique(tips$File))

  # Root tree to Euro D. magna
  tre.root <- root(tre.n1, outgroup=as.character(tips[Species=="Daphnia magna"]$File)[1], 
                   resolve.root=TRUE)

  # Color by continent
  cols=c(data.table(File=tre.root$tip.label) %>%
                  left_join(tips %>%
                  dplyr::select(File, color.sppcont)) %>%
                  dplyr::select(color.sppcont) %>%
                  unlist())

  # Plot tree
  tree <- groupClade(tre.root, .node=c(40))
  p <- {tree %>% 
    ggtree(., size=1.5, 
           root.position = 1, 
           branch.length = "none",
           aes(linetype=group)) %<+% 
    tips +
    geom_tippoint(aes(color = Species), size=10) +
    geom_nodelab(aes(x=branch, label=label), vjust=-0.5, size=6) +
    theme_tree(legend.position='right') +
    theme(text = element_text(size = 20, face = "bold.italic"), 
          legend.title = element_text(size=24, face="bold")) +
    labs(color="Species Complex", linetype="Clade") +
    scale_color_manual(values = c("Black", "steelblue", "deeppink4", "darkgreen"),
                       breaks = c("Daphnia magna", "Daphnia pulex",
                                  "Daphnia pulicaria", "Daphnia obtusa")) +
    scale_linetype_manual(values = c(1,4), labels = c("North America","Europe")) +
    geom_treescale(fontsize=6, family = "arial", linesize=2)}
  
  # Save tree
  pdf("/project/berglandlab/connor/figures/mtdna.12cds.iqtree.norm.sub.clade.pdf", width = 15, height = 15)
  p
  dev.off()

  # Finish tree i
  return(co.sum)

}
