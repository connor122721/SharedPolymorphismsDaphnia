# Make Figure 5
# Connor Murray 10.26.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(viridis)
require(scales)
library(DescTools)
library(patchwork)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Read in LD decay (r^2 ~ bp) data
ld <- data.table(readRDS("/scratch/csm6hg/data/linkage.r2genes.filt.conservative.rds") %>% 
             mutate(cont=case_when(cont=="euro" ~ "Euro D. pulex",
                                   cont=="nam" ~ "NAm. D. pulex")))

# Read in exon diversity data
exon.cat <- data.table(na.omit(readRDS(file = "/scratch/csm6hg/data/exon.diversity.stats.conservative.rds") %>% 
              pivot_longer(names_to = "stats", cols = c(fst:taj))))

# Plot LD decay
ld.plot <- {ld[dist>0][gene_comp=="TSP"] %>% 
    ggplot(., 
           aes(x=dist, 
               y=mean, 
               color=gene_comp)) + 
    facet_wrap(~cont, nrow = 2) +
    geom_line(size=2, alpha=0.85) +
    geom_line(data=ld[dist>0][gene_comp=="Not TSP"], 
              aes(x=dist, y=mean, color=gene_comp), size=2, alpha=0.85) +
    theme_bw() +
    annotation_logticks(scaled = TRUE) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1, 2000)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = c(0.005, 1)) +
    labs(x = "Distance (bp)", 
         color = "",
         y = expression(~bold("Gene Average "~italic(r)^2))) +
    scale_color_manual(values = c("TSP"="turquoise", "Not TSP"="deeppink4")) +
    theme(strip.text = element_text(face="bold.italic", size=22),
          legend.text = element_text(size=20), 
          legend.position = c(0.3, 0.7),
          legend.background = element_blank(),
          legend.title = element_text(face="bold", size=22),
          axis.text.x = element_text(face="bold", size=22),
          axis.text.y = element_text(face="bold", size=22),
          axis.title.x = element_text(face="bold", size=22),
          axis.title.y = element_text(face="bold", size=22),
          axis.title = element_text(face="bold", size=22))}

# Genetic diversity at exons
dive.plot <- {exon.cat[stats=="pi"] %>% 
    ggplot(aes(x=value,
               y=exon.size.cat,
               fill=classified)) +
    facet_wrap(~cont, nrow = 2) +
    geom_boxplot() + 
    theme_bw() +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_manual(values = c("TSP"="turquoise", "Not TSP"="deeppink4")) +
    labs(x = expression(~bold("log"[10]~"(Average Exon Genetic Diversity)")),
         y = "",
         fill = "") +
    theme(strip.text = element_text(face="bold.italic", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}

# Tajima's D at exons
tajima.plot <- {exon.cat[stats=="taj"] %>% 
    ggplot(aes(x=value,
               y=exon.size.cat,
               fill=classified)) +
    facet_wrap(~cont, nrow = 2) +
    geom_boxplot() + 
    theme_bw() +
    scale_fill_manual(values = c("TSP"="turquoise", "Not TSP"="deeppink4")) +
    labs(x = "Average Exon Tajima's D",
         y = "",
         fill = "") +
    theme(strip.text = element_text(face="bold.italic", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}

# FST at exons
fst.plot <- {exon.cat[stats=="fst"][!cont=="Across species"][value>1e-7] %>% 
    ggplot(aes(x=value,
               y=exon.size.cat,
               fill=classified)) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    facet_wrap(~cont, nrow = 2) +
    geom_boxplot() + 
    theme_bw() +
    scale_fill_manual(values = c("TSP"="turquoise", "Not TSP"="deeppink4")) +
    labs(x = expression(~bold("log"[10]~"(Average Exon "~italic("F"[ST])~")")),
         y = "",
         fill = "") +
    theme(strip.text = element_text(face="bold.italic", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}

# Compile diversity plots
mega.diver.plot <- (dive.plot + tajima.plot + fst.plot + 
                  plot_layout(guides = 'collect') & 
                  theme(legend.position = 'bottom'))

# Output
ggsave(plot = mega.diver.plot, filename = "/scratch/csm6hg/figs/Figure5new.pdf", width=16, height=8)
ggsave(plot = ld.plot, filename = "/scratch/csm6hg/figs/ldplot_conservative.pdf", width=10, height=8)


### Statistics ###
detach("package:DescTools", unload=T)

#### Beta1 
# NAM
t.test(ld[dist>0][gene_comp=="TSP"])

#### TAJIMA'S D
# NAM
wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Small"][classified=="TSP"][stats=="taj"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Small"][classified=="Not TSP"][stats=="taj"]$value)

wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Int"][classified=="TSP"][stats=="taj"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Int"][classified=="Not TSP"][stats=="taj"]$value)

wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Lar"][classified=="TSP"][stats=="taj"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Lar"][classified=="Not TSP"][stats=="taj"]$value)

# Euro
wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Small"][classified=="TSP"][stats=="taj"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Small"][classified=="Not TSP"][stats=="taj"]$value)

wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Int"][classified=="TSP"][stats=="taj"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Int"][classified=="Not TSP"][stats=="taj"]$value)

wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Lar"][classified=="TSP"][stats=="taj"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Lar"][classified=="Not TSP"][stats=="taj"]$value)

#### FST
# NAM
wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Small"][classified=="TSP"][stats=="fst"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Small"][classified=="Not TSP"][stats=="fst"]$value)

wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Int"][classified=="TSP"][stats=="fst"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Int"][classified=="Not TSP"][stats=="fst"]$value)

wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Lar"][classified=="TSP"][stats=="fst"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Lar"][classified=="Not TSP"][stats=="fst"]$value)

# Euro
wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Small"][classified=="TSP"][stats=="fst"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Small"][classified=="Not TSP"][stats=="fst"]$value)

wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Int"][classified=="TSP"][stats=="fst"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Int"][classified=="Not TSP"][stats=="fst"]$value)

wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Lar"][classified=="TSP"][stats=="fst"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Lar"][classified=="Not TSP"][stats=="fst"]$value)

#### THETA PI
# NAM
wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Small"][classified=="TSP"][stats=="pi"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Small"][classified=="Not TSP"][stats=="pi"]$value)

wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Int"][classified=="TSP"][stats=="pi"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Int"][classified=="Not TSP"][stats=="pi"]$value)

wilcox.test(x = exon.cat[cont %like% "America"][exon.size.cat %like% "Lar"][classified=="TSP"][stats=="pi"]$value, 
            y= exon.cat[cont %like% "America"][exon.size.cat %like% "Lar"][classified=="Not TSP"][stats=="pi"]$value)

# Euro
wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Small"][classified=="TSP"][stats=="pi"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Small"][classified=="Not TSP"][stats=="pi"]$value)

wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Int"][classified=="TSP"][stats=="pi"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Int"][classified=="Not TSP"][stats=="pi"]$value)

wilcox.test(x = exon.cat[cont %like% "Euro"][exon.size.cat %like% "Lar"][classified=="TSP"][stats=="pi"]$value, 
            y= exon.cat[cont %like% "Euro"][exon.size.cat %like% "Lar"][classified=="Not TSP"][stats=="pi"]$value)
