# Compile and output dsuite data 
# 5/19/2022
# Connor Murray
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)
library(viridis)
library(cowplot)

# vcftools QC output
setwd('/project/berglandlab/connor/dsuite4/data1')

# Extract all files
miss <- list.files(pattern = "*Dmin.txt$")
mylist <- lapply(miss, fread)
miss <- data.table(rbindlist(mylist))
colnames(miss)[5:7] <- c("Z", "p", "f4")

# Add false discovery rate correction
miss <- data.table(miss %>% 
            mutate(fdr.p=p.adjust(p, method = "fdr")))

# Aggregate statistics
miss1 <- data.table(miss %>% 
  filter(P1 %like% c("Daphnia.pulex.Europe") &
         P2 %like% c("Daphnia.pulicaria.Europe") &
         fdr.p < 0.05 &
         f4 > 0.05) %>% 
  mutate(P1_cont=paste(tstrsplit(P1, ".", fixed=T)[[1]],
                       tstrsplit(P1, ".", fixed=T)[[2]],
                       tstrsplit(P1, ".", fixed=T)[[3]]),
         P2_cont=paste(tstrsplit(P2, ".", fixed=T)[[1]],
                       tstrsplit(P2, ".", fixed=T)[[2]],
                       tstrsplit(P2, ".", fixed=T)[[3]]),
         P3_cont=paste(tstrsplit(P3, ".", fixed=T)[[1]],
                       tstrsplit(P3, ".", fixed=T)[[2]],
                       tstrsplit(P3, ".", fixed=T)[[3]])) %>%
  mutate(cont_comp=paste(P1_cont, P2_cont, P3_cont)) %>% 
  group_by(cont_comp, P1, P2, P3) %>% 
  summarize(d.mean=mean(Dstatistic),
            d.uci=quantile(Dstatistic, probs = 0.95),
            d.lci=quantile(Dstatistic, probs = 0.05),
            f.mean=mean(f4),
            f.uci=quantile(f4, probs = 0.95),
            f.lci=quantile(f4, probs = 0.05))) 

miss2 <- data.table(miss %>% 
             filter(P1 %like% c("Daphnia.pulex.NorthAmerica") &
                    P2 %like% c("Daphnia.pulicaria.NorthAmerica") &
                    fdr.p < 0.05 &
                    f4 > 0.05) %>% 
             mutate(P1_cont=paste(tstrsplit(P1, ".", fixed=T)[[1]],
                                           tstrsplit(P1, ".", fixed=T)[[2]],
                                           tstrsplit(P1, ".", fixed=T)[[3]]),
                             P2_cont=paste(tstrsplit(P2, ".", fixed=T)[[1]],
                                           tstrsplit(P2, ".", fixed=T)[[2]],
                                           tstrsplit(P2, ".", fixed=T)[[3]]),
                             P3_cont=paste(tstrsplit(P3, ".", fixed=T)[[1]],
                                           tstrsplit(P3, ".", fixed=T)[[2]],
                                           tstrsplit(P3, ".", fixed=T)[[3]])) %>%
                      mutate(cont_comp=paste(P1_cont, P2_cont, P3_cont)) %>% 
                      group_by(cont_comp, P1, P2, P3) %>% 
                      summarize(d.mean=mean(Dstatistic),
                                d.uci=quantile(Dstatistic, probs = 0.95),
                                d.lci=quantile(Dstatistic, probs = 0.05),
                                f.mean=mean(f4),
                                f.uci=quantile(f4, probs = 0.95),
                                f.lci=quantile(f4, probs = 0.05))) 

miss3 <- data.table(miss %>% 
                      filter(P1 %like% c("Daphnia.pulex.Europe") & 
                             P3 %like% c("Daphnia.pulicaria.Europe")) %>% 
                      mutate(P1_cont=paste(tstrsplit(P1, ".", fixed=T)[[1]],
                                           tstrsplit(P1, ".", fixed=T)[[2]],
                                           tstrsplit(P1, ".", fixed=T)[[3]]),
                             P2_cont=paste(tstrsplit(P2, ".", fixed=T)[[1]],
                                           tstrsplit(P2, ".", fixed=T)[[2]],
                                           tstrsplit(P2, ".", fixed=T)[[3]]),
                             P3_cont=paste(tstrsplit(P3, ".", fixed=T)[[1]],
                                           tstrsplit(P3, ".", fixed=T)[[2]],
                                           tstrsplit(P3, ".", fixed=T)[[3]])) %>%
                      mutate(cont_comp=paste(P1_cont, P2_cont, P3_cont)) %>% 
                      group_by(cont_comp, P1, P2, P3) %>% 
                      summarize(d.mean=mean(Dstatistic),
                                d.uci=quantile(Dstatistic, probs = 0.95),
                                d.lci=quantile(Dstatistic, probs = 0.05),
                                f.mean=mean(f4),
                                f.uci=quantile(f4, probs = 0.95),
                                f.lci=quantile(f4, probs = 0.05))) 

# dstatistics + f4-ratio plot
dstat1 <- {miss %>% 
  ggplot(.) + 
  geom_raster(
    aes(x=P3,
        y=P2,
        fill=Dstatistic)) +
  scale_fill_viridis(option="magma") +
  labs(x="Target", y="", color="Admixture statistic:") + 
  coord_flip() + 
  theme_classic() +
  #ylim(c(0,0.47)) +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.position = "bottom",
        legend.text = element_text(face="bold.italic",size=16),
        legend.title = element_text(face="bold", size=18),
        title = element_text(face="bold.italic", size=20),
        axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold.italic", size=16),
        axis.title.x = element_text(face="bold.italic", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))}

dstat2 <- {miss2 %>% 
    ggplot(.) + 
    geom_pointrange(
      aes(x=P2,
          y=d.mean,
          ymin=d.lci,
          ymax=d.uci,
          color="D-statistic"), 
      position = position_jitter(height=0,width=0.3),
      size=2) + 
    geom_pointrange(
      aes(x=P2,
          y=f.mean,
          ymin=f.lci,
          ymax=f.uci,
          color="f4-ratio"), 
      position = position_jitter(height=0,width=0.3),
      size=2) +
    scale_color_manual(values=c("f4-ratio"="#E69F00", "D-statistic"="#0072B2")) +
    labs(x="", y="", color="Admixture statistic:") + 
    coord_flip() + 
    theme_classic() +
    ylim(c(0,0.47)) +
    theme(strip.text = element_text(face="bold.italic", size=16),
          legend.position = "bottom",
          legend.text = element_text(face="bold.italic",size=16),
          legend.title = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          axis.text.x = element_text(face="bold", size=16),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold.italic", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}

dstat3 <- {miss3 %>% 
    ggplot(.) + 
    geom_pointrange(
      aes(x=P2,
          y=d.mean,
          ymin=d.lci,
          ymax=d.uci,
          color="D-statistic"), 
      position = position_jitter(height=0,width=0.3),
      size=2) + 
    geom_pointrange(
      aes(x=P2,
          y=f.mean,
          ymin=f.lci,
          ymax=f.uci,
          color="f4-ratio"), 
      position = position_jitter(height=0,width=0.3),
      size=2) +
    scale_color_manual(values=c("f4-ratio"="#E69F00", "D-statistic"="#0072B2")) +
    labs(x="", y="", color="Admixture statistic:") + 
    coord_flip() + 
    theme_classic() +
    ylim(c(0,0.47)) +
    theme(strip.text = element_text(face="bold.italic", size=16),
          legend.position = "bottom",
          legend.text = element_text(face="bold.italic",size=16),
          legend.title = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          axis.text.x = element_text(face="bold", size=16),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold.italic", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}

pdf("../../figures/dstat_pulex_2.pdf", width=15, height=10)
plot_grid(dstat1, dstat2, nrow=1, rel_widths = c(2,1))
dev.off()