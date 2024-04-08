# Compile and output dsuite data 
# 5/25/2022
# Connor Murray
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)
library(viridis)
library(cowplot)

# vcftools QC output
setwd('/project/berglandlab/connor/dsuite4/data3')

# Extract all files
miss <- list.files(pattern = "*tree.txt$")
miss <- miss[!miss %like% "pulicaria"]
mylist <- lapply(miss, fread)
miss <- data.table(rbindlist(mylist), miss)
colnames(miss)[5:7] <- c("Z", "p", "f4")
miss <- data.table(miss %>% 
          mutate(start=as.integer(tstrsplit(miss, "_")[[3]]),
                 stop=as.integer(tstrsplit(miss, "_")[[4]])))

# Add false discovery rate correction
miss <- data.table(miss %>% 
            mutate(fdr.p=p.adjust(p, method = "fdr")))

pdf("../../figures/dsuite.pulicEuro.P2.15k.pdf", width = 12, height = 10)

miss[!P2=="Daphnia.pulicaria.Europe"] %>% 
  ggplot(., aes(x=(start+stop)/2,
                y=-log10(fdr.p),
                color=Dstatistic)) +
  facet_wrap(~P2, scales="free") +
  geom_point() +
  #geom_smooth(color="steelblue4",linetype=1) +
  #geom_line() +
  scale_color_viridis(option="inferno") +
  labs(x="Window start", y="-log10(p-value)", color="D-statistic") + 
  theme_bw() +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.position = "bottom",
        legend.text = element_text(face="bold",size=10),
        legend.title = element_text(face="bold.italic", size=18),
        title = element_text(face="bold.italic", size=20),
        axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold.italic", size=18),
        axis.title = element_text(face="bold", size=20))

dev.off()

# Aggregate statistics
miss1 <- data.table(miss %>% 
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
  group_by(P2, P3) %>% 
  summarize(d=max(Dstatistic),
            f=max(f4),
            p=min(fdr.p)))

# dstatistics + f4-ratio plot
dstat1 <- {miss1 %>% 
  filter(p<0.05) %>% 
  ggplot(.) + 
  geom_raster(
    aes(x=P2,
        y=P3,
        fill=f, 
        alpha=-log10(p)
        )) +
  scale_fill_viridis(option="inferno") +
  scale_x_discrete(limits=c("Daphnia.pulicaria.Europe",
                            "Daphnia.pulexcaria.NorthAmerica")) +
  scale_y_discrete(limits=c("Daphnia.pulex.Europe", 
                            "Daphnia.pulicaria.Europe",
                            "Daphnia.pulex.NorthAmerica",
                            "Daphnia.pulexcaria.NorthAmerica",
                            "Daphnia.pulicaria.NorthAmerica")) +
  labs(x="P2", y="P3", fill="f4-ratio", alpha="-log10(P-value)") + 
  theme_classic() +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.position = "top",
        legend.text = element_text(face="bold",size=14),
        legend.title = element_text(face="bold", size=18),
        title = element_text(face="bold.italic", size=20),
        axis.text.x = element_text(face="bold.italic", size=16),
        axis.text.y = element_text(face="bold.italic", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))}
dstat2 <- {miss1 %>% 
    filter(p <0.05) %>% 
    ggplot(.) + 
    geom_raster(
      aes(x=P2,
          y=P3,
          fill=d, 
          alpha=-log10(p)
      )) +
    scale_fill_viridis(option="viridis") +
    scale_x_discrete(limits=c("Daphnia.pulicaria.Europe",
                              "Daphnia.pulexcaria.NorthAmerica")) +
    scale_y_discrete(limits=c("Daphnia.pulex.Europe", 
                              "Daphnia.pulicaria.Europe",
                              "Daphnia.pulex.NorthAmerica",
                              "Daphnia.pulexcaria.NorthAmerica",
                              "Daphnia.pulicaria.NorthAmerica")) +
    labs(x="P2", y="P3", fill="D-statistic", alpha="-log10(P-value)") + 
    theme_classic() +
    theme(strip.text = element_text(face="bold.italic", size=16),
          legend.position = "top",
          legend.text = element_text(face="bold",size=14),
          legend.title = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          axis.text.x = element_text(face="bold.italic", size=16),
          axis.text.y = element_text(face="bold.italic", size=16),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))}

pdf("../../figures/dstat_pulex_raster_f4.pdf")
dstat1
dev.off()

pdf("../../figures/dstat_pulex_raster_d_.pdf")
dstat2
dev.off()