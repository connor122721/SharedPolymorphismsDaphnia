# Trikinetics experiment Daphnia 10.27/11.2/11.27.2023
# Analyzed: 12.11.2023-1.1.2024
# Madison Karram & Connor Murray

# Library
library(data.table)
library(tidyverse)
library(ggridges)
library(patchwork)
library(cowplot)

# working directory
setwd("C:/Users/Conno/Desktop/gitmaster/TransSpecificPolymorphisms_Daphnia/Figures/Figure_5_Blue_wavelength_opsin_gene/Scripts/Trikinetics_Experiment/")

# Combined monitor data to one table
dt <- data.table(read.csv(file = "data/agg.clean.full.trikinetics.csv"))

# Common theme element
themei <- (
  theme_bw() +
    theme(strip.text = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          title = element_text(face="bold", size=20),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20)))

# Summarize
dt1 <- data.table(dt %>% 
         group_by(clone, id, light_trt, genotype) %>% 
         summarize(mean=mean(meani, na.rm = T),
                   med=median(meani, na.rm = T),
                   act=mean(active, na.rm = T)))

# Order
dt$light_trt <- factor(dt$light_trt, levels = c('White', 'Blue', 'Dark'))

# Add ID column
dt$ID <- paste(dt$clone, dt$light_trt, dt$id, sep="_")

# Decay of activity over time
act1.decay <- {
  dt %>% 
    ggplot(aes(x=hour,
               y=meani,
               color=genotype)) +
    geom_smooth(alpha=0.4, method = "loess") +
    scale_color_manual(breaks = c("1|1", "1|2", "2|2"), 
                        values=c("darkgoldenrod", "darkcyan", "darkred")) +
    facet_wrap(~light_trt, 
               scales = "free_x") +
    labs(x="Experimental Hour",
         y="Mean Activity Level",
         fill="Clone") +
    themei
}

# Save output
ggsave(mega, 
       filename = "trikinetics.decay.pdf", 
       height = 10, 
       width=14)

# Mixed models
library(lme4)
library(lmerTest)

### Mean Activity Over Expeirmental Period ###

# Genotype x light on activity
m.b0 <- lmer(meani ~ 1 + hour + (1|ID) +(1|exp), data=dt[!light_trt=="Dark"])
m.b1 <- lmer(meani ~ genotype + hour + (1|ID) + (1|exp), data=dt[!light_trt=="Dark"])
m.b2 <- lmer(meani ~ genotype*light_trt + hour + (1|ID) + (1|exp), data=dt[!light_trt=="Dark"])

anova(m.b0, m.b1, m.b2)
summary(m.b2)

# PostHoc testing
posthoc <- difflsmeans(model = m.b2, 
                       test.effs="genotype*light_trt", 
                       ddf="Kenward-Roger")

# Adjust p-values
post <- data.table(test=rownames(posthoc[-c(1:4),]), 
                   posthoc[-c(1:4),]) %>% 
  mutate(p.adj=p.adjust(p = `Pr(>|t|)`, method = "holm"))

plot(posthoc, which = "genotype:light_trt")
