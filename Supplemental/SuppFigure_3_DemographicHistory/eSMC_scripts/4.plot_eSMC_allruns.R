# Collect and plot eSMC output
# 11.20.2024
# Connor Murray
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
require(scales)

# Output
setwd('/project/berglandlab/connor/')

# Mutation rate
mu=5.69e-09

# Generations per year for msmc
gen=5

# Metadata
fin <- data.table(read.csv("backup_project/metadata/samples.fin.9.8.22.csv"))

# Get distribution for beta
beta <- rbind(data.table(readRDS("data/beta_eSMC_gen5.rds"), gen=5),
              data.table(readRDS("data/beta_eSMC_gen10.rds"), gen=10))

# Get eSMC results
esmc <- rbind(data.table(readRDS("data/eSMC_gen5.rds"), gen=5),
              data.table(readRDS("data/eSMC_gen10.rds"), gen=10))

# All eSMC files
files <- rbindlist(lapply(list.files(path = "eSMC_output_11_15_24_gen10*", pattern = "rds", full.names = T), 
                          readRDS), fill = T) %>%
  right_join(fin %>% 
               select(Sample, Species, Continent, cont))

# Output esmc
es <- na.omit(data.table(files %>%
              group_by(Sample, chromosome) %>%
              mutate(time_index = seq_along(x))))


#saveRDS(esmc, file = "../data/eSMC_gen5.rds")

# Summary stats of beta
confint(lm(b~1, data=beta[cont=="Daphnia.pulex.NorthAmerica"][gen==5]), level = 0.95)
confint(lm(b~1, data=beta[cont=="Daphnia.pulex.NorthAmerica"][gen==10]), level = 0.95)
mean(beta[cont=="Daphnia.pulex.NorthAmerica"][gen==5]$b)
mean(beta[cont=="Daphnia.pulex.NorthAmerica"][gen==10]$b)

confint(lm(b~1, data=beta[cont=="Daphnia.pulex.Europe"][gen==5]), level = 0.95)
confint(lm(b~1, data=beta[cont=="Daphnia.pulex.Europe"][gen==10]), level = 0.95)
mean(beta[cont=="Daphnia.pulex.Europe"][gen==5]$b)
mean(beta[cont=="Daphnia.pulex.Europe"][gen==10]$b)

# Output msmc
msmc_sep <- data.table(read.csv(file = "backup_project/chapter1/msmc/msmsc2.mlg.pulex.csv") %>%
                         mutate(x1 = left_time_boundary/mu*gen,
                                y1 = (1/lambda)/(2*mu)) %>%
                         group_by(time_index, cont) %>%
                         summarize(x1=mean(x1, na.rm = T),
                                   meany1=mean(y1, na.rm = T),
                                   uciy1=quantile(y1, na.rm = T, probs = 0.95),
                                   lciy1=quantile(y1, na.rm = T, probs = 0.05), 
                                   gen=gen) %>%
                         select(cont, x1, meany1, uciy1, lciy1, gen),
                       analysis="MSMC2")

ms <- data.table(read.csv(file = "backup_project/chapter1/msmc/msmsc2.mlg.pulex.csv") %>%
          mutate(x1 = left_time_boundary/mu*gen,
                 y1 = (1/lambda)/(2*mu)))

# Get intervals
confint(lm(y1~1, data = ms[cont=="Daphnia.pulex.Europe"][x1 <= 1e6 & x1 > 1e3]), level = 0.95)
confint(lm(y1~1, data = ms[cont=="Daphnia.pulex.NorthAmerica"][x1 <= 1e6 & x1 > 1e3]), level = 0.95)

confint(lm(y~1, data = es[cont=="Daphnia.pulex.Europe"][x <= 1e6 & x > 1e3]), level = 0.95)
confint(lm(y~1, data = es[cont=="Daphnia.pulex.NorthAmerica"][x <= 1e6 & x > 1e3]), level = 0.95)

# Make mega demo history object
mega <- data.table(rbind(na.omit(esmc[!cont %like% "pulexcaria"]) ,
                         msmc_sep[!cont %like% "pulexcaria"], fill=T))

# Through time MSMC
require(scales)
plot <- {
  mega %>%
    ggplot(., aes(x = x1,
                  y = meany1,
                  ymin = lciy1,
                  ymax = uciy1,
                  color = cont,
                  fill = cont)) +
    geom_ribbon(alpha=0.5) +
    geom_line(size=2.2) +
    facet_wrap(~analysis+gen, nrow = 1) +
    annotation_logticks(scaled = TRUE) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits=c(10000, 5e6)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits=c(1000, 2e7)) +
    theme_bw() +
    scale_color_manual(values=c("Daphnia.pulex.Europe"="#661100",
                                "Daphnia.pulex.NorthAmerica"="#0072B2")) +
    scale_fill_manual(values=c("Daphnia.pulex.Europe"="#661100",
                               "Daphnia.pulex.NorthAmerica"="#0072B2")) +
    coord_fixed() +
    labs(x="Years ago",
         y="Effective population size (Ne)",
         color="",
         linetype="") +
    theme(legend.text = element_text(face="bold", size=22),
          legend.title = element_text(face="bold", size=22),
          legend.background = element_blank(),
          strip.text = element_text(face="bold", size=18),
          legend.position= "bottom",
          axis.text.x = element_text(face="bold", size=22),
          axis.text.y = element_text(face="bold", size=22),
          axis.title.x = element_text(face="bold", size=24),
          axis.title.y = element_text(face="bold", size=24))
}

# Output plot
pdf("/project/berglandlab/connor/SuppFig3.esmc.pdf", width=8, height=6)
plot
dev.off()
