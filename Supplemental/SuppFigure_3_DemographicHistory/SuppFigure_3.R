# Visualize output MSMC & SMC++ - SuppFigure 3
# Connor Murray 12.13.22
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
require(scales)

# SMC++ Output
setwd('/project/berglandlab/connor/msmc')

# Mutation rate
mu=5.69e-09

# Generations per year
gen=1

# Output smcpp
smcpp_sep <- data.table(read.csv("smcpp.phase.csv") %>% 
      select(Sample=label, cont, x1=x, y1=y, time) %>% 
      filter(!y1>8e6 & !y1<10000) %>% 
      group_by(time, cont) %>% 
      summarize(x1=mean(x1, na.rm = T),
                meany1=mean(y1, na.rm = T),
                uciy1=quantile(y1, na.rm = T, probs = 0.95),
                lciy1=quantile(y1, na.rm = T, probs = 0.05)) %>% 
      select(cont, x1, meany1, uciy1, lciy1),
      analysis="SMC++")

# Output smcpp
msmc_sep <- data.table(read.csv(file = "msmsc2.mlg.pulex.csv") %>% 
      mutate(x1 = left_time_boundary/mu*gen, 
             y1 = (1/lambda)/(2*mu)) %>% 
      group_by(time_index, cont) %>% 
      summarize(x1=mean(x1, na.rm = T),
                meany1=mean(y1, na.rm = T),
                sdy1=sd(y1, na.rm = T),
                dfy1=length(y1)) %>%
      mutate(sey1=sdy1/sqrt(dfy1)) %>% 
      mutate(lciy1=meany1-(1.96*sey1),
             uciy1=meany1+(1.96*sey1)) %>% 
      select(cont, x1, meany1, uciy1, lciy1),
      analysis="MSMC2")

# Through time MSMC
require(scales)
plot <- {ggplot() +
  geom_ribbon(data=msmc_sep[!cont %like% "pulexcaria"],
              aes(x = x1,
                  ymin = lciy1,
                  ymax = uciy1,
                  fill = cont,
                  linetype = analysis), 
              alpha=0.5) +
  geom_line(data=msmc_sep[!cont %like% "pulexcaria"],
            aes(x = x1,
                y = meany1,
                color = cont,
                linetype = analysis), size=2.2) +
  geom_ribbon(data=smcpp_sep[!cont %like% "pulexcaria"],
              aes(x = x1,
                  ymin = lciy1,
                  ymax = uciy1,
                  fill = cont,
                  linetype = analysis), 
              alpha=0.5) +
  geom_line(data=smcpp_sep[!cont %like% "pulexcaria"],
            aes(x = x1,
                y = meany1,
                color = cont,
                linetype = analysis), size=2.2) +
  facet_grid(~analysis) +
  annotation_logticks(scaled = TRUE) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(5000, 1e6)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(7000, 2e7)) +
  theme_bw() +
  scale_color_manual(values=c("Daphnia.pulex.Europe"="#661100", 
                              "Daphnia.pulex.NorthAmerica"="#0072B2")) +
  scale_fill_manual(values=c("Daphnia.pulex.Europe"="#661100", 
                              "Daphnia.pulex.NorthAmerica"="#0072B2")) +
  labs(x=paste("Years ago (g=", gen, ", u=", mu, ")", sep=""), 
       y="Effective population size (Ne)",
       color="",
       linetype="") +
  theme(legend.text = element_text(face="bold", size=22),
        legend.title = element_text(face="bold", size=22),
        legend.background = element_blank(),
        strip.text = element_text(face="bold.italic", size=22),
        legend.position= "bottom",
        axis.text.x = element_text(face="bold", size=22),
        axis.text.y = element_text(face="bold", size=22),
        axis.title.x = element_text(face="bold", size=24),
        axis.title.y = element_text(face="bold", size=24))}

# Output plot
pdf("/project/berglandlab/connor/figures/SuppFigure_2.pdf", width=14, height=10)
plot
dev.off()

# We use MSMC2 and remove most recent and very ancient estimates
msc <- data.table(read.csv(file = "msmsc2.mlg.pulex.csv") %>% 
             mutate(x1 = left_time_boundary/mu*gen, 
                    y1 = (1/lambda)/(2*mu)) %>% 
             group_by(time_index, cont) %>% 
             summarize(x1=mean(x1, na.rm = T),
                       meany1=mean(y1, na.rm = T),
                       sdy1=sd(y1, na.rm = T),
                       dfy1=length(y1)) %>% 
             filter(time_index>0 & time_index<25))

# Calculate the mean and standard error
l.model.nam <- lm(meany1 ~ 1, msc[cont=="Daphnia.pulex.NorthAmerica"])
l.model.euro <- lm(meany1 ~ 1, msc[cont=="Daphnia.pulex.Europe"])

# Calculate the confidence interval
confint(l.model.nam, level=0.95)
confint(l.model.euro, level=0.95)
