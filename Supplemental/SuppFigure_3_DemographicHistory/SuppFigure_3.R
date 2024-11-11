# Collect and plot Ne software output
# 11.5.2024
# Connor Murray 

# Libraries
library(data.table)
library(tidyverse)
require(scales)

# eSMC Output
setwd('/project/berglandlab/connor/daphnia_eSMC_output/')

# Read demographic history object
mega <- fread("../demo_history_Daphnia.txt")
 
# Ne over time multi-software outcome
nePlot <- {
  
  mega %>% 
    ggplot(., aes(x = x1, 
                  y = meany1,
                  ymin = lciy1,
                  ymax = uciy1,
                  color = cont,
                  fill = cont, 
                  linetype = analysis)) +
    geom_ribbon(alpha=0.5) +
    geom_line(size=2.2) +
    facet_wrap(~analysis, nrow = 1) +
    annotation_logticks(scaled = TRUE) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5000, 5e6)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = c(7000, 2e7)) +
    theme_bw() +
    scale_color_manual(values=c("Daphnia.pulex.Europe"="#661100", 
                                "Daphnia.pulex.NorthAmerica"="#0072B2")) +
    scale_fill_manual(values=c("Daphnia.pulex.Europe"="#661100", 
                               "Daphnia.pulex.NorthAmerica"="#0072B2")) +
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
ggsave(plot = nePlot, filename = "/project/berglandlab/connor/SuppFig3.new.pdf", width = 10, height =6)
