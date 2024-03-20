# Bam flagstats
# Connor Murray 6.27.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(fuzzyjoin)

# Working directory
setwd("/project/berglandlab/connor/")

# Sample metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Remove dorthe's extra samples

# Read in SNP LD info
files <- system("ls -f -R /project/berglandlab/connor/BACKUP_scratch/all_bam/bamqcout/*.stats", intern = TRUE)
listi <- lapply(files, fread, skip=0, select=c(1,4), header=F, fill=T)
setattr(listi, 'names', list.files(path = "/project/berglandlab/connor/BACKUP_scratch/all_bam/bamqcout/", pattern = ".stats"))

# Bind list 
dt <- data.table(rbindlist(listi, use.names = T, idcol = T) %>% 
          mutate(Sample=str_remove_all(str_remove_all(tstrsplit(.id, ".", fixed=T)[[1]], "_finalmap_RG"), "_finalmap_mdup")) %>% 
          filter(Sample %in% unique(fin$Sample),
                 V4 %in% c("mapped", "properly")) %>% 
          pivot_wider(values_from = "V1", names_from = "V4") %>% 
          group_by(Sample) %>% 
          mutate(prop_paired=100*(properly/mapped)) %>% 
          left_join(fin %>% select(Sample, cont, Origin), by=c("Sample")))

# Mean/median of mapped and properly paired
dt1 <- dt %>% 
  group_by(cont) %>% 
  summarize(mean_map=mean(mapped),
            med_map=median(mapped),
            mean_properly=mean(properly),
            med_properly=median(properly))

# Pulex
# Calculate the mean and standard error
l.model.nam <- lm(mapped ~ 1, dt[cont=="Daphnia.pulex.NorthAmerica"])
l.model.euro <- lm(mapped ~ 1, dt[cont=="Daphnia.pulex.Europe"])
l.model.nam2 <- lm(properly ~ 1, dt[cont=="Daphnia.pulex.NorthAmerica"])
l.model.euro2 <- lm(properly ~ 1, dt[cont=="Daphnia.pulex.Europe"])

# Calculate the confidence interval
confint(l.model.nam, level=0.95)
confint(l.model.euro, level=0.95)

confint(l.model.nam2, level=0.95)
confint(l.model.euro2, level=0.95)

# Pulicaria
l.model.nampul <- lm(mapped ~ 1, dt[cont=="Daphnia.pulicaria.NorthAmerica"])
l.model.europul <- lm(mapped ~ 1, dt[cont=="Daphnia.pulicaria.Europe"])
l.model.nampul2 <- lm(properly ~ 1, dt[cont=="Daphnia.pulicaria.NorthAmerica"])
l.model.europul2 <- lm(properly ~ 1, dt[cont=="Daphnia.pulicaria.Europe"])

confint(l.model.nampul, level=0.95)
confint(l.model.europul, level=0.95)

confint(l.model.nampul2, level=0.95)
confint(l.model.europul2, level=0.95)

# Pulex x Pulicaria hybrids
l.model.namhyb <- lm(mapped ~ 1, dt[cont=="Daphnia.pulexcaria.NorthAmerica"])
l.model.namhyb2 <- lm(properly ~ 1, dt[cont=="Daphnia.pulexcaria.NorthAmerica"])

confint(l.model.namhyb, level=0.95)
confint(l.model.namhyb2, level=0.95)

# Obtusa
l.model.namob <- lm(mapped ~ 1, dt[cont=="Daphnia.obtusa.NorthAmerica"])
l.model.euroob <- lm(mapped ~ 1, dt[cont=="Daphnia.obtusa.Europe"])
l.model.namob2 <- lm(properly ~ 1, dt[cont=="Daphnia.obtusa.NorthAmerica"])
l.model.euroob2 <- lm(properly ~ 1, dt[cont=="Daphnia.obtusa.Europe"])

confint(l.model.namob, level=0.95)
confint(l.model.euroob, level=0.95)

confint(l.model.namob2, level=0.95)
confint(l.model.euroob2, level=0.95)

# Plot
dt %>% 
  ggplot(aes(x=mapped,
             y=cont, 
             color=cont)) +
  geom_jitter()
