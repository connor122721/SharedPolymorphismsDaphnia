# Code to generate site frequency spectra from moments and calculate shared polymorphisms
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Setup
library(tidyverse)
library(data.table)

# Working directory
setwd("/scratch/csm6hg/phd/moments")

# Moments results
ld.dt <- rbind(data.table(readRDS("/scratch/csm6hg/phd/data/ld.dt.tot.20x20.rescale90.rds"), run="90"),
               data.table(readRDS("/scratch/csm6hg/phd/data/ld.dt.tot.20x20.rescale80.rds"), run="80"),
               data.table(readRDS("/scratch/csm6hg/phd/data/ld.dt.tot.20x20.rescale50.rds"), run="50"),
               data.table(readRDS("/scratch/csm6hg/phd/data/ld.dt.tot.20x20.rescale10.rds"), run="10"))

# Average values for SFS
ld.dt1 <- data.table(ld.dt %>% 
              group_by(row, col, model, vcf, size, run) %>% 
              summarize(count_mean = mean(count, na.rm = T),
                        count_med = median(count, na.rm = T)))

# Denominator of all alleles
denom = sum(ld.dt1[model == "empirical"][run=="90"]$count_mean)

# Get SAP per replicate
ld.dt.sap <- data.table(ld.dt %>% 
                filter(row!="V1" & col !=1) %>% 
                group_by(model, rep, size, run) %>% 
                summarize(sap = (sum(count)/denom)*100))

# Calculate the confidence interval
sap.model_Mig <- lm(sap ~ 1, ld.dt.sap[model=="split_mig"][size==20][run=="90"][sap<3 & sap>=0])
confint(sap.model_Mig, level=0.95)

sap.model_Mig <- lm(sap ~ 1, ld.dt.sap[model=="split_mig"][size==20][run=="80"][sap<3 & sap>=0])
confint(sap.model_Mig, level=0.95)

sap.model_Mig <- lm(sap ~ 1, ld.dt.sap[model=="split_mig"][size==20][run=="50"][sap<3 & sap>=0])
confint(sap.model_Mig, level=0.95)

sap.model_noMig <- lm(sap ~ 1, ld.dt.sap[model=="split_no_mig"][size==20][run=="90"][sap<3 & sap>=0])
confint(sap.model_noMig, level=0.95)

sap.model_noMig <- lm(sap ~ 1, ld.dt.sap[model=="split_no_mig"][size==20][run=="80"][sap<3 & sap>=0])
confint(sap.model_noMig, level=0.95)

sap.model_noMig <- lm(sap ~ 1, ld.dt.sap[model=="split_no_mig"][size==20][run=="50"][sap<3 & sap>=0])
confint(sap.model_noMig, level=0.95)

# Plotting SFS
SFS <- {
  ld.dt.sap[!model %in% c("from_ests", "empirical")][sap < 2 & sap >0] %>% 
    ggplot(., 
           aes(x=run, 
               y=sap, 
               color=model)) +
    geom_jitter() + 
    theme_classic() +
    theme(axis.text.x=element_text(face="bold", size=16),
          axis.text.y=element_text(face="bold", size=16),
          text = element_text(face="bold", size=18), 
          legend.text = element_text(face="bold", size=16),
          strip.text=element_text(face="bold", size=16)) 
  
}

# Merge list of summary stats
files.sum <- list.files(full.names = TRUE, pattern = '*.txt$', path = "../../output_eSMC/")
files.sum <- files.sum[!files.sum %like% "eSMC_gen10"]
sap <- lapply(files.sum, function(x){fread(x,header=F)})
names(sap) <- files.sum

# Merge
sap.dt <- data.table(rbindlist(sap, use.names = T, idcol = T, fill = T)) %>% 
  mutate(rep=tstrsplit(.id, ".", fixed=T)[[7]],
         vcf="Filtered")

colnames(sap.dt) <- c("file", "filt", "ns", "model",
                      "AIC", "BIC", "sap", "NNam", "NEuro", "s", "rep", "vcf")

# saveRDS(sap.dt, file = "/scratch/csm6hg/phd/data/sap.dt.tot.20x20.7000.2000.rds")
# sap.dt <- data.table(readRDS("/scratch/csm6hg/data/sap.dt.tot.20x20.rds"))

# saveRDS(ld.dt, file = "/scratch/csm6hg/data/ld.dt.tot.20x20.rds")
# ld.dt <- data.table(readRDS("/scratch/csm6hg/data/ld.dt.tot.20x20.rds"))

# Average SAP and get CI
sap.model_mig <- lm(sap ~ 1, sap.dt[model=="sfs_split_mig_model"][ns==20])
sap.model_emp <- lm(sap ~ 1, sap.dt[model=="sfs_empirical"][ns==20])
sap.model_noMig <- lm(sap ~ 1, sap.dt[model=="sfs_split_no_mig_model"][ns==20])

# Calculate the confidence interval
confint(sap.model_mig, level=0.95)
confint(sap.model_emp, level=0.95)
confint(sap.model_noMig, level=0.95)

sap.dt2 <- sap.dt %>% 
  group_by(model,ns) %>% 
  summarize(m=mean(sap))

# Function to calculate standardized residuals
calculate_standardized_residuals <- function(data) {
  
  #data=ld.dt1
  
  # Calculate the predicted site frequency spectrum (predSFS) for each model
  predSFS <- data %>%
    filter(model != "empirical") %>%
    group_by(row, col, model) %>%
    summarize(predicted_sfs = sum(count_mean))
  
  # Merge the predicted and observed site frequency spectra
  merged_data <- data[model=="empirical"] %>%
    left_join(predSFS, by = c("col", "row"))
  
  # Calculate the standardized residuals
  merged_data <- merged_data %>%
    mutate(residuals = ifelse(model.x == "empirical", 
                              (count_mean - predicted_sfs)/sqrt(count_mean), NA_real_))
  
  # Return the data with the standardized residuals
  return(merged_data)
}

# Resiudals
resid <- calculate_standardized_residuals(data = ld.dt1) 
hist(resid[model.y=="split_no_mig"]$residuals)
hist(resid[model.y=="split_mig"]$residuals)

res.model_Mig <- lm(residuals ~ 1, resid[model.y=="split_mig"])
mean(resid[model.y=="split_mig"]$residuals, na.rm = T)
confint(res.model_Mig, level=0.95)

res.model_noMig <- lm(residuals ~ 1, resid[model.y=="split_no_mig"])
mean(resid[model.y=="split_no_mig"]$residuals, na.rm = T)
confint(res.model_noMig, level=0.95)

bl <- colorRampPalette(c("navy", "royalblue", "lightskyblue"))(300)                      
re <- colorRampPalette(c("mistyrose", "red2", "darkred"))(400)

# Plotting residuals
resid.plot <- {
  
  resid[!count_mean==0] %>% 
    #mutate(residuals=datawizard::rescale(residuals, to=c(0,1))) %>% 
    ggplot(., aes(
      x=row,
      y=col,
      fill=residuals)) +
    geom_tile() +
    coord_flip() +
    facet_grid(~ model.y,
               labeller = labeller(model.y = c("empirical"="Empirical",
                                             "split_mig"="Split + migration",
                                             "split_no_mig"="Split"),
                                   vcf = c("Filtered"="Genome-wide"))) +
    theme_classic() +
    theme(axis.text.x=element_text(face="bold", size=16),
          axis.text.y=element_text(face="bold", size=16),
          text = element_text(face="bold", size=18), 
          legend.text = element_text(face="bold", size=16),
          strip.text=element_text(face="bold", size=16)) +
    scale_fill_gradientn(colours=c(bl, "white", re), na.value = "grey98",
                         limits=c(-41,55)) +
    scale_y_continuous(expand=c(0, 0)) +
    xlab(expression(paste(bold("Euro."), bolditalic(" D. pulex")))) +
    ylab(expression(paste(bold("NAm."), bolditalic(" D. pulex")))) +
    labs(fill="Residuals)")
  
}

# Make cowplot
library(cowplot)
fig <- plot_grid(SFS, resid.plot,
                 labels = c("A.", "B."), 
                 label_size = 20, nrow = 2)

# Output
ggsave("/scratch/csm6hg/figs/Fig3_sfs.pdf", fig, width = 12, height=8)
