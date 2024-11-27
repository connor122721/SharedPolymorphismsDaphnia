# Code to generate site frequency spectra from moments and calculate shared polymorphisms
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Setup
library(tidyverse)
library(data.table)

# Working directory
setwd("/scratch/csm6hg/phd/moments/sfss/sfs_csvs_rescale80")

# Read in SFS function
load_sfs <- function(name, min=1, scale=1) {
  
  # Read in SFS from CSV file
  sfs <- fread(name, header=FALSE)
  print(name)
  
  # Mask absent allele corner
  sfs[1, 1] <- 0
  
  # Get size of SFS. This currently only supports square SFSs.
  size <- sfs %>% count %>% as.numeric
  
  # "Pivot" data frame from CSV into data frame that can be plotted 
  sfs_pivoted <- sfs %>%
    pivot_longer(1:size, names_to="row", values_to="count") %>% 
    mutate(col = rep(1:size, each=size))
  
  # Reorder factor levels so that the columns of the SFS get plotted in the correct
  # order.
  sfs_pivoted$row <- factor(as.factor(sfs_pivoted$row),
                            paste0("V", seq(size)))
  
  # Scale counts in SFS
  sfs_pivoted <- sfs_pivoted %>%
    mutate(count = count * scale)
  
  # Round values less min down to 0
  sfs_pivoted <- sfs_pivoted %>%
    mutate(count = ifelse(count <= min, 0, count))
  
  # Return
  sfs_pivoted
}

# List of SFSs
files <- list.files(full.names = TRUE, pattern = '*.csv$')
files <- files[!files %like% "_100."]

# Read in SFSs
ld.sfs <- lapply(files, load_sfs)
names(ld.sfs) <- files
ld.dt <- data.table(rbindlist(ld.sfs, use.names = T, idcol = T)) %>% 
  mutate(model=case_when(.id %like% "no_mig" ~ "split_no_mig",
                         .id %like% "empir" ~ "empirical",
                         .id %like% "from_es" ~ "from_ests",
                         TRUE ~ "split_mig"),
         size=case_when(.id %like% "_100" ~ 100,
                        .id %like% "_20" ~ 20),
         rep=tstrsplit(.id, ".", fixed=T)[[3]],
         vcf="Filtered")

# saveRDS(ld.dt, file = "/scratch/csm6hg/phd/data/ld.dt.tot.20x20.rescale80.rds")
# ld.dt <- data.table(readRDS("/scratch/csm6hg/phd/data/ld.dt.tot.20x20.rescale50.rds"))

# Average values for SFS
ld.dt1 <- data.table(ld.dt %>% 
              group_by(row, col, model, vcf, size) %>% 
              summarize(count_mean = mean(count, na.rm = T),
                        count_med = median(count, na.rm = T)))

# Denominator of all alleles
denom = sum(ld.dt1[model == "empirical"]$count_mean)

# SAP calculation - confidence intervals
sum(ld.dt1[model == "empirical"][!row=="V1"][!col==1]$count_mean)/denom

sum(ld.dt1[model == "split_mig"][!row=="V1"][!col==1]$count_mean)/denom

sum(ld.dt1[model == "split_no_mig"][!row=="V1"][!col==1]$count_mean)/denom

# Get SAP per replicate
ld.dt.sap <- data.table(ld.dt %>% 
                filter(row!="V1" & col !=1) %>% 
                group_by(model, rep, size) %>% 
                summarize(sap=(sum(count)/denom)*100))

# Calculate the confidence interval
sap.model_Mig <- lm(sap ~ 1, ld.dt.sap[model=="split_mig"][size==20][sap<3 & sap>0])
confint(sap.model_Mig, level=0.95)

sap.model_noMig <- lm(sap ~ 1, ld.dt.sap[model=="split_no_mig"][size==20][sap<3 & sap>=0])
confint(sap.model_noMig, level=0.95)

# Plotting SFS
SFS <- {ld.dt1[!model=="from_ests"] %>% 
    ggplot(., 
           aes(x=row, 
               y=col, 
               fill=log10(count_mean))) +
    geom_tile() + 
    facet_grid(vcf ~ model,
               labeller = labeller(sfs = c("empirical"="Empirical",
                                           "from_ests"="From parameter estimates",
                                           "split_mig"="Split-with-migration model",
                                           "split_no_mig"="Split-without-migration model"),
                                   vcf = c("Filtered"="Genome-wide"))) +
    theme_classic() +
    theme(axis.text.x=element_text(face="bold", size=16),
          axis.text.y=element_text(face="bold", size=16),
          text = element_text(face="bold", size=18), 
          legend.text = element_text(face="bold", size=16),
          strip.text=element_text(face="bold", size=16)) +
    scico::scale_fill_scico(palette = "vik", na.value="white", begin = 1, end = 0) +
    coord_fixed() +
    scale_y_continuous(expand=c(0, 0)) +
    xlab(expression(paste(bold("Euro."), bolditalic(" D. pulex")))) +
    ylab(expression(paste(bold("NAm."), bolditalic(" D. pulex")))) +
    labs(fill=bquote(bold(10^SNPs)))
  
}

# Merge list of summary stats
files.sum <- list.files(full.names = TRUE, pattern = '*.txt$', path = "../../output_eSMC")
files.sum <- files.sum[files.sum %like% "rescaled80"]
sap <- lapply(files.sum, function(x){fread(x,header=F,fill=T)})
names(sap) <- files.sum

# Merge
sap.dt <- data.table(rbindlist(sap, use.names = T, idcol = T)) %>% 
  mutate(rep=tstrsplit(.id, ".", fixed=T)[[7]],
         vcf="Filtered")

colnames(sap.dt) <- c("file", "filt", "ns", "model",
                      "AIC", "BIC", "sap", "sap_ests", "extra", "rep", "vcf")

# saveRDS(sap.dt, file = "/scratch/csm6hg/phd/data/sap.dt.tot.20x20.rescale80.rds")
# sap.dt <- data.table(readRDS("/scratch/csm6hg/data/sap.dt.tot.20x20.rds"))

# Average SAP and get CI
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
  
  resid[!count_mean==0][!model.x=="from_ests"][!model.y=="from_ests"] %>% 
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
