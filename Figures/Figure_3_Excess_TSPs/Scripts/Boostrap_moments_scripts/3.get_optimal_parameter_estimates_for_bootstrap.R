# Gather optimal parameters across every run
# 2.14.22024
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)

# Working directory
setwd("/scratch/csm6hg/phd/moments")

# Genetic parameters for VCF filtered for SNPs 
mu <- 5.69e-9
L_filt <- 12e7
L_filt_new <- 26610786

# Load demographic inference results
files <- list.files(path="output", pattern = ".txt", full.names = T)
filt <- data.table(rbindlist(lapply(files, read.table), use.names = T, idcol = T, fill = T) %>% 
  mutate(vcf = "filt"))

# rename columns
colnames(filt) <- c("file", "model", "size", 
                    "theta", "ll", "NA_pop_est",
                    "EU_pop_est", "t_est", 
                    "m_est", "vcf")

# Get optimal estimate values
inf_output <- filt
opt_inf_output <- inf_output %>%
  filter(size == 20) %>%
  group_by(model, vcf, file) %>%
  filter(ll == max(ll)) %>%
  select(vcf, model, NA_pop_est, EU_pop_est, t_est, m_est, theta, file) %>%
  arrange(vcf, model, file)

# Get optimal estimate values, converted back into raw moments' units
opt_inf_op_raw <- data.table(opt_inf_output %>%
  mutate(coeff = ifelse(vcf == "BUSCO", 
                        as.numeric(theta) / (4 * mu * L_BUSCO), 
                        as.numeric(theta) / (4 * mu * L_filt)),
         NA_pop_est = NA_pop_est / coeff,
         EU_pop_est = EU_pop_est / coeff,
         t_est = t_est / (2 * coeff),
         m_est = m_est * 2 * coeff))

# Avergaes for testing
aPop=2e5
opt_inf_op_raw %>% 
  group_by(model) %>% 
  summarize(na_m=mean(NA_pop_est)*aPop,
            eu_m=mean(EU_pop_est)*aPop,
            t_m=mean(t_est)*2*aPop,
            m_m=mean(m_est)/(2*aPop))

# Adjust for testing different historic NEs
scaling=0.8
opt_inf_op_raw_rescale1 <- data.table(opt_inf_output %>%
  mutate(coeff = ifelse(vcf == "BUSCO", 
                        as.numeric(theta) / (4 * mu * L_BUSCO), 
                        as.numeric(theta) / (4 * mu * L_filt)),
         NA_pop_est = (NA_pop_est / coeff)*scaling,
         EU_pop_est = (EU_pop_est / coeff)*scaling,
         t_est = ((t_est / (2 * coeff))),
         m_est = ((m_est * 2 * coeff))))

# Adjust for testing different historic NEs
scaling=0.9
opt_inf_op_raw_rescale2 <- data.table(opt_inf_output %>%
  mutate(coeff = ifelse(vcf == "BUSCO", 
                        as.numeric(theta) / (4 * mu * L_BUSCO), 
                        as.numeric(theta) / (4 * mu * L_filt)),
         NA_pop_est = (NA_pop_est / coeff)*scaling,
         EU_pop_est = (EU_pop_est / coeff)*scaling,
         t_est = ((t_est / (2 * coeff))),
         m_est = ((m_est * 2 * coeff))))

# Output table of optimized values
write.table(opt_inf_op_raw, 
            file = "optimizedRuns_moments.txt",sep = "\t", 
            append = F, quote = F, row.names = F)

# Output table of optimized values - rescaled
write.table(opt_inf_op_raw_rescale1, 
            file = "optimizedRuns_moments_rescaleNe80.txt",sep = "\t", 
            append = F, quote = F, row.names = F)

# Output table of optimized values - rescaled
write.table(opt_inf_op_raw_rescale2, 
            file = "optimizedRuns_moments_rescaleNe90.txt",sep = "\t", 
            append = F, quote = F, row.names = F)

# Convert output using McCoy et al. 2014's ancestral population-based method.
# This is what we ultimately used in the paper.
opt_inf_op_raw %>% 
  mutate(N_anc_est = 200000 / (EU_pop_est),
         NA_pop_est = (NA_pop_est) * N_anc_est,
         EU_pop_est = (EU_pop_est) * N_anc_est,
         t_est = t_est * 2 * N_anc_est,
         m_est = m_est / (2 * N_anc_est)) %>%
  select(vcf, model, NA_pop_est, EU_pop_est, t_est, m_est, theta, file) %>%
  arrange(vcf, model, file) 
