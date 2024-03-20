# Gather optimal parameters across every run
# 2.14.22024
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)

# Working directory
setwd("/scratch/csm6hg/moments")

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
                        theta / (4 * mu * L_BUSCO), 
                        theta / (4 * mu * L_filt)),
         NA_pop_est = NA_pop_est / coeff,
         EU_pop_est = EU_pop_est / coeff,
         t_est = t_est / (2 * coeff),
         m_est = m_est * 2 * coeff))

# Output table of optimized values
write.table(opt_inf_op_raw, file = "optimizedRuns_moments.txt",sep = "\t", 
            append = F, quote = F, row.names = F)

# Convert output to values corresponding to "new" genome length estimates
opt_inf_output %>%
  mutate(new_coeff = ifelse(vcf == "BUSCO", 
                        L_BUSCO / L_BUSCO_new, 
                        L_filt / L_filt_new),
         NA_pop_est = NA_pop_est * new_coeff,
         EU_pop_est = EU_pop_est * new_coeff,
         t_est = t_est * new_coeff,
         m_est = m_est * new_coeff)

# Convert output without genome length inference by fixing NA_pop_est to 700,000
opt_inf_op_raw %>% 
  mutate(coeff = 700000 / NA_pop_est,
         NA_pop_est = NA_pop_est * coeff,
         EU_pop_est = EU_pop_est * coeff,
         t_est = t_est * 2 * coeff,
         m_est = m_est / (2 * coeff)) %>%
  select(vcf, model, NA_pop_est, EU_pop_est, t_est, m_est, theta, file) %>%
  arrange(vcf, model, file)

# Convert output using McCoy et al. 2014's ancestral population-based method.
# This is what we ultimately used in the paper.
opt_inf_op_raw %>% 
  mutate(N_anc_est = 200000 / EU_pop_est,
         NA_pop_est = NA_pop_est * N_anc_est,
         EU_pop_est = EU_pop_est * N_anc_est,
         t_est = t_est * 2 * N_anc_est,
         m_est = m_est / (2 * N_anc_est)) %>%
  select(vcf, model, NA_pop_est, EU_pop_est, t_est, m_est, theta, file) %>%
  arrange(vcf, model, file) %>% pull(t_est)
