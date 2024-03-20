# Shared polymorphism for VCF files
# 9.11.2022
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(viridis)
library(doParallel)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
out <-  c("new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samps <- seqGetData(genofile, var.name = "sample.id")
fin <- fin[Sample %in% samps]

# Comparisons
dt <- data.table(expand.grid(unique(fin$cont), unique(fin$cont)))

# Remove pairwise duplicates
compare <- data.table(dt[!duplicated(t(apply(dt, 1, sort))),])
rownames(compare) <- 1:dim(compare)[1]
colnames(compare) <- c("Species1", "Species2")

# Reduce number of comparisons
compare <- compare[!Species1 %like% c("obtusa")][!Species2 %like% c("obtusa")] 
#compare <- compare[!Species1==Species2]

# Read in BUSCO SNPs
snps <- data.table(fread("metadata/busco_snps_new"))

# Read in total SNPs
#snps <- data.table(fread("metadata/snps_new"))

# Restrict to Genome-wide snps
seqSetFilter(genofile, variant.sel = unique(snps$variant.id))

# Register cores
doParallel::registerDoParallel(cores = 15)

# Classification function for snps
snp_set <- function(x) {
  if(x[1] == 1 & x[2] == 1) return("invariant")
  if(x[1] == 0 & x[2] == 0) return("invariant")
  if(x[1] == 1 & x[2] == 0) return("fixed")
  if(x[1] == 0 & x[2] == 1) return("fixed")
  if((x[1]>0 & x[1]<1) & (x[2]==1 | x[2]==0)) return("polymorphic_1")
  if((x[2]>0 & x[2]<1) & (x[1]==1 | x[1]==0)) return("polymorphic_2")
  if((x[1]>0 & x[1]<1) & (x[2]>0 & x[2]<1)) return("shared_poly")
}

# Run proportion snps shared across comparison groups -2DSFS
sfs2d.boot <- function(i, num.sampi, booti) { 
  #i=9; num.sampi=4; booti=1
  
  # Progress message
  message(paste("Species comp:", i, compare[i]$Species1, compare[i]$Species2, sep=" "))
  
  # Reset any previous filters
  seqResetFilter(genofile)
  message("Getting allele count information")
  
  # Species comparisons
  compi <- c(as.character(compare[i]$Species1), 
             as.character(compare[i]$Species2))
  
  # Samples
  samp1 <- as.vector(sample(fin[cont %in% 
                                  as.character(compare[i]$Species1)]$Sample,
                            num.sampi))
  
  samp2 <- as.vector(sample(fin[cont %in% 
                                  as.character(compare[i]$Species2)]$Sample, 
                            num.sampi))
  
  # Per Species 1 individual
  seqSetFilter(genofile, 
               sample.id=c(samp1),
               variant.sel=unique(snps$variant.id))
  
  # Alelle counts for alternative allele
  snp.dt1 <- data.table(variant.id=seqGetData(genofile, "variant.id"), 
                        ac = seqAlleleCount(genofile, ref.allele=0L),
                        Species = compi[1],
                        num.samples = num.sampi) %>% 
    mutate(af = ac/(num.samples*2)) %>% 
    mutate(maf = case_when(af > 0.5 ~ 1-af,
                           af <= 0.5 ~ af)) %>% 
    mutate(af.bin=cut_interval(af, n = num.sampi*2)) %>% 
    mutate(af.bin=case_when(af.bin=="[0,0.125]" & af==0 ~ "[0]",
                            af.bin=="[0,0.135]" & !af==0 ~ "(0,0.125]",
                            af.bin=="(0.875,1]" & af==1 ~ "[1]",
                            af.bin=="(0.875,1]" & !af==1 ~ "(0.875,1]",
                            TRUE ~ as.character(af.bin)))
  
  # Per Species 2 individual
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=c(samp2),
               variant.sel=unique(snps$variant.id))
  
  # Alelle counts for alternative allele
  snp.dt2 <- data.table(variant.id=seqGetData(genofile, "variant.id"), 
                        ac = seqAlleleCount(genofile, ref.allele=0L),
                        Species = compi[2],
                        num.samples = num.sampi) %>% 
    mutate(af = ac/(num.samples*2)) %>% 
    mutate(maf = case_when(af > 0.5 ~ 1-af,
                           af <= 0.5 ~ af)) %>% 
    mutate(af.bin=cut_interval(af, n = num.sampi*2)) %>% 
    mutate(af.bin=case_when(af.bin=="[0,0.125]" & af==0 ~ "[0]",
                            af.bin=="[0,0.135]" & !af==0 ~ "(0,0.125]",
                            af.bin=="(0.875,1]" & af==1 ~ "[1]",
                            af.bin=="(0.875,1]" & !af==1 ~ "(0.875,1]",
                            TRUE ~ as.character(af.bin)))
  
  # Rbind both tables
  dt <- rbind(snp.dt1, snp.dt2)
  
  # Apply classification function
  dt1 <- dt[,list(classified = snp_set(af)),
            list(variant.id)]
  
  # Shared alleles of interest
  seqSetFilter(genofile, variant.id = c(dt1[classified == "fixed"]$variant.id))
  
  # Go through every individual  
  sdtt1 <- foreach(indie=1:length(fin$Sample), .combine = "rbind", .errorhandling = "remove") %dopar% { 
    
      # Progress message
      print(paste("Comparison:", indie, sep=" "))
      
      # Reset filter
      seqResetFilter(genofile, sample = T, variant = F)
      
      # Iteratively go through each individual within every pop
      seqSetFilter(genofile, 
                   sample.id = c(fin[indie]$Sample),
                   variant.sel = c(dt1[classified == "fixed"]$variant.id))
      
      # Extract dosage information
      dos <- seqGetData(genofile, "$dosage")
      
      # Summarize
      indie_dt <- data.table(t(table(dos)), 
                             Sample = fin[indie]$Sample,
                             Species = fin[indie]$Species,
                             Continent = fin[indie]$Continent,
                             cont = fin[indie]$cont,
                             Spp1_comp = compi[1],
                             Spp2_comp = compi[2],
                             total_fixed_snps = length(unique(dt1[classified == "fixed"]$variant.id)),
                             total_shared_snps = length(unique(dt1[classified == "shared_poly"]$variant.id)),
                             total_invariant_snps = length(unique(dt1[classified == "invariant"]$variant.id)),
                             total_poly1_snps = length(unique(dt1[classified == "polymorphic_1"]$variant.id)),
                             total_poly2_snps = length(unique(dt1[classified == "polymorphic_2"]$variant.id)),
                             num.samples = num.sampi,
                             id = indie)
      
      # Finish individual
      return(indie_dt)
      
    }
  
  # Summarize output
  summ <- data.table(sdtt1 %>%
                       group_by(Dosage=dos, Species, Continent, Spp1_comp, Spp2_comp) %>% 
                       summarize(N = mean(N, na.rm = T),
                                 mean.total_fixed_snps = mean(total_fixed_snps, na.rm = T),
                                 mean.total_shared_snps = mean(total_shared_snps, na.rm = T),
                                 mean.total_invariant_snps = mean(total_invariant_snps, na.rm = T),
                                 mean.total_poly1_snps = mean(total_poly1_snps, na.rm = T),
                                 mean.total_poly2_snps = mean(total_poly2_snps, na.rm = T),
                                 num.samples = unique(num.samples)))
  
  # Write output
  saveRDS(summ, 
            file=paste("/project/berglandlab/connor/fixed/data/fixed.100reps.4samps", 
                       compare[i]$Species1, 
                       compare[i]$Species2, 
                       boot,
                       i,
                       "rds", sep="."))
  
  # Return
  return(summ)
}

# Run bootstrap
o.fin <- foreach(boot = 1:100, .combine="rbind", .errorhandling = "remove") %do% {
  #boot=1
  
  # Progress message
  print(paste("Bootstrap:", boot, sep=" "))
  
  # Run comparison groups
  o <- foreach(k = 1:length(compare$Species1), .combine="rbind", .errorhandling = "remove") %do% { 
    out <- sfs2d.boot(i=k, num.sampi = 4, booti=boot)
    
    # Finish comprison
    return(out)
  }
  
  # Finish bootstrap
  return(o)
}

# Write output
saveRDS(o.fin, file="/project/berglandlab/connor/fixed/data/fixed.100reps.4samps.new.rds")
