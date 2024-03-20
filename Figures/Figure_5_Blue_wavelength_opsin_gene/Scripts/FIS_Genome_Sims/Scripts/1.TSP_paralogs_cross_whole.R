# Trans species polymorphism - FIS for gene
# 6.9.2023
# ijob -c 10 --mem=100G -p largemem -A berglandlab_standard
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)

# Executable in command line
arg <- commandArgs(TRUE)
SC1 <- arg[1] # Parent 1
SC2 <- arg[2] # Parent 2
focalSC = arg[3] # Focal group
SC1_ind <- arg[4] # Parent 1 ind
SC2_ind <- arg[5] # Parent 2 ind
iter <- arg[6] # Number file

#SC1="A"; SC2="C"; focalSC="1"; iter=10

# Working directory
setwd("/project/berglandlab/connor/paralogs/")

# Register cores
doParallel::registerDoParallel(cores = 15)

# Metadata
fin <- data.table(read.csv("../metadata/samples.9.8.22.csv"))

# Executable in command line
out <- c("../new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samps <- seqGetData(genofile, var.name = "sample.id")
fin <- fin[Sample %in% samps]

# Read in total SNPs
tot_snps <- data.table(fread("../metadata/snps_new"))
#tot_snps = tot_snps[gene %in% c("Daphnia11806-RA", "Daphnia00056-RA")]

# Simple classification function for snps
snp_hets <- function(x1, x2) {
  if( x1 == 0 & x2 == 0 ) return("invariant")
  if( x1 == 1 & x2 == 1 ) return("invariant")
  if( x1 == 1 & x2 == 0 ) return("fixed")
  if( x1 == 0 & x2 == 1 ) return("fixed")
  if( x1 == x2 ) return("het")
  return("other")
}

# SNP metadata
fin.tot2 <- data.table(readRDS(file="../data/classified_snps_filt.rds"))
#fin.tot2 = fin.tot2[gene %in% c("Daphnia11806-RA", "Daphnia00056-RA")]

# Use classsification script
classify_het_parents <- function(classi) {
  # classi="shared_poly"
  
  # Parent 1
  SC_1 <- unique(fin[SC %in% SC1][WildSequenced==0][Nonindependent==0]$Sample)
  
  # Parent 1
  SC_2 <- unique(fin[SC %in% SC2][WildSequenced==0][Nonindependent==0]$Sample)
  
  # Which snp set are we working on
  fin.tot_var <- fin.tot2[classified %in% classi]
  
  # See if all TSPs are hets in Parent 1
  scA_out1 <- foreach(l=1:length(SC_1), .combine = "rbind") %dopar% {
    
    # Reset filters
    seqResetFilter(genofile)
    seqSetFilter(genofile, 
                 sample.id=SC_1[l], 
                 variant.sel=fin.tot_var$variant.id)
    
    # Alelle counts for alternative allele
    snp.dtl <- data.table(data.table(variant.id=seqGetData(genofile, "variant.id"),
                                     ac = seqAlleleCount(genofile, ref.allele = 0L),
                                     Sample = SC_1[l],
                                     num = seqGetData(genofile, "$num_allele"),
                                     SC = SC1) %>% 
                            mutate(af = ac/2))

    # Finish
    return(snp.dtl)
  }
  
  # See if all TSPs are hets in Parent 2
  scC_out1 <- foreach(l=1:length(SC_2), .combine = "rbind") %dopar% {
    
    # Reset filters
    seqResetFilter(genofile)
    seqSetFilter(genofile, 
                 sample.id=SC_2[l], 
                 variant.sel=fin.tot_var$variant.id)
    
    # Alelle counts for alternative allele
    snp.dtl <- data.table(data.table(variant.id=seqGetData(genofile, "variant.id"),
                                     ac = seqAlleleCount(genofile, ref.allele=0L),
                                     Sample = SC_2[l],
                                     num = seqGetData(genofile, "$num_allele"),
                                     SC = SC2) %>% 
                            mutate(af = ac/2))
    
    # Finish
    return(snp.dtl)
    
  }
  
  # Rbind Parents 1 and 2
  sc1 <- data.table(rbind(scA_out1, scC_out1))
  
  # Filter biallelic sites
  out1 <- sc1[num==2][Sample %in% c(SC_1, SC_2)]
  
  # Find het variants across most SCs
  out2 <- data.table(out1 %>% 
    group_by(SC, af, variant.id) %>% 
    summarize(n=length(variant.id)) %>% 
    mutate(full=case_when(SC=="A" & n>=length(SC_1)*(2/3) ~ "Y",
                          SC=="C" & n>=length(SC_2)*(2/3) ~ "Y",
                          TRUE~"N")))
  
  # Filter out sites
  out2 <- out2[full=="Y"] %>% select(-c(n))
  out2 <- data.table(out2 %>% group_by(variant.id) %>% mutate(n.variant=length(variant.id)))
  
  # Call SNPs
  out3 <- out2[n.variant==2,list(classified=snp_hets(x1=af[SC==SC1],x2=af[SC==SC2])),list(variant.id)]
  table(out3$classified)
  
  ### TOTAL AF ###
  
  # Reset filters
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=SC_1, 
               variant.sel=fin.tot_var$variant.id)
  
  # Alelle counts for alternative allele
  snp.a <- data.table(data.table(variant.id=seqGetData(genofile, "variant.id"),
                                 ac = seqAlleleCount(genofile, ref.allele=0L),
                                 num = seqGetData(genofile, "$num_allele"),
                                 SC = SC1,
                                 parent = "parent") %>% 
                        mutate(af = ac/(length(SC_1)*2)))
  
  # Reset filters
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=SC_2, 
               variant.sel=fin.tot_var$variant.id)
  
  # Alelle counts for alternative allele
  snp.c <- data.table(data.table(variant.id=seqGetData(genofile, "variant.id"),
                                 ac = seqAlleleCount(genofile, ref.allele=0L),
                                 num = seqGetData(genofile, "$num_allele"),
                                 SC = SC2,
                                 parent = "parent") %>% 
                        mutate(af = ac/(length(SC_2)*2)))
  
  # Focal list of samples
  focalList <- c(fin[SC %in% focalSC | cont %in% focalSC | AxCF1Hybrid %in% focalSC]$Sample)
  
  # Reset filters
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=focalList, 
               variant.sel=fin.tot_var$variant.id)
  
  # Alelle counts for alternative allele
  snp.f <- data.table(data.table(variant.id=seqGetData(genofile, "variant.id"),
                                 ac = seqAlleleCount(genofile, ref.allele=0L),
                                 num = seqGetData(genofile, "$num_allele"),
                                 SC = "F1",
                                 parent = "F1") %>% 
                        mutate(af = ac/(length(focalList)*2)))
  
  # rbind
  snp.f1 <- data.table(rbind(snp.a, snp.c, snp.f))
  
  # Finish
  return(out3)
}

# Classification for both TSP and Genome-wide
out.shared <- classify_het_parents(classi=c("shared_poly"))
out.genome <- classify_het_parents(classi=c(unique(fin.tot2[classified %like% "polymorphic"]$classified)))

# Summarize
out.shared1 <- data.table(out.shared, snpset="TSP")
out.genome1 <- data.table(out.genome, snpset="Not TSP")
out1 <- data.table(rbind(out.shared1, out.genome1))
table(out1$classified, out1$snpset)

# Output
# saveRDS(out1, file = "parents_A_C_classified_shared.rds")
# out1 <- data.table(readRDS("parents_A_C_classified_shared.rds"))

# Focal list of samples
focalList <- c(fin[SC %in% focalSC | cont %in% focalSC | AxCF1Hybrid %in% focalSC]$Sample)

# See if all TSPs are hets in Focal sample
scf1_out1 <- foreach(l=1:length(focalList), .combine = "rbind") %dopar% {

  # Go through each SNPset
  ggg <- foreach(v=1:2, .combine = "rbind") %do% {
  # v=1; l=12
    
  # Reset filter to biallelic heterozygote TSPs
  print(paste("F1 sample:", l, v, sep=" "))
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=focalList[l], 
               variant.sel=out1[classified=="het"][snpset==unique(out1$snpset)[v]]$variant.id)
  
  # Get depth information
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  # Alelle frequency at heterozygous sites
  snp.dtf <- data.table(variant.id = seqGetData(genofile, "variant.id"),
                        length = seqGetData(genofile, "$num_allele"),
                        all = seqGetData(genofile, "$dosage")[1,],
                        all_num = seqAlleleCount(genofile, ref.allele = 0L),
                        dos.alt = data.table(data.table(t(ad$data)) %>%
                                               filter(row_number() %% 2 == 1))$V1,
                        dos.ref = as.numeric(data.table(t(dp))$V1),
                        Sample = focalList[l],
                        iteration = l,
                        snpset = unique(out1$snpset)[v]) %>% 
    mutate(af.dos = dos.alt/dos.ref)
  
  return(snp.dtf)
  }
  
  return(ggg)
}

# Merge HWE with other dosage data
fin.dt <- data.table(scf1_out1, focalSC=1, SC1=SC1, SC2=SC2)

# Add metadata
fin.dt1 <- data.table(fin.dt %>%
                left_join(tot_snps, by=c("variant.id")))

# Distribution of genotype frequencies in F1s
f1_daphnia <- {
  
  fin.dt1[all %in% c(0,1,2)] %>% 
    group_by(all, variant.id, snpset, col, gene) %>% 
    summarize(mean=n()) %>% 
    ggplot(., 
           aes(x=as.factor(all), 
               y=(mean/48)*100, 
               color=snpset, 
               group=variant.id)) +
      geom_line() +
      geom_point() + 
      facet_wrap(gene~snpset) +
      ylim(c(0,100)) +
      theme_bw() +
      labs(x="Allele dosage",
           y="F1 genotype frequency (%)") +
      theme(legend.position = "none", 
            strip.text = element_text(face="bold", size=18),
            legend.text = element_text(face="bold", size=18),
            axis.text.x = element_text(face="bold", size=18),
            axis.text.y = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=20),
            axis.title.y = element_text(face="bold", size=20),
            axis.title = element_text(face="bold", size=20))

}

ggsave("Daphnia11806-RA_Daphnia00056-RA_F1genotype.pdf", f1_daphnia, width = 10, height=8)

# Output w/ metadata
saveRDS(fin.dt1, file = "AxC_tot_het_withmeta_Daphnia11806-RA_Daphnia00056-RA.rds")
