# Trans species polymorphism - FIS for gene
# 3.12.2024
# ijob -c 10 --mem=50G -p largemem -A berglandlab_standard
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)

# Executables
SC1="C"
focalSC="1"
iter=0

# Working directory
setwd("/project/berglandlab/connor/chapter1/paralogs/")

# Register cores
doParallel::registerDoParallel(cores = 10)

# Metadata
fin <- data.table(read.csv("../../metadata/samples.9.8.22.csv"))
samps <- fread("../../../MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
samps[LabGenerated==1][SC=="selfedC"]  # CxC lab
samps[LabGenerated==1][AxCF1Hybrid==1] # AxC lab

# Load GDS
genofile <- seqOpen("../../../Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

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
fin.tot2 <- data.table(readRDS(file="/scratch/csm6hg/data/classified_snps_filt.rds"))
#fin.tot2 = fin.tot2[gene %in% c("Daphnia11806-RA", "Daphnia00056-RA")]
fin.tot2 = fin.tot2[classified %in% c("polymorphic_1","polymorphic_2","shared_poly")]

# Restrict to conservative SNPs
con <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
fin.tot2 <- fin.tot2[ch %in% unique(con$ch)]

# SNP Metadata
snp.dt <- data.table(variant.id = seqGetData(genofile, "variant.id"),
                     position = seqGetData(genofile, "position"),
                     chrom = seqGetData(genofile, "chromosome"))
snp.dt <- snp.dt[position %in% fin.tot2$position]
#snp.dt <- snp.dt[chrom=="Scaffold_1931_HRSCAF_2197"][position %in% 6350219:6351796]
#snp.dt <- snp.dt[position %in% fin.tot2[gene=="Daphnia11806-RA"]$position]

# Parent 1
SC_1 <- unique(samps[SC %in% SC1][WildSequenced==0][Nonindependent==0]$clone)
  
# See if all TSPs are hets in Parent 1
p1 <- foreach(l=1:length(SC_1), .combine = "rbind") %dopar% {
    
    # Reset filters: l=2
    seqResetFilter(genofile)
    seqSetFilter(genofile, 
                 sample.id=SC_1[l], 
                 variant.sel=snp.dt$variant.id)
    
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
  
# Filter biallelic sites
out1 <- p1[num==2]
  
# Find het variants across most SCs
out2 <- data.table(out1 %>% 
    group_by(SC, af, variant.id) %>% 
    summarize(n=length(variant.id)) %>% 
    mutate(full=case_when(n>=length(SC_1)*(2/3) ~ "Y",
                          TRUE~"N")))
  
# Extract all het. sites in C parents
out3 <- out2[full=="Y"][af==0.5] %>% 
    mutate(classified="het") %>% 
    select(variant.id,classified)

# Focal list of samples
focalList <- c(samps[LabGenerated==1][SC=="selfedC"]$clone)

# See if all TSPs are hets in Focal sample
f1 <- foreach(l=1:length(focalList), .combine = "rbind") %dopar% {

  # Reset filter to biallelic heterozygote TSPs: l=4
  print(paste("F1 sample:", l, sep=" "))
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=focalList[l], 
               variant.sel=out3$variant.id)
  
  # Get depth information
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  # Alelle frequency at heterozygous sites
  snp.dtf <- data.table(variant.id = seqGetData(genofile, "variant.id"),
                        position = seqGetData(genofile, "position"),
                        chrom = seqGetData(genofile, "chromosome"),
                        length = seqGetData(genofile, "$num_allele"),
                        all = seqGetData(genofile, "$dosage")[1,],
                        all_num = seqAlleleCount(genofile, ref.allele = 0L),
                        dos.alt = data.table(data.table(t(ad$data)) %>%
                                               filter(row_number() %% 2 == 1))$V1,
                        dos.ref = as.numeric(data.table(t(dp))$V1),
                        Sample = focalList[l],
                        iteration = l) %>% 
    mutate(af.dos = dos.alt/dos.ref)
  
return(snp.dtf)
}

# Add metadata
fin.dt1 <- data.table(f1 %>%
             left_join(fin.tot2, 
                       by=c("position", "chrom")))

# Plotting data
fin.dt2 <- data.table(fin.dt1[all %in% c(0,1,2)] %>% 
  mutate(classified=case_when(classified=="shared_poly"~"TSP",
                              TRUE ~ "Not TSP")) %>% 
  group_by(all, gene, classified) %>% 
  summarize(mean=n()) %>% 
  group_by(classified, gene) %>% 
  mutate(sumObs=sum(mean),
         prop=mean/sumObs*100))

# Gene-wise average
gene <- data.table(fin.dt2 %>% 
  group_by(all, classified) %>% 
  summarize(meanProp=mean(prop),
            medProp=median(prop)))

# Distribution of genotype frequencies in F1s
f1_daphnia <- {
  
  fin.dt2[!gene=="Daphnia11806-RA"] %>% 
    ggplot(., 
           aes(x=as.factor(all), 
               y=prop, 
               color=classified,
               group=gene)) +
      geom_line(alpha=0.05) +
      ylim(c(0,100)) +
      geom_line(data=gene, 
                aes(x=as.factor(all), y=meanProp,color=classified, group=1), size=2) +
      geom_line(data=fin.dt2[gene=="Daphnia11806-RA"], color="blue") +
      facet_wrap(~classified) +
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

ggsave("/scratch/csm6hg/figs/CxC_F1genotype_dist.pdf", f1_daphnia, width = 10, height=6)

# Output w/ metadata
saveRDS(fin.dt1, file = "/scratch/csm6hg/data/CxC_tot_het_withmet.rds")
