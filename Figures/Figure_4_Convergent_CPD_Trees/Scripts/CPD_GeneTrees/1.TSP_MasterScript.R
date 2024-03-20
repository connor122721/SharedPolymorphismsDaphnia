# Trans species polymorphism
# 2.18.2023
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(DescTools)
library(readxl)
require(scales)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
out <- c("new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Read in total SNPs
tot_snps <- data.table(fread("metadata/snps_new"))

# Read in BUSCO SNPs
snps <- data.table(fread("metadata/busco_snps_new"))

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samps <- seqGetData(genofile, var.name = "sample.id")
fin <- fin[Sample %in% samps]

# Classification function for snps
snp_set <- function(x1, x2, thr) {
  if(x1==1 & x2==1) return("invariant")
  if(x1==0 & x2==0) return("invariant")
  if(x1==1 & x2==0) return("fixed")
  if(x1==0 & x2==1) return("fixed")
  if((x1>thr & x1<(1-thr)) & (x2>(1-thr) | x2<thr)) return("polymorphic_1")
  if((x2>thr & x2<(1-thr)) & (x1>(1-thr) | x1<thr)) return("polymorphic_2")
  if((x1>thr & x1<(1-thr)) & (x2>thr & x2<(1-thr))) return("shared_poly")
  return("filtered")
}

# Extract shared SNPs given snp set and get derived allele frequency
shared_snps <- function(x, s1, s2, compi, genofile) {
  #x=snps;s1=samp1;s2=samp2;compi=compi;genofile=genofile
  
  # Per Species 1 individual
  seqResetFilter(genofile)
  seqSetFilter(genofile,
               sample.id=c(s1),
               variant.sel=unique(x$variant.id))
  
  # Alelle counts for alternative allele
  snp.dt1 <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                        ac = seqAlleleCount(genofile, ref.allele=0L),
                        Species = compi[1],
                        num.samps = length(unique(s1))) %>%
    mutate(af = ac/(length(unique(s1))*2)) %>%
    mutate(maf = case_when(af > 0.5 ~ 1-af,
                           af <= 0.5 ~ af)) %>%
    mutate(af.bin=RoundTo(maf, 0.05, "floor"))
  
  # Per Species 2 individual
  seqResetFilter(genofile)
  seqSetFilter(genofile,
               sample.id=c(s2),
               variant.sel=unique(x$variant.id))
  
  # Alelle counts for alternative allele
  snp.dt2 <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                        ac = seqAlleleCount(genofile, ref.allele=0L),
                        Species = compi[2],
                        num.samps = length(unique(s2))) %>%
    mutate(af = ac/(length(unique(s2))*2)) %>%
    mutate(maf = case_when(af > 0.5 ~ 1-af,
                           af <= 0.5 ~ af)) %>%
    mutate(af.bin=RoundTo(maf, 0.05, "floor"))
  
  # Rbind both tables
  dt.return <- data.table(rbind(snp.dt1, snp.dt2))
  
}

# Include mainland european samples
ccc <- data.table(fread("/project/berglandlab/connor/candgene/samples_phase_var1", header=F))
colnames(ccc) <- c("Sample", "Continent")

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(ccc[Continent=="NorthAmerica"]$Sample)
samp2 <- as.vector(ccc[Continent=="Europe"]$Sample)

# Samples to make trees from
samps <- fin[Sample %in% c(samp1, samp2)]

### BUSCO SNPS ###
dt.busco <- shared_snps(x=snps, s1=samp1, s2=samp2, compi=compi, genofile=genofile)

# Apply classification function
dt1.busco <- data.table(dt.busco[,list(classified = snp_set(x1=af[Species==compi[1]], 
                                                            x2=af[Species==compi[2]], 
                                                            thr=0.01)),
                                 list(variant.id)] %>%
                          left_join(dt.busco) %>%
                          left_join(snps))
prop.table(table(dt1.busco$classified))

# Add snp totals - shared SNPs
fin.busco2 <- data.table(dt1.busco %>%
                           dplyr::select(-c(ac,af,n,af.bin,num.samps)) %>%
                           pivot_wider(names_from=Species, values_from=maf))

### Genome-wide SNPS ###
dt.fin <- shared_snps(x=tot_snps, s1=samp1, s2=samp2, compi=compi, genofile=genofile)

# Apply classification function
dt1.fin <- data.table(dt.fin[,list(classified = snp_set(x1=af[Species==compi[1]], 
                                                        x2=af[Species==compi[2]], 
                                                        thr=0.01)),
                             list(variant.id)] %>%
                        left_join(dt.fin) %>%
                        left_join(tot_snps))
prop.table(table(dt1.fin$classified))

# Add snp totals - shared SNPs
fin.tot2 <- data.table(dt1.fin %>%
                         dplyr::select(-c(ac,af,n,af.bin,num.samps)) %>%
                         pivot_wider(names_from=Species, values_from=maf))

# Rbind SNP table
fin <- data.table(rbind(data.table(fin.busco2, dt="BUSCO genes"), 
                        data.table(fin.tot2, dt="Genome-wide")))

# Simplify anotation
busco <- fin.busco2
busco[,simpleAnnot:=col]
busco[col%in%c("downstream_gene_variant", "upstream_gene_variant"), simpleAnnot:="Inter"]
busco[grepl("intergenic", col), simpleAnnot:="Inter"]
busco[grepl("stop", col), simpleAnnot:="Stop"]
busco[grepl("start", col), simpleAnnot:="Start"]
busco[grepl("splice", col), simpleAnnot:="Splice"]
busco[grepl("missense", col), simpleAnnot:="NS"]
busco[grepl("synonymous", col), simpleAnnot:="Syn"]
busco[grepl("initiator", col), simpleAnnot:="Start"]
busco[grepl("non_coding", col), simpleAnnot:="Inter"]
busco[grepl("intron", col), simpleAnnot:="Intron"]
busco[grepl("5_prime", col), simpleAnnot:="5'UTR"]
busco[grepl("3_prime", col), simpleAnnot:="3'UTR"]

tot <- fin.tot2
tot[,simpleAnnot:=col]
tot[col%in%c("downstream_gene_variant", "upstream_gene_variant"), simpleAnnot:="Inter"]
tot[grepl("intergenic", col), simpleAnnot:="Inter"]
tot[grepl("stop", col), simpleAnnot:="Stop"]
tot[grepl("start", col), simpleAnnot:="Start"]
tot[grepl("splice", col), simpleAnnot:="Splice"]
tot[grepl("missense", col), simpleAnnot:="NS"]
tot[grepl("synonymous", col), simpleAnnot:="Syn"]
tot[grepl("initiator", col), simpleAnnot:="Start"]
tot[grepl("non_coding", col), simpleAnnot:="Inter"]
tot[grepl("intron", col), simpleAnnot:="Intron"]
tot[grepl("5_prime", col), simpleAnnot:="5'UTR"]
tot[grepl("3_prime", col), simpleAnnot:="3'UTR"]

saveRDS(tot, "candgene/classified_snps_filt_exon.rds")
saveRDS(busco, "candgene/busco_classified_snps_filt_exon.rds")
#tot <- data.table(readRDS("candgene/classified_snps_filt_exon.rds"))
