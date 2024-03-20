# PopGenome - Diversity & Divergence
# 11.16.2022
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(PopGenome)
library(tidyverse)
library(viridis)
library(doParallel)
library(readxl)
library(cowplot)
require(scales)

# Working directory
setwd("/project/berglandlab/connor/")

# Register cores
#doParallel::registerDoParallel(cores = 10)

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
# out <- c("new_vcf2/daphnia.filt.mlg.genome.11.18.22.vcf.gz")
out <- c("new_vcf2/daphnia.filtered.chr.busco.vcf.gz")


# Read in total SNPs BUSCO Genes
snps <- data.table(fread("metadata/busco_snps_new"))

# Read in total SNPs
tot <- data.table(readRDS("data/classified_snps_filt.rds"))

# Chromosome metadata
chrom.meta <- data.table(rbind(fread("metadata/D84Agoodscaffstouse.12chromosomes.bed")))
names(chrom.meta) <- c("chrom", "start", "stop")

# Bed files
bed <- data.table(rbind(fread("data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed")))
names(bed) <- c("chrom", "start", "stop")

# Gets denominator for statistic
bed <- bed[chrom %in% unique(chrom.meta$chrom)]
bed <- bed[,len:=stop-start]
bed.tot <- bed[, list(total=sum(len)), list(chrom)]

# Names of species
continent <- c(unique(fin$cont))

# Reads in vcf with bootstrap sampling
vcf.fun <- function(lower.maf, upper.maf) {
 
    # Go through each chromosome
    pp <- foreach(chr.num = 1:dim(chrom.meta)[1], .combine="rbind", .errorhandling = "remove") %do% {
      # lower.maf=0; upper.maf=0.5; chr.num=1
      
      message(paste("Chromosome:", chr.num, sep=" "))
      
      # Load tabindex vcf file - load in by chromosome
      vcf <- readVCF(filename = out, 
                     tid = unique(chrom.meta$chrom)[chr.num], 
                     gffpath = "new_vcf2/Daphnia.aed.0.6.gff",
                     samplenames = unique(fin$Sample),
                     frompos = 1, 
                     topos = chrom.meta[chr.num]$stop, 
                     numcols = 10000, 
                     include.unknown = TRUE,
                     approx = FALSE)
      
      # Filter out sites based on MAF
      genofile <- set.filter(vcf, missing.freqs = FALSE, minor.freqs = TRUE, 
                             maf.lower.bound = lower.maf, maf.upper.bound = upper.maf)
      
      genofile <- set.populations(genofile, 
                  list("pulex_euro"=c(fin[cont=="Daphnia.pulex.Europe"]$Sample),
                       "pulex_nam"=c(fin[cont=="Daphnia.pulex.NorthAmerica"]$Sample),
                       "pulic_euro"=c(fin[cont=="Daphnia.pulicaria.Europe"]$Sample),
                       "pulic_nam"=c(fin[cont=="Daphnia.pulicaria.NorthAmerica"]$Sample),
                       "ob_nam"=c(fin[cont=="Daphnia.obtusa.NorthAmerica"]$Sample)))
      
      # Outgroup
      genome_syn <- set.outgroup(genofile,
                     c(fin[cont=="Daphnia.obtusa.Europe"]$Sample))
      
      # Diversity and neutrality statistics
      genome_syn <- set.synnonsyn(genome_syn, ref.chr="totalHiCwithallbestgapclosed.fa")
      genome_syn <- diversity.stats.between(genome_syn, subsites="syn")
      genome_syn <- diversity.stats(genome_syn, subsites="syn")
      
      # QC
      # genome_syn@region.data@outgroup
      
      # Extract data
      dxy_with <- data.table(genome_syn@nuc.diversity.within/genome_syn@n.biallelic.sites)
      dxy_bet <- data.table(genome_syn@nuc.diversity.between/genome_syn@n.biallelic.sites)
      
      # Denominator for pi
      denom.pi <- (bed.tot[chrom==chrom[chr.num]]$total - 
                     sum(table(genome_syn@region.stats@minor.allele.freqs[[1]][1,]>lower.maf)[2]))
      
      # Some replicates will not have any sites > MAF 0.25
      if(is.na(denom.pi) =="TRUE") {
        denom.pi = genome_syn@n.sites
      }
      
      # Compile data
      total <- data.table(dxy.nampul = dxy_with[[2]],
                          dxy.eurpul = dxy_with[[1]],
                          dxy.betpul = dxy_bet[[1]],
                          chrom = chrom.meta[chr.num]$chrom,
                          n.tot.Sites = genome_syn@n.sites,
                          denom.pi = denom.pi,
                          nMAF.Above = sum(table(genome_syn@region.stats@minor.allele.freqs[[1]][1,]>upper.maf)[2]), 
                          nMAF.Below = sum(table(genome_syn@region.stats@minor.allele.freqs[[1]][1,]<=upper.maf)[2]),
                          bedSites = bed.tot[chrom==chrom.meta[chr.num]$chrom]$total,
                          biall.Sites = genome_syn@n.biallelic.sites, 
                          Sites.inc = table(genome_syn@region.data@included)[2],
                          n.samp = length(unique(fin$Sample)),
                          up.maf = upper.maf,
                          low.maf = lower.maf)
      
      # Quality control - if replicate does not have any sites > MAF upper
      if(is.na(total$nMAF.Below) == "TRUE") {
        total$nMAF.Below = genome_syn@n.biallelic.sites
        total$Sites.inc = table(genome_syn@region.data@included)[1]
      }
      
      # Finish by chrom
      return(total)
    }
}

# dxy function
dt <- vcf.fun(lower.maf=0.01, upper.maf = 0.45)

#saveRDS(pp, file="data/dxy_synonymous_0.01_0.45_maf_genome.rds")

# dxy data for synonymous sites
dt <- data.table(readRDS("data/dxy_synonymous_0.01_0.45_maf_genome.rds"))

# Neutral expectation of number of TSPs (upper bound)
dAB=mean(dt$dxy.betpul) # Between species dxy
dA=mean(dt$dxy.nampul) # within NAm. D. pulex dxy
dB=mean(dt$dxy.eurpul) # within Euro D. pulex dxy
num.snps=sum(dt$biall.Sites) # Total number of BUSCO Gene SNPs

e1=exp(-((dAB-max(dA,dB))/dA))
e2=exp(-((dAB-max(dA,dB))/dB))
prop.tsp=e1*e2

# Expected Number of TSPs
(prop.tsp*num.snps)

# Actual Number of TSPs - Synonymous above freq.
emp=length(unique(tot[col=="synonymous_variant"][classified=="shared_poly"]$variant.id))

# Enrichment
log2(emp/(prop.tsp*num.snps))
