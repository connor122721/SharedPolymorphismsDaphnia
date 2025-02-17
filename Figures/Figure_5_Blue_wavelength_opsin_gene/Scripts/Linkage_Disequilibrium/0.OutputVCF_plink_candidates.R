# Trans species polymorphism - haplotypes
# 9.13.2022
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(viridis)
library(doParallel)
library(readxl)
library(cowplot)
library(adegenet)
library(seqinr)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
out <- c("new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Read in total SNPs
tot_snps <- data.table(fread("metadata/snps_new"))

# Register cores
#doParallel::registerDoParallel(cores = 10)

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samps <- seqGetData(genofile, var.name = "sample.id")
fin <- fin[Sample %in% samps]

# Read in gene annotations
panth <- data.table(read_excel("../daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Gene analyses
pro <- read.fasta("/project/berglandlab/Karen/genomefiles/Daphnia.proteins.aed.0.6.fasta", seqtype="AA")
pro <- data.table(gene=getName(pro), AA.length=getLength(pro))
pro <- pro[,splice:=tstrsplit(gene, "-")[[2]]]

# Read in gtf file
gene.gtf <- data.table(fread("../daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))

# SNP metadata - TSPs
tot <- data.table(readRDS(file="/project/berglandlab/connor/data/classified_snps_filt.rds"))

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(fin[cont %in% as.character(compi[1])]$Sample)
samp2 <- as.vector(fin[cont %in% as.character(compi[2])]$Sample)

# Output VCF
write.table(data.table(fin[cont=="Daphnia.pulex.Europe"]$Sample, 
            fin[cont=="Daphnia.pulex.Europe"]$Sample), 
            file = "/project/berglandlab/connor/linkage/euro_pulex.txt",
            quote = F, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(data.table(fin[cont=="Daphnia.pulex.NorthAmerica"]$Sample, 
                       fin[cont=="Daphnia.pulex.NorthAmerica"]$Sample), 
            file = "/project/berglandlab/connor/linkage/nam_pulex.txt",
            quote = F, row.names = FALSE, col.names = FALSE, sep = "\t")
