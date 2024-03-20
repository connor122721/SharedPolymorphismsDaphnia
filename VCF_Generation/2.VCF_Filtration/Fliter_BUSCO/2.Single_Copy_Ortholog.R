# Compile single copy ortholog list
# 9.5.2022
# ijob -c 1 --mem=40G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(doParallel)
library(SNPRelate)
library(ggforce)
library(readxl)

# Working directory
setwd("/project/berglandlab/connor/metadata")
#doParallel::registerDoParallel(cores = 10)

# Read in GTF
gtf <- data.table(fread("Daphnia.aed.0.6.gtf") %>% 
                    mutate(transcript=tstrsplit(V9, ";")[[1]],
                           gene=tstrsplit(V9, ";")[[2]]) %>% 
                    mutate(transcript=str_remove_all(transcript, pattern = "transcript_id "),
                           gene=str_remove_all(gene, pattern = "gene_id ")) %>% 
                    select(!c(V6,V7,V8,V9)))
colnames(gtf)[1:5] <- c("chrom", "id", "ann", "start", "stop")

# Data
dt <- data.table(fread("busco_daphnia.tsv", fill = T))

# Restrict to SCO's
dt.tot <- data.table(dt[Status == "Complete"] %>% 
                       mutate(Gene = tstrsplit(Sequence, "-")[[1]],
                              Splice = tstrsplit(Sequence, "-")[[2]]))

# Restrict to highest score Splice
genes <- unique(dt.tot[!Splice=="RA"]$Gene)
dt.tot1 <- data.table(dt.tot %>%
                        group_by(Gene) %>% 
                        slice_max(order_by = Score, n=1))
colnames(dt.tot1)[1] <- c("Busco_id")

# Write table
#write.table(dt.tot1, file = "/project/berglandlab/connor/metadata/single_copy_orthologs", quote = F,row.names = F)

# Create SNP metadata
#out <-  c("../new_vcf2/daphnia.filt.qual.miss.rep.ann.gds")
out <- c("../new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")
#out <-  c("../new_vcf2/daphnia.whatshap.shapeit.ann.gds")
#out <-  c("../new_vcf2/daphnia.whatshap.ann.gds")

#seqVCF2GDS(vcf.fn = "../new_vcf2/daphnia.whatshap.ann.vcf.gz", out.fn = out, parallel = T)

# Open GDS
genofile <- seqOpen(out)

# Create gds metadata
snp.dt <- data.table(chrom=seqGetData(genofile, var.name = "chromosome"),
                     variant.id=seqGetData(genofile, var.name = "variant.id"),
                     position=seqGetData(genofile, var.name = "position"))

# Samples
samps <- seqGetData(genofile, var.name = "sample.id")

# Gets annotations
tmp <- seqGetData(genofile, "annotation/info/ANN") 
tmp.af <- data.table(variant.id= snp.dt$variant.id)
len1 <- tmp$length
len2 <- tmp$data

# Add overlapping annotations
snp.dt1 <- data.table(len=rep(len1, times=len1),
                      ann=len2,
                      id=rep(snp.dt$variant.id, times=len1))

# Extracting data between the 2nd and third | symbol
annot <- tstrsplit(snp.dt1$ann, "\\|")
snp.dt1 <- data.table(snp.dt1 %>% 
               mutate(class = annot[[2]],
                      gene = annot[[7]],
                      aa.change = annot[[10]],
                      protein = annot[[11]]))

# Collapsing additional annotations to original SNP vector length
snp.dt1.an <- snp.dt1[,list(n = unique(length(class), 
                                       length(col),
                                       length(gene), 
                                       length(aa.change), 
                                       length(protein)),
                            col = paste(class, collapse=","),
                            gene = paste(gene, collapse=","),
                            aa.change = paste(aa.change, collapse=","),
                            protein = paste(protein, collapse=",")),
                      list(variant.id=id)]

# Restrict to first gene in annotation column
snp.dt1.an <- data.table(snp.dt1.an %>% 
  mutate(col = tstrsplit(col,"\\,")[[1]],
         gene = tstrsplit(gene, "\\,")[[1]],
         aa.change = tstrsplit(aa.change, "\\,")[[1]],
         protein = tstrsplit(protein, "\\,")[[1]]))

# Extract gene name
unique(snp.dt1.an$gene)
length(unique(snp.dt1.an$gene))

# Read in master bed file
bed <- data.table(readRDS("/project/berglandlab/connor/data/filtered.dep.miss.rep.chr.strict.good.sites"))

# Restrict to unfiltered sites
snp.dt1.an.filt <- data.table(snp.dt1.an %>% 
  left_join(snp.dt, by=c("variant.id")) %>% 
  right_join(bed, by=c("chrom", "position", "variant.id")))

# Extract SNPS
genome <- data.table(snp.dt1.an.filt)

# Extract SNPS within BUSCO genes
busco <- data.table(snp.dt1.an.filt[gene %in% unique(dt$Sequence)])

# Write genome SNP table
write.table(genome, 
            file = "/project/berglandlab/connor/metadata/snps_new_strict", 
            quote = F, 
            sep = "\t",
            row.names = F)

# Write BUSCO SNP table
write.table(busco, 
            file = "/project/berglandlab/connor/metadata/busco_snps_new_strict", 
            quote = F, 
            sep = "\t",
            row.names = F)

#pdf("../figures/snp.numbers.busco.pdf", width = 10, height = 8)

# Pie chart of SNP categories
busco[col %in% c("synonymous_variant", "missense_variant")] %>% 
  group_by(chrom, col) %>% 
  summarize(value=n()) %>% 
  ggplot(., 
       aes(x=chrom, y=value, fill=col)) +
  geom_bar(stat="identity", width=1) +
  theme_minimal() + 
  labs(x="Chromosome",
       y="Number of SNPs",
       fill="Annotation") +
  theme(axis.text.x = element_text(face="bold", size=18, angle=-40),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))

#dev.off()

# Metadata
fin <- data.table(read.csv("../metadata/samples.fin.9.8.22.csv"))

# Run BUSCO SNP PCA
busco <- data.table(fread("../metadata/busco_snps_new"))
ccm_pca <- snpgdsPCA(genofile, 
                     autosome.only = F, 
                     sample.id = as.character(fin[!Species=="Daphnia obtusa"]$Sample),
                     num.thread = 1, 
                     snp.id = unique(busco$variant.id),
                     maf=0.01)

#saveRDS(ccm_pca, file = "pca.filt.busco.rds")

# Skree plot
plot(ccm_pca$varprop[1:20]*100)

# Merge by Sample
pca <- data.table(sample=ccm_pca$sample.id, PC1=ccm_pca$eigenvect[,1],
                  PC2=ccm_pca$eigenvect[,2], PC3=ccm_pca$eigenvect[,3],
                  PC4=ccm_pca$eigenvect[,4], PC5=ccm_pca$eigenvect[,5])

# Remove missing samples
pca <- data.table(pca[sample %in% fin$Sample])

# Merge PCA and metadata
pca <- data.table(merge(pca, fin,
                        by.x="sample", by.y="Sample"))

# Grouping PCA
fin.pca <- data.table(pca %>%
                        group_by(Species) %>%
                        summarise(count=n(),
                                  PC1.avg=mean(PC1),
                                  PC2.avg=mean(PC2),
                                  PC3.avg=mean(PC3)))

#pdf("../figures/pca12.filt.qual.miss.rep.dep.chr.ann.busco.nooutgroup.pdf", width = 10, height = 8)

# PC 1/2 Plot
pca %>%
  ggplot(., aes(x=PC1, 
                y=PC2, 
                fill=Continent,
                shape=Species)) +
  #facet_wrap(~Species+Continent) +
  geom_point(size=6, alpha=0.6) +
  geom_mark_ellipse() +
  theme_minimal() + 
  labs(x=paste("PC1 (", round(ccm_pca$varprop[[1]], digits=3)*100, " %)", sep=""),
       y=paste("PC2 (", round(ccm_pca$varprop[[2]], digits=3)*100, " %)", sep="")) +
  scale_fill_brewer(name = "Species complex", palette = "Set1") +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.text = element_text(size=16, face="bold.italic"),
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))

dev.off()

# Read in BUSCO SNPs
snps <- data.table(fread("busco_snps_new"))

# Restrict to BUSCO snps
seqResetFilter(genofile)
seqSetFilter(genofile, 
             variant.sel = unique(snps$variant.id))

# Create BUSCO SNP gds
seqGDS2VCF(genofile, vcf.fn = "../new_vcf2/daphnia.filtered.chr.busco.vcf")

# Read in annotations
panth <- data.table(read_excel("../../daphnia_ref/Daphnia_annotation_PANTHER.xls", ))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Restrict to BUSCO genes
pp <- data.table(panth[qseqid %in% snps$gene] %>% 
  select(gene=qseqid, bio_func, mol_func, cell_func))

# Write BUSCO SNP table
saveRDS(pp, file = "/project/berglandlab/connor/metadata/busco_panther")

# Read in 1% threshold TSPs
tsp <- data.table(read.table("../data/tsp_thresh_0.01", header = T))

# Restrict to single copy ortholog
sco <- panth[qseqid %in% c(tsp$gene)]

# Gene names
gene.names <- as.list(sco$qseqid) %>% unlist()

# GO list
go.list <- str_split(as.list(sco$GO), ";") 
names(go.list) <- gene.names

# Make annotations
annotation <- makeAnnotation(annotationList = go.list) 

# Hit list
mockHitlist = go.list %>% sapply(function(x){'GO:0031304' %in% x}) %>% 
  {go.list[.]} %>% 
  names

# Run over representation analysis
oraOut <- ora(annotation = annotation, pAdjust = "FDR")

# View output
out.dt <- data.table(oraOut$results)

### GENOME LENGTH FOR DAVID ###
gene.gtf <- data.table(fread("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf"))
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

library(GenomicRanges)
library(GenomicFeatures)
library(patchwork)

# Filter bed file 
bed <- fread("/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed", header = T)
colnames(bed) <- c("seqnames", "start", "end")

26610786/81967509
1603681/2552329

# Filter bed
bed.i <- makeGRangesFromDataFrame(bed, keep.extra.columns = T)

# Annotating CNVs w/metadata
dt.genome <- makeGRangesFromDataFrame(gene.gtf[sec=="transcript"] %>% 
                                      dplyr::select(chrom, start, stop), 
                                    keep.extra.columns = T)

# Filter out sites
genome <- data.frame(setdiff(x=dt.genome, y=bed.i, ignore.strand=T))
sum(genome$width)

# BUSCO genes
dt.busco <- makeGRangesFromDataFrame(gene.gtf[gene %in% dt.tot1$Sequence][sec=="transcript"]%>% 
                              dplyr::select(chrom, start, stop), 
                            keep.extra.columns = T)

# Filter out sites
busco <- data.frame(setdiff(x=dt.busco, y=bed.i, ignore.strand=T))
sum(busco$width)

# Genome-wide genes
gene.gtf[sec=="transcript"] %>% 
  mutate(length=stop-start) %>% 
  summarize(len=sum(length))
