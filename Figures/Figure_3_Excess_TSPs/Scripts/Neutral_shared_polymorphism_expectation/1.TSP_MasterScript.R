# Trans species polymorphism
# 9.26.2022
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(viridis)
library(car)
library(DescTools)
library(readxl)
library(cowplot)
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
shared_snps <- function(x,s1,s2,compi,genofile) {
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

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(fin[cont %in% as.character(compi[1])]$Sample)
samp2 <- as.vector(fin[cont %in% as.character(compi[2])]$Sample)

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

#saveRDS(tot, "data/classified_snps_filt.rds")
#saveRDS(busco, "data/busco_classified_snps_filt.rds")
tot <- data.table(readRDS("data/classified_snps_filt.rds"))

# Restrict to conservative TSP set
con.tsp <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
tot <- tot[ch %in% con.tsp$ch]

### 2D SFS ###
test <- data.table(fin[ch %in% con.tsp$ch] %>%
            mutate(af.bin.x=RoundTo(Daphnia.pulex.NorthAmerica, 0.01, "floor"),
                   af.bin.y=RoundTo(Daphnia.pulex.Europe, 0.01, "floor")) %>% 
            group_by(af.bin.x, af.bin.y, col, dt) %>%
            summarize(nSNPs=length(unique(variant.id))) %>% 
            mutate(col=case_when(col=="missense_variant" ~ "Non-synonymous",
                                 col=="synonymous_variant" ~ "Synonymous",
                                 col %in% c("intergenic_region", 
                                            "downstream_gene_variant",
                                            "upstream_gene_variant") ~ "Intergenic",
                                 TRUE ~ col)) %>% 
            filter(col %in% c("Intergenic", "Synonymous", "Non-synonymous", 
                              "Intron_variant", "5_prime_UTR_variant",
                              "3_prime_UTR_variant")))

# saveRDS(test, "data/classified_snps_filt_transform.rds")

# 2Dsfs
sfs2d <- {test[dt=="Genome-wide"][!c(af.bin.x==0 & af.bin.y==0)] %>%
  ggplot(., 
         aes(x=af.bin.x,
             y=af.bin.y,
             fill=log10(nSNPs))) +
  geom_raster() +
  theme_classic() +
  labs(x="North America D. pulex MAF", 
       y="Europe D. pulex MAF",
       fill="SNP Counts^10") +
  #facet_wrap(~dt, nrow = 2) +
  #paletteer::scale_fill_paletteer_c("gameofthrones::targaryen") +
  scale_fill_viridis(option = "magma", begin = 0, end = 1) +
  theme(title = element_text(face="bold", size=20),
        legend.text = element_text(face="bold", size=26),
        legend.title = element_text(face="bold", size=24),
        legend.background = element_blank(),
        axis.text.x = element_text(face="bold", size=24),
        strip.text = element_text(face = "bold", size=16),
        axis.text.y = element_text(face="bold", size=24),
        axis.title.x = element_text(face="bold.italic", size=30),
        axis.title.y = element_text(face="bold.italic", size=30))}

### SFS ###
toti <- rbind(tot[classified=="shared_poly"],
              tot[Daphnia.pulex.NorthAmerica > 0.01 & classified=="polymorphic_1"],
              tot[Daphnia.pulex.Europe > 0.01 & classified=="polymorphic_2"])

buscoi <- rbind(busco[classified=="shared_poly"],
                busco[Daphnia.pulex.NorthAmerica > 0.01 & classified=="polymorphic_1"],
                busco[Daphnia.pulex.Europe > 0.01 & classified=="polymorphic_2"])

num <- data.table(toti %>% 
         group_by(classified) %>%
         summarize(n.snps.tot=n()))

# Colors
cols <- c("Poly Euro."="#999999", "Poly NAm."="#E69F00", "TSP"="#56B4E9")

# Minor alleles 
maf <- unique(c(seq(from=0, to=0.1, length.out=10),
              seq(from=0.1, to=0.4, length.out=5),
              0.45, 0.5, 0.51))

# Go through each minor allele
m <- foreach(i=1:length(maf), .combine = "rbind", .errorhandling = "remove") %do% {
  
  # Message
  print(paste("MAF filter", maf[i], i, sep=" "))
  
  # Filter
  buscoi <- data.table(rbind(busco[Daphnia.pulex.NorthAmerica >= maf[i] & 
                                   Daphnia.pulex.Europe >= maf[i] &
                                   classified=="shared_poly"],
                             busco[classified %in% c("polymorphic_1", "polymorphic_2") &
                                   Daphnia.pulex.NorthAmerica >= maf[i] | 
                                   Daphnia.pulex.Europe >= maf[i]]) %>% 
                  mutate(classified=case_when(classified %in% c("polymorphic_1", 
                                                                "polymorphic_2") ~ "Poly.",
                                              classified == "shared_poly" ~ "TSP")) %>% 
                  filter(simpleAnnot %in% c("NS", "Syn")))
  
  print(paste("Number of BUSCO TSP snps:", length(unique(buscoi[classified=="TSP"]$variant.id))))
  print(paste("Number of BUSCO Poly. snps:", length(unique(buscoi[classified=="Poly."]$variant.id))))
  
  toti <- data.table(rbind(tot[Daphnia.pulex.NorthAmerica >= maf[i] & 
                               Daphnia.pulex.Europe >= maf[i] & 
                                 classified=="shared_poly"],
                           tot[classified %in% c("polymorphic_1", "polymorphic_2")][
                               Daphnia.pulex.NorthAmerica >= maf[i] | 
                               Daphnia.pulex.Europe >= maf[i]]) %>% 
                  mutate(classified=case_when(classified %in% c("polymorphic_1", 
                                                                "polymorphic_2") ~ "Poly.",
                                              classified == "shared_poly" ~ "TSP")) %>% 
                  filter(simpleAnnot %in% c("NS", "Syn")))
  
  print(paste("Number of Genome-wide TSP snps:", length(unique(toti[classified=="TSP"]$variant.id))))
  print(paste("Number of Genome-wide Poly. snps:", length(unique(toti[classified=="Poly."]$variant.id))))
  
  # Using synonymous and non synonymous sites
  j <- table(buscoi$simpleAnnot,
             buscoi$classified)[c(2,1), c(1,2)] %>% 
    fisher.test()
  
  jk <- table(buscoi$simpleAnnot,
              buscoi$classified)[c(2,1), c(1,2)] 
  
  print(paste("BUSCO Odds ratio:", 
              j$estimate))
  
  k <- table(toti$simpleAnnot,
             toti$classified)[c(2,1), c(1,2)] %>% 
    fisher.test()
  
  kk <- table(toti$simpleAnnot,
             toti$classified)[c(2,1), c(1,2)]
  
  print(paste("Genome-wide Odds ratio:", 
              k$estimate))
  
  # Summarize odds
  jj <- data.table(rbind(data.table(snpset="BUSCO Genes",
                   maf=maf[i],
                   maf_upper=maf[i+1],
                   odd=j$estimate,
                   pval=j$p.value,
                   lower.p=j$conf.int[1],
                   upper.p=j$conf.int[2],
                   Polysyn=jk[1],
                   TSPsyn=jk[3],
                   Polyns=jk[2],
                   TSPns=jk[4]),
        data.table(snpset="Genome-wide",
                   maf=maf[i],
                   maf_upper=maf[i+1],
                   odd=k$estimate,
                   pval=k$p.value,
                   lower.p=k$conf.int[1],
                   upper.p=k$conf.int[2],
                   Polysyn=kk[1],
                   TSPsyn=kk[3],
                   Polyns=kk[2],
                   TSPns=kk[4])) %>%
        mutate(alpha = 1-(1/odd)) %>% 
        mutate(p_adjust = p.adjust(pval, method = "fdr")) %>% 
        mutate(signif = ifelse(.$p_adjust < 0.05, "Sig.","Not Sig.")))
  
  print(paste("alpha B:", 
              jj$alpha))
  print("")
  
  
  # Finish
  return(jj)
}

# Ods ratio calculation
ods_calc <- function(dtf, poly, maf.l, snpy) {
  #dtf=busco; poly="polymorphic_2"; maf.l=0; snpy="BUSCO Genes"
  
  # Filter
  dtf.tmp <- data.table(rbind(dtf[Daphnia.pulex.NorthAmerica >= maf.l & 
                                    Daphnia.pulex.Europe >= maf.l][classified=="shared_poly"],
                              dtf[classified %in% c(poly)][Daphnia.pulex.NorthAmerica >= maf.l | 
                                    Daphnia.pulex.Europe >= maf.l]) %>% 
                          mutate(classified=case_when(classified %in% c("polymorphic_1", 
                                                                        "polymorphic_2") ~ "Poly.",
                                                      classified == "shared_poly" ~ "TSP")) %>% 
                          filter(simpleAnnot %in% c("NS", "Syn")))
  
  j <- table(dtf.tmp$simpleAnnot,
             dtf.tmp$classified)[c(2,1), c(1,2)] %>% 
    fisher.test()
  
  jk <- table(dtf.tmp$simpleAnnot,
              dtf.tmp$classified)[c(2,1), c(1,2)] 
  
  data.table(snpset=snpy,
             maf=maf.l,
             poly=poly,
             odd=j$estimate,
             pval=j$p.value,
             lower.p=j$conf.int[1],
             upper.p=j$conf.int[2],
             Polysyn=jk[1],
             TSPsyn=jk[3],
             Polyns=jk[2],
             TSPns=jk[4])
}

oddy <- data.table(rbind(ods_calc(dtf=busco, poly = "polymorphic_1", maf.l=0, snpy="BUSCO Genes"),
           ods_calc(dtf=busco, poly = "polymorphic_2", maf.l=0, snpy="BUSCO Genes"),
           ods_calc(dtf=tot, poly = "polymorphic_1", maf.l=0, snpy="Genome-wide"),
           ods_calc(dtf=tot, poly = "polymorphic_2", maf.l=0, snpy="Genome-wide"))) %>% 
  mutate(poly=case_when(poly=="polymorphic_1" ~ "NAm. D. pulex",
                        poly=="polymorphic_2" ~ "European D. pulex"))

#saveRDS(oddy, file="/scratch/csm6hg/data/odds_ratio_conservative.rds")

# Odds ratio function
odds_ratio <- function(poly_syn, poly_ns, tsp_syn, tsp_ns) {
  # poly_syn=33; tsp_syn=1; poly_ns=18; tsp_ns=1
  
  fish <- matrix(c(poly_syn, poly_ns, 
           tsp_syn, tsp_ns), 
         nrow = 2, ncol = 2) %>% 
  fisher.test()
  
  data.table(odd=fish$estimate,
             pval=fish$p.value,
             lower.p=fish$conf.int[1],
             upper.p=fish$conf.int[2])
}

# Calculate Odds ratio by gene 
gene.odd <- data.table(gene %>% 
              group_by(snpset, gene) %>% 
              summarize(odds_ratio(poly_syn = Polysyn, poly_ns = Polyns,
                                   tsp_syn = TSPsyn, tsp_ns = TSPns)))
                         
# Combine 
gene.fin <- data.table(gene %>% 
             full_join(gene.odd, by=c('gene', 'snpset')) %>% 
             left_join(pro, by=c("gene")))
