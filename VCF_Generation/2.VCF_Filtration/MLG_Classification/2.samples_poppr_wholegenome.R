# Classifying clonal lineages
# Connor Murray 9.22.2021

# module load intel/18.0 intelmpi/18.0 R/3.6.3; R
# module load  gcc/7.1.0  openmpi/3.1.4;
# ijob -A berglandlab --mem=50G -p standard -c 10
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(poppr)
library(vcfR)
library(dartR)
library(parallel)
library(adegenet)

# Working directory
setwd("/project/berglandlab/connor/new_vcf2")

# Metadata
fin.dt <- data.table(read.csv(file = "../metadata/samples.1.31.22.csv", header = T))
                  
# Genlight objects
gen <- list.files(pattern = "*.daphnia.filtered.chr.busco.genlight.rds$", full.names = T)

# Register cores
#doParallel::registerDoParallel(cores = 10)

# Cbind chromosomal genlights
dt <- foreach(i=1:12, .combine = "cbind", .errorhandling = "remove") %do% {

  # Start message
  print(paste("Starting chromosome", i, sep=" "))

  # Merge all genlight chunks
  vcf <- c(gen[i] %>% 
              map(readRDS))[[1]]
  
  print(length(unique(vcf@ind.names)))
  
  # Remove unwanted individuals
  vcf2 <- vcf[indNames(vcf) %in% as.character(fin.dt$Sample)]
  vcf2@ploidy <- as.integer(2)
  vcf2<-as.matrix(vcf2)
  
  # Finish message
  print(paste("Finished chromosome:", i, sep=" "))

  # Return function
  return(vcf2)
  
}
dt1 <- new("genlight", dt)
       
# Save genlight object
saveRDS(dt1, file = "daphnia.filtered.chr.busco.genlight.rds")
dt1 <- readRDS("daphnia.filtered.chr.busco.genlight.rds")

# All diploid individuals
ploidy(dt1) <- 2

# Convert to snpclone object
xsnp <- as.snpclone(dt1, parallel = T, n.cores = 10)

# Save genlight object
saveRDS(xsnp, file = "daphnia.filtered.chr.busco.genlight.snpclone.rds")
xsnp <- readRDS("daphnia.filtered.chr.busco.genlight.snpclone.rds")

# Extract individuals and MLGs
snp.dt <- data.table(data.table(Sample=indNames(xsnp),
                                mlg=mll(xsnp)) %>% 
             left_join(fin.dt %>% select(Sample, Species, Continent, SC), 
            by=c("Sample")))

# Add populations to snpclone object
pop(xsnp) <- paste(str_replace(snp.dt$Species,pattern = " ", replacement = "."), 
                   snp.dt$Continent, sep=".")

# Exract mlgs
dt <- data.table(n=t(mlg.table(xsnp, plot=F)),
                 mlg=colnames(mlg.table(xsnp, plot=F))) 

dt <- data.table(dt %>% pivot_longer(cols = colnames(dt)[1:7],
                                     names_prefix = "n.", 
                                     names_to="pop",
                                     values_to="n"))

# Rank order based on genome counts 
# dt$mlg <- factor(dt$mlg, levels=as.factor(data.table(dt %>% arrange(desc(n)))$mlg))

# Number of clones per mlg
a <- dt[n>1] %>% 
  ggplot(.) +
  #facet_wrap(~pop, scales="free_x") +
  geom_col(aes(x=reorder(mlg, -n), y=log10(n), fill=pop), 
           position = "identity") +
  labs(x="", 
       y="",
       fill="Population",
       title="Default MLG threshold") +
  theme_classic() +
  scale_y_continuous(breaks=c(log10(2),1,log10(25),log10(50),2,
                              log10(200),log10(300)),
                     labels=c(2,10,25,50,100,200,300)) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.position = c(0.6, 0.8), 
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

# Gather thresholds to test
thresholds <- mlg.filter(xsnp, 
                         distance = bitwise.dist, 
                         stats = "thresholds",
                         threshold = 8, 
                         threads = 10)

# We can use these thresholds to find an appropriate cutoff
pcut <- cutoff_predictor(thresholds)
cutoff

# Apply filter
mlg.filter(xsnp, distance = bitwise.dist) <- pcut
nmll(xsnp)

# Exract mlgs
dt <- data.table(n=t(mlg.table(xsnp, plot=F)),
                 mlg=colnames(mlg.table(xsnp, plot=F))) 

dt <- data.table(dt %>% pivot_longer(cols = colnames(dt)[1:7],
                                     names_prefix = "n.", 
                                     names_to="pop",
                                     values_to="n"))

# Reset mlg filter
mll(xsnp) <- "original" 

# Different distances will give different results
nres <- filter_stats(xsnp,
                     distance = bitwise.dist,
                     threshold =  1e+06 + .Machine$double.eps^0.5,
                     plot = TRUE,
                     missing = "ignore", 
                     hist = NULL,
                     nclone=c(162))

# Best threshold
print(thresh <- cutoff_predictor(nres$farthest$THRESHOLDS))
mlg.filter(xsnp, missing = "ignore") <- thresh
xsnp
  
# Exract mlgs
dt <- data.table(n=t(mlg.table(xsnp, plot=F)),
                 mlg=colnames(mlg.table(xsnp, plot=F))) 

dt <- data.table(dt %>% pivot_longer(cols = colnames(dt)[1:7],
                                     names_prefix = "n.", 
                                     names_to="pop",
                                     values_to="n"))

# Rank order based on genome counts 
#dt$mlg <- factor(dt$mlg, levels=as.factor(data.table(dt %>% arrange(desc(n.Total)))$mlg))

# Number of clones per mlg
b <- dt[n>1] %>% 
  ggplot() +
  geom_col(aes(x=reorder(mlg, -n), y=log10(n), fill=pop), 
           position = "identity") +
  labs(x="Unique MLG", 
       y="", 
       title="0.006 threshold: 162 MLGs") +
  theme_classic() +
  scale_y_continuous(breaks=c(log10(2),1,log10(25),log10(50),2,
                              log10(200),log10(300)),
                     labels=c(2,10,25,50,100,200,300)) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.position = "none", 
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

# Multiple plots
plot_grid(a,b, ncol = 1, align = "hv")

# Clonal diversity
monstat <- diversity_stats(mlg.table(xsnp))
diversity_ci(mlg.table(xsnp), n = 10000L, ci=97.5, center = T, raw = FALSE)

# Resampling affects expected MLGs
rarecurve(mlg.table(xsnp, plot=F, 
                    sublist = c("Daphnia.pulex.Europe", "Daphnia.pulex.NorthAmerica")),
          ylab="Number of expected MLGs", 
          sample=min(rowSums(mlg.table(xsnp, plot=F, 
                    sublist = c("Daphnia.pulex.Europe", "Daphnia.pulex.NorthAmerica")))),
          font = 2, cex = 1)

# Distance matrix
dist.snp <- bitwise.dist(xsnp)

length(dt[pop=="Daphnia.pulex.Europe"][n>0]$mlg)
length(dt[pop=="Daphnia.pulex.NorthAmerica"][n>0]$mlg)

# Minimum spanning networks
dm <- poppr.msn(xsnp, 
                distmat = dist.snp,
                sublist = c("Daphnia.pulex.Europe"))

dm <- poppr.msn(xsnp, 
                distmat = dist.snp,
                sublist = c("Daphnia.pulex.NorthAmerica"))

dm <- poppr.msn(xsnp, 
                distmat = dist.snp,
                sublist = c("Daphnia.pulicaria.Europe", "Daphnia.pulicaria.NorthAmerica"))

dm <- poppr.msn(xsnp, 
                distmat = dist.snp,
                sublist = c("Daphnia.pulicaria.Europe", "Daphnia.pulex.NorthAmerica"))

dm <- poppr.msn(xsnp, 
                distmat = dist.snp,
                sublist = c("Daphnia.obtusa.Europe", 
                            "Daphnia.obtusa.NorthAmerica", 
                            "Daphnia.pulex.NorthAmerica"))

# Unique preidentified superclones from Europe
euro.sc <- unique(snp.dt[Continent=="Europe"][!SC==""]$SC)

# Compare to european clonal classifications
euro.valif <- foreach(i=1:length(euro.sc), .combine = "rbind") %do% {
  
    print(i)
    p <- snp.dt[Continent=="Europe"][SC==euro.sc[i]]
    p.list <- data.table(unique.mlg=length(unique(p$mlg)),
                         n.min.mlg=min(table(p$mlg)),
                         n.max.mlg=max(table(p$mlg)),
                         n.mean.mlg=mean(table(p$mlg)),
                         SC=euro.sc[i],
                         SC.size=length(p$SC),
                         iteration=i)
    
  return(p.list)
    
}

a <- euro.valif %>% 
  group_by(SC) %>% 
  summarise(prop.right=n.max.mlg/SC.size,
            unique.mlg,
            SC.size) %>% 
  filter(SC.size>1, unique.mlg>1) %>% 
  ggplot(aes(x=reorder(SC, -SC.size), y=prop.right*100)) +
  geom_col(aes(fill=unique.mlg)) +
  geom_label(aes(y=0, label=SC.size)) +
  scale_fill_viridis(option = "inferno") + 
  labs(x="Super clone", 
       y="Proportion of samples within major MLG (%)",
       fill="Number of MLGs") +
  theme_classic() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.background = element_blank(),
        axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))


b <- euro.valif %>% 
  group_by(SC) %>% 
  summarise(prop.right=n.max.mlg/SC.size,
            unique.mlg,
            SC.size) %>% 
  filter(SC.size>1, unique.mlg>1) %>% 
  ggplot(aes(x="", y=prop.right*100, label=SC)) +
  geom_boxplot() +
  geom_label_repel(aes(fill=unique.mlg), color="white") +
  scale_fill_viridis(option = "inferno") + 
  ylim(0,100) +
  labs(x="Super clone", 
       y="",
       fill="Number of MLGs") +
  theme_classic() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=15),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

plot_grid(a, b, rel_widths = c(2,1))
