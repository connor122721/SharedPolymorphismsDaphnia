# Connor Murray
# 11.4.2024
# Run eSMC on Daphnia dataset

# Install
# git clone https://github.com/TPPSellinger/eSMC
# path="/home/csm6hg/eSMC/eSMC_2.0.5.tar.gz" # Path to the dowloaded eSMC package
# devtools::install_local(path)
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Library
library(eSMC)
library(data.table)

# Working directory
setwd("/project/berglandlab/connor/backup_project/")
#doParallel::registerDoParallel(4)

# From consule
args<-commandArgs(TRUE)
samplei=args[1] # samplei="SRR7592114"

##############
# Parameters # 
##############

# Please Fill all value
mu=5.69*10^(-9)   # Mutation rate per position per generation 
gen=15 # sexual generations
scaling=1/gen # Ratio of recombination 
r=8*10^(-8)*scaling   # recombination rate per position per generation 
rho=r/mu # ratio recombination/mutation 
M=2 # Number of haplotypes
###################################
# Set one to TRUE (rest to FALSE) #
###################################
ER=F # True to estimate recombination rate
SF=F # True to estimate selfing rate
SB=T # True to estimate germination rate

# Set boundaries
BoxB=c(0.05, 1) #  min and max value of germination rate 
Boxs=c(0, 0.99) #  min and max value of selfing rate 

############
# Get data #
############

# example 
NC=12 # Number of analysed scaffold

# We build the data for each chromosome 
full_data=list()

# We build the data of the 12 scaffolds
chroms=fread("metadata/goodChrom.txt", header = F)

# Merge all multihetsep files
for(scaffold in 1:NC){
  #scaffold=1
  filename=paste(chroms[scaffold],"/",samplei,".",chroms[scaffold],".","phase.filt.samp", sep="")
  data=Get_real_data(path="/project/berglandlab/connor/backup_project/chapter1/msmc",M,filename,delim="\t")
  full_data[[scaffold]]=data
}

# Run eSMC - takes ~10mins each chrom
result=eSMC(n=30,rho=rho,data,BoxB=BoxB,Boxs=Boxs,SB=SB,SF=SF,Rho=ER,Check=F,NC=1) 
    
# Assemble summary stats
dti <- data.table(x = result$Tc*(result$mu/mu)*gen,
                  y = result$Xi*(result$mu/(2*mu)),
                  L = result$L, 
                  Tc = result$Tc,
                  gens = gen,
                  mu = result$mu, 
                  Sample = samplei, 
                  chromosome = chromi)

# Output
saveRDS(dti, 
        file = paste("/project/berglandlab/connor/daphnia_eSMC_output/",
        samplei, ".", chromi, ".rds", sep=""))
