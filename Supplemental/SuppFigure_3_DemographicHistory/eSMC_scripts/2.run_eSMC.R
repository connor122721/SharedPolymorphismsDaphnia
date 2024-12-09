# Connor Murray
# 11.4.2024
# Run eSMC on Daphnia dataset

# Install
# git clone https://github.com/TPPSellinger/eSMC
# path="/home/csm6hg/eSMC/eSMC_2.0.5.tar.gz" # Path to the dowloaded eSMC package
#devtools::install_local(path)

# Library
library(eSMC)
library(foreach)
library(data.table)

# Working directory
setwd("/project/berglandlab/connor/backup_project/")
#doParallel::registerDoParallel(4)

# From consule
args<-commandArgs(TRUE)
filei=args[1]
samplei=args[2]
chromi=args[3]

##############
# Parameters # 
##############

# Please Fill all value
mu=4.33*10^(-9)   # Mutation rate per position per generation 
r=8*10^(-8)   # recombination rate per position per generation 
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
gen=5 # generation

# We build the data for each chromosome 

# Read in all scaffolds for sample
data = Get_real_data(M = M, filename = filei, delim = "\t")
    
# Run eSMC - takes ~10mins each chrom
result=eSMC(n=30, 
            rho=rho, 
            data, 
            BoxB=BoxB, 
            Boxs=Boxs, 
            SB=SB,
            SF=SF, 
            Rho=ER, 
            Check=F, 
            NC=1) 
    
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
