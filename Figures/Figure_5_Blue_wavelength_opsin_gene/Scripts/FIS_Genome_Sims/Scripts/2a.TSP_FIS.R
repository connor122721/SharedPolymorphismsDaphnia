# Trans species polymorphism - FIS
# 3.18.2024
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)
  library(tidyverse)
  registerDoMC(10)

### Read in data
  dt1 <- data.table(readRDS("/scratch/csm6hg/data/AxC_lab_tot_het_withmet.rds"), exp="Lab AxC F1s")
  dt2 <- data.table(readRDS("/scratch/csm6hg/data/AxC_Nolab_subsamp_tot_het_withmet.rds"), exp="Wild AxC F1s Subsampled")
  dt3 <- data.table(readRDS("/scratch/csm6hg/data/AxC_Nolab_tot_het_withmet.rds"), exp="Wild AxC F1s")
  dt4 <- data.table(readRDS("/scratch/csm6hg/data/CxC_tot_het_withmet.rds"), exp="Lab Selfed C F1s")
  
  dat <- data.table(rbind(dt1, dt2, dt3, dt4))
  
  dat[,geno:="hom"]
  dat[af.dos>.01 & af.dos<.99, geno:="het"]
  dat[,rd:=(dos.ref + dos.alt)]
  dat[,posBin:=round(position, -4)]
  
  # Add metadata 
  fin <- data.table(read.csv("/project/berglandlab/connor/metadata/samples.9.8.22.csv"))
  dat <- data.table(dat %>% left_join(fin %>% select(c(id,WildSequenced,year,pondID)), by=c("Sample"="id")))
 
  # Restrict to conservative TSP set
  con.tsp <- data.table(readRDS("/scratch/csm6hg/data/unchanged_SNPs_across_assembly.RDS"))
  dat <- data.table(dat[ch %in% con.tsp$ch] %>% 
              mutate(snpset=case_when(classified=="shared_poly"~"TSP",
                                      classified %in% c("polymorphic_1","polymorphic_2") ~ "Not TSP")))
 
### Simulation functions
  
  # Calculate FIS given dosage
  fis <- function(nHet, nTot, p=.5) {
    numerator <- (nHet/nTot)
    denominator <- (nTot/(nTot-1))*(2*p*(1-p) - (nHet/nTot)/(2*nTot))
    1 - numerator/denominator
  }

  # Simulation given read depth
  sim_geno <- function(nAlt, nRef, rdSim=T) {
    # nAlt=10;nRef=0
    # first, sample genotype
    geno <- rbinom(1, 1, .5)
    if(geno==0) return("hom")
    if(geno==1 & rdSim==T) {
      af <- rbinom(1, nAlt+nRef, .5)/(nAlt+nRef)
      if(af>.01 & af<.99) {
        return("het")
      } else {
        return("hom")
          }
      }
    if(geno==1 & rdSim==F) return("het")
  }
  
  ### Summarize via gene class version
  setkey(dat, variant.id.y)

  # Go through 100 independent simulations
  o <- foreach(i=0:100) %dopar% {
    message(paste(i, i, sep=" / "))
    if(i==0) {
      
      # Empirical Data
      dat.ag <- dat[rd>0 & !col%in%c("intergenic_region", "downstream_gene_variant", "upstream_gene_variant"),
              list(FIS=fis(nHet=sum(geno=="het"),
                            nTot=sum(geno=="het") + sum(geno=="hom"),
                            p=.5),
                    rd=round(mean(rd), -1),
                    perm=i, rdSim=T),
               list(variant.id.y, gene, exp)]

    } else {
      
      # Simulate w/read depth
      dat.sim.a <- dat[rd>0 & !col%in%c("intergenic_region", "downstream_gene_variant", "upstream_gene_variant"),
                        list(sim.geno=sim_geno(nAlt=dos.alt, nRef=dos.ref), 
                             rd=dos.alt+dos.ref), 
                        list(variant.id.y, Sample, gene, exp)]
      
      # Output blueopsin sims
      #saveRDS(dat.sim.a[gene=='Daphnia11806-RA'], file = "/scratch/csm6hg/data/Crosses.daphnia11806-ra.sims.rds")
      dat.ag <- dat.sim.a[rd>0,list(FIS=fis(nHet=sum(sim.geno=="het"),
                                    nTot=sum(sim.geno=="het") + sum(sim.geno %like% "hom"),
                            p=.5),
                    rd=round(mean(rd), -1),
                    perm=i, rdSim=T),
               list(variant.id.y, gene, exp)]
    }
    return(dat.ag)
  }
  
  # Bind output
  o.b <- rbindlist(o, fill=T)
  
  # Summarize
  o.ag <- o.b[,list(FIS=mean(FIS, na.rm = T), 
                  med=median(FIS, na.rm = T),
                  sd=sd(FIS, na.rm = T),
                  lci=quantile(FIS, probs = 0.025, na.rm = T),
                  uci=quantile(FIS, probs = 0.975, na.rm = T),
                  .N), 
            list(gene, perm, exp)]
  o.ag[,se:=sd/sqrt(N)]
  
  saveRDS(o.ag, file = "/scratch/csm6hg/data/Crosses_geneclass_conservative_100boot_raw.rds")
  
  # Simulation given read depth - and all genotypes are heterozygous
  
  sim_geno_het <- function(nAlt, nRef, rdSim=T) {
    # nAlt=10;nRef=0
    # first, sample genotype
    geno <- 1
    if(geno==0) return("hom ref")
    if(geno==2) return("hom alt")
    if(geno==1 & rdSim==T) {
      af <- rbinom(1, nAlt+nRef, .5)/(nAlt+nRef)
      if(af>.01 & af<.99) {
        return("het")
      } else {
        if(af<0.5)
          return("hom ref")
        else{
          return("hom alt")
        }
      }
    }
    if(geno==1 & rdSim==F) return("het")
  }
  
  # Simulate given read depth and at a 
  o.hom <- foreach(i=0:100) %dopar% {
    message(paste(i, i, sep=" / "))
    
    # Simulate w/read depth
    dat.sim.a <- dat[rd>0 & !col%in%c("intergenic_region", "downstream_gene_variant", "upstream_gene_variant"),
                     list(sim.geno=sim_geno_het(nAlt=dos.alt, nRef=dos.ref), 
                          rd=dos.alt+dos.ref,
                          perm=i), 
                     list(variant.id.y, Sample, snpset, gene, exp)]
  }
  
  o.test <- rbindlist(o.hom)
  
  blop.sum <- o.test[gene=="Daphnia11806-RA"] %>% 
    group_by(variant.id.y, snpset, exp) %>% 
    mutate(num=n()) %>%
    group_by(variant.id.y, snpset, exp, sim.geno) %>% 
    summarize(prop=n()/num)
  
  
  blop <- data.table(o.test[gene=="Daphnia11806-RA"] %>% 
                       group_by(variant.id.y, snpset, sim.geno, exp) %>% summarize(y=n()) %>% left_join(blop.sum) %>% 
                       mutate(y.prop=y/num*100,
                              snpset=case_when(snpset=="TSP"~"Sim. TSP",
                                               snpset=="Not TSP"~"Sim. Not TSP"),
                              variant.id=paste("test12", snpset),
                              all=case_when(sim.geno=="het"~1,
                                            sim.geno=="hom ref"~0,
                                            sim.geno=="hom alt"~2)))
  
  # Output
  saveRDS(blop, file = "/scratch/csm6hg/data/AxC_geneclass_conservative_100boot_Wild_sim_fullhet.rds")
  
  # Plot
  blop %>%
    ggplot(., 
           aes(x=as.factor(all), 
               y=y.prop, 
               color=snpset,
               group=variant.id)) +
    geom_line(size=1.5) +
    geom_point() 
  
  
  # Simulation given read depth
  sim_geno <- function(nAlt, nRef, rdSim=T) {
    # nAlt=10;nRef=0
    # first, sample genotype
    geno <- rbinom(1, 2, .5)
    if(geno==0) return("hom ref")
    if(geno==2) return("hom alt")
    if(geno==1 & rdSim==T) {
      af <- rbinom(1, nAlt+nRef, .5)/(nAlt+nRef)
      if(af>.01 & af<.99) {
        return("het")
      } else {
        if(af<0.5)
          return("hom ref")
        else{
          return("hom alt")
        }
      }
    }
    if(geno==1 & rdSim==F) return("het")
  }
  
  
  
  

### Standard version - using genome-wide SNPs ###
  setkey(dat, variant.id)

    o <- foreach(i=0:10)%dopar%{
      message(paste(j, i, sep=" / "))
      if(i==0) {
        dat.ag <-
        dat[rd>0,list(FIS=fis(nHet=sum(geno=="het"),
                              nTot=sum(geno=="het") + sum(geno=="hom"),
                              p=.5),
                      rd=round(mean(rd), -1),
                      perm=i, rdSim=T, boot=j),
                 list(variant.id, snpset, col)]
      } else {
        dat.sim.a <- dat[rd>0,list(sim.geno=sim_geno(nAlt=dos.alt, nRef=dos.ref), rd=dos.alt+dos.ref), list(variant.id, Sample, snpset)]
        dat.ag.a <-
        dat.sim.a[rd>0,list(FIS=fis(nHet=sum(sim.geno=="het"),
                                  nTot=sum(sim.geno=="het") + sum(sim.geno=="hom"),
                              p=.5),
                      rd=round(mean(rd), -1),
                      perm=i, rdSim=T, boot=j ),
                 list(variant.id, snpset, col)]


       dat.sim.b <- dat[rd>0, list(sim.geno=sim_geno(nAlt=dos.alt, nRef=dos.ref, rdSim=F), rd=dos.alt+dos.ref), list(variant.id, Sample, snpset)]
       dat.ag.b <-
       dat.sim.b[rd>0,list(FIS=fis(nHet=sum(sim.geno=="het"),
                                 nTot=sum(sim.geno=="het") + sum(sim.geno=="hom"),
                             p=.5),
                     rd=round(mean(rd), -1),
                     perm=i, rdSim=F, boot=j),
                list(variant.id, snpset, col)]

        dat.ag <- rbind(dat.ag.a, dat.ag.b)
      }
      return(dat.ag)
    }
    o <- rbindlist(o, fill=T)
    return(o)
  }
  o <- rbindlist(o, fill=T)

  o.ag <- o[,list(t=t.test(FIS~snpset)$statistic,
                  FIS_tsp=mean(FIS[snpset=="TSP"]),
                  FIS_gw=mean(FIS[snpset!="TSP"])), list(perm, rdSim, boot)]

  ggplot(data=o.ag, aes(x=as.factor(perm), y=t)) + geom_boxplot() + facet_wrap(~rdSim)


  o.ag2 <- o[,list(FIS=mean(FIS), sd=sd(FIS), .N), list(perm, snpset, rdSim, boot)]
  o.ag2[,se:=sd/sqrt(N)]

  ggplot(data=o.ag2, aes(x=snpset, y=FIS, color=as.factor(perm==0), group=interaction(perm, boot))) +
  geom_line() + geom_point() +
  geom_segment(aes(x=snpset, xend=snpset, y=FIS-2*se, yend=FIS+2*se)) + facet_grid(as.factor(perm==0)~rdSim)



