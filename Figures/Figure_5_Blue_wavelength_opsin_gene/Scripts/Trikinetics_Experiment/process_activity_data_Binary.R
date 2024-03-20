# Trans species polymorphism - FIS
# 3.12.2024
# ijob -c 1 --mem=20G -p standard -A berglandlab
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(lubridate)
library(lme4)
library(emmeans)
library(ggplot2)

setwd("/scratch/csm6hg/data")

### raw data
  exp1 <- as.data.table(readRDS("1027newdata.rds"))
  exp2 <- as.data.table(readRDS("112newdata.rds"))
  exp3 <- as.data.table(readRDS("1205newdata.rds"))
  
  exp1[,exp:="10_27_2023"]
  exp2[,exp:="11_02_2023"]
  exp3[,exp:="12_05_2023"]

### some renaming
  setnames(exp1, c("V3", "V2"), c("time", "day"))
  setnames(exp2, c("V3", "V2"), c("time", "day"))
  setnames(exp3, c("V3", "V2"), c("time", "day"))

  dat <- rbindlist(list(exp1, exp2, exp3), fill=T)
  #dat[,time_hms:=hms(time)]
  dat[,day:=dmy(day)]
  dat[,hour:= as.numeric(tstrsplit(time, ":")[[1]])]
  dat[,min:=  as.numeric(tstrsplit(time, ":")[[2]])]
  dat[,sec:=  as.numeric(tstrsplit(time, ":")[[3]])]
  dat[,time_sec:=hour*60*60 + min*60 + sec]
  dat[,jday:=yday(day)]
  dat[,experiment_time:=jday+time_sec/(24*60*60)]
  dat[,ID:=paste(file_category, position, light_trt, clone, exp, sep=";")]

### determine first observation & elapsed time
  dat.first <- dat[,list(start_time=min(experiment_time), .N), list(ID)]

  dat <- merge(dat, dat.first, by="ID")
  dat[,elapsed_time:=experiment_time - start_time]
  dat[,elapsed_time_hours:=elapsed_time*24]

### trim outliers & extra time
  dat <- dat[elapsed_time_hours>=.5 & elapsed_time_hours<=12]

### does everything look okay? yes.
  ggplot(data=dat, aes(x=activity)) + geom_histogram()

### basic summary plot (!!!)
  dat.activity.ag <- dat[,list(mean=mean(activity>0), .N), list(etm=round(elapsed_time_hours*6)/6, light_trt, genotype)]
  dat.activity.ag[,se:=1.96*sqrt(mean*(1-mean)/N)]
  dat.activity.ag[,uci:=mean-se]
  dat.activity.ag[,lci:=mean+se]
  dat.activity.ag[,light_trt:=factor(light_trt, levels=c("white", "blue", "dark"))]

  save(dat.activity.ag, dat, file="/Users/alanbergland/Documents/GitHub/misc/connor/activity/activity_data.Rdata")

  ggplot(data=dat.activity.ag, aes(x=etm, y=mean, color=genotype, group=genotype)) +
  geom_smooth(method="glm") +
  facet_grid(~light_trt) + theme_bw()


### time model - binomial
  dat[,actbin:=activity>0]
  t1 <- glmer(actbin~elapsed_time_hours + (1|ID) + (1|clone) + (1|exp), data=dat, family=binomial(), nAGQ = 0)
  t2 <- glmer(actbin~elapsed_time_hours + light_trt + (1|ID) + (1|clone) + (1|exp), data=dat, family=binomial(), nAGQ = 0)
  t3 <- glmer(actbin~elapsed_time_hours + light_trt + genotype + (1|ID) + (1|clone) + (1|exp), data=dat, family=binomial(), nAGQ = 0)
  t4 <- glmer(actbin~elapsed_time_hours + light_trt * genotype + (1|ID) + (1|clone) + (1|exp), data=dat, family=binomial(), nAGQ = 0)
  t5 <- glmer(actbin~elapsed_time_hours * light_trt * genotype + (1|ID) + (1|clone) + (1|exp), data=dat, family=binomial(), nAGQ = 0)
  anova(t1, t2, t3, t4, t5)
  save(t1, t2, t3, t4, t5, dat, file="/Users/alanbergland/Documents/GitHub/misc/connor/activity/lmer_repeated_measures.Rdata")
  save(t5, file="/Users/alanbergland/Documents/GitHub/misc/connor/activity/lmer_repeated_measures_t5.Rdata")
  Anova(t5, type="III")

### summarize across time
  dat.sum <- dat[,list(mean_activity=log(1+mean(activity)),
                                            bin_activity=as.integer(sum(activity>0)),
                                            sum_activity=sum(activity),
                                            .N, maxTime=max(elapsed_time_hours)),
                  list(light_trt, genotype, exp, ID, clone)]

  table(is.integer(dat.sum$bin_activity*dat.sum$N))

### binomial linear mixed model
  dat.sum[,frac:=bin_activity/N]
  hist(dat.sum$bin_activity, breaks=100)
  t1.glmer <- glmer(frac~light_trt+(1|exp)+(1|clone),          data=dat.sum, weights=dat.sum$N, family=binomial(), nAGQ = 0)
  t2.glmer <- glmer(frac~light_trt+genotype+(1|exp)+(1|clone), data=dat.sum, weights=dat.sum$N, family=binomial(), nAGQ = 0)
  t3.glmer <- glmer(frac~light_trt*genotype+(1|exp)+(1|clone), data=dat.sum, weights=dat.sum$N, family=binomial(), nAGQ = 0)
  anova(t1.glmer, t2.glmer, t3.glmer)
  plot(t3.glmer)
  Anova(t3.glmer)
  save(t1.glmer, t2.glmer, t3.glmer, file="/Users/alanbergland/Documents/GitHub/misc/connor/activity/lmer_summarized.Rdata")
