################################
#### Run Segmentation Model ####
################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_behavior")

set.seed(2)

library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')


setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

#load and manipulate data
dat<- read.csv("CRW_HC_sim.csv", as.is = T)  #for hard-clustering sim
# dat<- read.csv("CRW_MM_sim.csv", as.is = T)  #for mixed-membership sim
dat$dt<- 3600
dat$id<- 1

dat.list<- df.to.list(dat=dat)
names(dat)[3:4]<- c("dist","rel.angle")

#filter data for tstep of interest
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval

#define bin number and limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(dat[dat$dt == 3600,]$dist, na.rm = T)
upper90.thresh=as.numeric(quantile(dat[dat$dt == 3600,]$dist, 0.90, na.rm=T)) 
dist.bin.lims=seq(from=0, to=upper90.thresh, length.out = 5)
dist.bin.lims=c(dist.bin.lims, max.dist)  #6 bins

#assign bins to obs
for (i in 1:length(behav.list)) {
  behav.list[[i]]<- behav.list[[i]] %>% assign.dist.bin(dist.bin.lims = dist.bin.lims,
                                                        max.dist = max.dist) %>%
    assign.rel_angle.bin(angle.bin.lims = angle.bin.lims)
}

behav.list2<- lapply(behav.list, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment



## Run RJMCMC

#prior
alpha = 1

ngibbs = 40000
dat.res<- behavior_segment(data = behav.list2, ngibbs = ngibbs, nbins = c(6,8), alpha = alpha)
#takes 9 min for 40000 iterations


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, brk.cols = 99)  #brk.cols is max matrix cols

## Heatmaps
plot.heatmap(data = behav.list, nbins = c(6,8), brkpts = brkpts, dat.res = dat.res, type = "behav")


## Compare True vs Modeled Breakpoints
model.brkpts<- na.omit(t(brkpts[-1])) %>% as.vector()
true.brkpts<- which(diff(as.numeric(as.factor(dat$true.behav))) != 0)  #for hard-clustering sim
# true.brkpts<- which(diff(as.numeric(as.factor(dat$behav_coarse))) != 0)  #for mixed-membership sim
all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts), type = rep(c("True","Model"),
                                                                        c(length(true.brkpts),
                                                                          length(model.brkpts))))

ggplot(all.brkpts, aes(x=brks, y=type)) +
  geom_point(size=2) +
  theme_bw() +
  labs(x="Time", y="Type") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10))

dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF


#export results for hard-clustering and mixed-membership simulations
# write.csv(dat_out, "CRW_HC_tsegs.csv", row.names = F)
# write.csv(dat_out, "CRW_MM_tsegs.csv", row.names = F)