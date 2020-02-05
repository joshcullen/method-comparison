#######################
#### Run LDA Model ####
#######################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_LDA_behavior")

set.seed(1)

library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)


source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

#get data
dat<- read.csv("CRW_HC_tsegs.csv", as.is = T)  #hard clustering sim
# dat<- read.csv("CRW_MM_tsegs.csv", as.is = T)  #mixed-membership sim
dat.list<- df.to.list(dat)  #for later behavioral assignment
nbins<- c(6,8)  #number of bins per param (in order)
dat_red<- dat %>% dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- get.summary.stats_behav(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
ndata.types=length(nbins)

#run Gibbs sampler
res=LDA_behavior_gibbs(dat=obs, gamma1=gamma1, alpha=alpha,
                       ngibbs=ngibbs, nmaxclust=nmaxclust,
                       nburn=nburn, ndata.types=ndata.types)

#Check traceplot of log marginal likelihood
plot(res$loglikel, type='l')

#Extract and plot proportions of behaviors per time segment
theta.post<- res$theta[(nburn+1):ngibbs,]  #extract samples from posterior
theta.estim<- theta.post %>% apply(2, mean) %>% matrix(nrow(obs), nmaxclust) #calc mean of posterior
boxplot(theta.estim, xlab="Behavior", ylab="Proportion of Total Behavior")

#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
round(apply(theta.estim, 2, sum)/nrow(theta.estim), digits = 3)


## Viz histograms from model
behav.res<- get_behav_hist(res = res, dat_red = dat_red)
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors

#Plot histograms of proportion data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(param ~ behav, scales = "fixed")



## Viz behavior over time
#Assign behaviors (via theta) to each time segment
theta.estim<- apply(theta.estim[,1:3], 1, function(x) x/sum(x)) %>% t()  #normalize probs for only first 3 behaviors being used
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)

#DEFINE BEHAVIOR ORDER BASED ON HISTOGRAMS
names(theta.estim)<- c("id", "tseg","Resting","Exploratory","Transit")  
nobs<- data.frame(id = obs$id, tseg = obs$tseg, n = apply(obs[,3:8], 1, sum)) #calc obs per tseg using SL bins (more reliable than TA)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim2<- aug_behav_df(dat = dat[-1,] %>% mutate(date=1:nrow(dat[-1,])),
                            theta.estim = theta.estim, nobs = nobs)

#Change into long format
theta.estim.long<- theta.estim2 %>% gather(key, value, -id, -tseg, -time1, -date)
names(theta.estim.long)[5:6]<- c("behavior","prop")
theta.estim.long$behavior<- factor(theta.estim.long$behavior,
                                   levels = c("Resting","Exploratory","Transit"))


#generate long form of true behavior for HARD CLUSTERING SIM
true.behavior<- matrix(0, 2500, 3) %>% data.frame(., time1 = 1:2500)
names(true.behavior)[1:3]<- c("Resting","Exploratory","Transit")
ind<- factor(dat$true.behav[-1], levels = c("Resting","Exploratory","Transit")) %>% as.numeric()
for(i in 1:nrow(true.behavior)) {
  true.behavior[i, ind[i]]<- 1
}
true.behavior.long<- true.behavior %>% gather(key, value, -time1)
names(true.behavior.long)[2:3]<- c("behavior","prop")
true.behavior.long$behavior<- factor(true.behavior.long$behavior,
                                     levels = c("Resting","Exploratory","Transit"))



#generate long form of true behavior for MIXED-MEMBERSHIP SIM
true.behavior<- matrix(0, 5000, 3) %>% data.frame(., time1 = 1:5000)
names(true.behavior)[1:3]<- c("Resting","Exploratory","Transit")
tseg<- rep(1:50, each = 100)
tmp<- data.frame(behav = as.numeric(factor(dat$behav_fine[-1],
                                           levels = c("Resting","Exploratory","Transit"))),
                 tseg = tseg)
for(i in 1:(length(dat$behav_fine[-1])/100)) {
  tmp1<- tmp %>% filter(tseg == i) %>% dplyr::select(behav) %>% table()/100
  mat<- matrix(0, 1, 3)
  mat[,as.numeric(names(tmp1))]<- tmp1
  mat1<- matrix(mat, 100, 3, byrow = T)
  true.behavior[which(tmp$tseg == i), 1:3]<- mat1
}
true.behavior.long<- true.behavior %>% gather(key, value, -time1)
names(true.behavior.long)[2:3]<- c("behavior","prop")
true.behavior.long$behavior<- factor(true.behavior.long$behavior,
                                     levels = c("Resting","Exploratory","Transit"))

ggplot(theta.estim.long) +
  geom_line(aes(x=time1, y=prop, color = behavior), size = 1) +
  geom_line(data = true.behavior.long, aes(x=time1, y=prop, color=behavior), size = 0.5) +
  labs(x = "\nObservation", y = "Proportion of Behavior\n") +
  scale_color_viridis_d("Behavior", guide=F) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12)) +
  facet_wrap(~behavior)


#assign  behavior from sim to data
dat2<- assign_behav(dat.list = dat.list, theta.estim2 = theta.estim2)
dat2$behav<- factor(dat2$behav, levels = c("Resting","Exploratory","Transit"))
dat2[,c("behav","prop")]<- dat2[c(2501,1:2500),c("behav","prop")]  #FOR HC SIM
# dat2[,c("behav","prop")]<- dat2[c(5001,1:5000),c("behav","prop")]  #FOR MM SIM

ggplot() +
  geom_path(data = dat2, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat2[-1,], aes(x, y, fill=behav), size=2.5, pch=21,
             alpha=dat2$prop[-nrow(dat2)]) +
  scale_fill_viridis_d("Behavior") +
  geom_point(data = dat2[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat2[nrow(dat2),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  coord_equal()


#export results
# write.csv(dat2, "Modeled HC Sim Tracks w Behav.csv", row.names = F)  #for hard-clustering sim
# write.csv(dat2, "Modeled MM Sim Tracks w Behav.csv", row.names = F)  #for mixed-membership sim