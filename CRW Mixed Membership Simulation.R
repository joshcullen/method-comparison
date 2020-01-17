#####################################
#####################################
#### Mixed Membership Simulation ####
#####################################
#####################################

library(tidyverse)
library(circular)

source('Simulation Functions.R')

### Simulate full track ###

#define behaviors and randomly sample 50 (for 50 time segments)
#weight probs so that behavior 1 (Resting) occurs 50%, behavior 2 (Exploratory) occurs 35%, and behavior 3 (Transit) occurs 15%

set.seed(2)

#create vector of dominant behaviors per time segment (50 segments)
behav<- sample(c(1,2,3), 50, replace = TRUE, prob = c(0.5, 0.35, 0.15))
table(behav) #check freq

#randomly choose 3 segments of each behavior to be 'pure' instead of mixed
behav1.pure<- sample(which(behav==1), 3, replace = FALSE)
behav2.pure<- sample(which(behav==2), 3, replace = FALSE)
behav3.pure<- sample(which(behav==3), 3, replace = FALSE)

pure<- c(behav1.pure, behav2.pure, behav3.pure) %>% sort()


#create vector of behaviors within each time segment (duration of 100 steps)
behav.full<- vector("list", length(behav))
for (i in 1:length(behav)) {
  if (i %in% pure) {
    behav.full[[i]]<- rep(behav[i], 100)
  }else if (behav[i] == 1) {
    behav.full[[i]]<- sample(c(1,2,3), 100, replace = TRUE, prob = c(0.8, 0.1, 0.1))
  } else if (behav[i] == 2) {
    behav.full[[i]]<- sample(c(1,2,3), 100, replace = TRUE, prob = c(0.1, 0.8, 0.1))
  } else if (behav[i] == 3) {
    behav.full[[i]]<- sample(c(1,2,3), 100, replace = TRUE, prob = c(0.1, 0.1, 0.8))
  }
}
behav.full<- unlist(behav.full)

SL.params<- data.frame(shape=c(0.25, 2, 10), scale = c(1, 1, 1))
TA.params<- data.frame(mu=c(pi, pi, 0), rho = c(0.8, 0, 0.8))


#simulate track
#n=duration of each observation (behav.full), behav is a vector of behaviors, SL.params and TA.params are DFs of the necessary params from which to generate distributions for SL and TA from gamma and wrapped cauchy distribs, and Z0 is the initial location

track<- CRW.sim(n=1, behav = behav.full, SL.params = SL.params, TA.params = TA.params, Z0=c(0,0))
names(track)[5]<- "behav_fine"
track$behav_fine<- factor(track$behav)
levels(track$behav_fine)<- c("Resting","Exploratory","Transit")
track<- track %>% mutate(behav_coarse = c(NA, rep(behav, each=100)) %>% factor())
levels(track$behav_coarse)<- c("Resting","Exploratory","Transit")

true.brkpts<- which(diff(behav) != 0) * 100


# plot that puppy
ggplot(data = track[-1,], aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_coarse), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))


#compare distributions of SL and TA among fine-scale behaviors
ggplot(track, aes(SL)) +
  geom_density(aes(fill=behav_fine), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()

ggplot(track, aes(TA)) +
  geom_density(aes(fill=behav_fine), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()


#compare distributions of SL and TA among coarse-scale behaviors
ggplot(track, aes(SL)) +
  geom_density(aes(fill=behav_coarse), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()

ggplot(track, aes(TA)) +
  geom_density(aes(fill=behav_coarse), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()




################################
#### Run Segmentation Model ####
################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_segmentation_behavior")


library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')


#load and manipulate data
dat<- track
dat$dt<- 3600
dat$id<- 1

dat.list<- df.to.list(dat=dat)
names(dat)[3:4]<- c("dist","rel.angle")
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval

## Run RJMCMC
ngibbs = 40000
dat.res<- behavior_segment(dat = behav.list, ngibbs = ngibbs)
#takes 17 min for 40000 iterations


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, brk.cols = 99)  #brk.cols is max matrix cols

model.brkpts<- na.omit(t(brkpts[-1])) %>% as.vector()


## Heatmaps
plot.heatmap(data = behav.list, nbins = c(6,8), brkpts = brkpts, dat.res = dat.res, type = "behav")


## Compare True vs Modeled Breakpoints
all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts), type = rep(c("True","Model"),
                                                                        c(length(true.brkpts),
                                                                          length(model.brkpts))))

ggplot(all.brkpts, aes(x=brks, y=type)) +
  geom_point(size=2) +
  theme_bw() +
  labs(x="Time", y="Type") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10))

dat_out<- map(behav.list, assign.time.seg) %>% map_dfr(`[`)  #assign time seg and make as DF






#######################
#### Run LDA Model ####
#######################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_LDA_behavior")


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


#get data
dat.list<- df.to.list(dat_out)  #for later behavioral assignment
obs<- get.summary.stats_behav(dat_out)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
ind1=grep('y1',colnames(obs))
ind2=grep('y2',colnames(obs))
nmaxclust=max(length(ind1),length(ind2))-1  #max possible is 1 fewer than largest number of bins

#run Gibbs sampler
res=LDA_behavior_gibbs(dat=obs, gamma1=gamma1, alpha=alpha,
                       ngibbs=ngibbs, nmaxclust=nmaxclust,
                       nburn=nburn)

#Check traceplot of log marginal likelihood
plot(res$loglikel, type='l')

#Extract and plot proportions of behaviors per time segment
theta.post<- res$theta[(nburn+1):ngibbs,]  #extract samples from posterior
theta.estim<- theta.post %>% apply(2, mean) %>% matrix(nrow(obs), nmaxclust) #calc mean of posterior
boxplot(theta.estim, xlab="Behavior", ylab="Proportion of Total Behavior")

#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
round(apply(theta.estim, 2, sum)/nrow(theta.estim), digits = 3)

#true proportions of behavioral states
table(behav)/50


## Viz histograms from model
behav.res<- get_behav_hist(res)
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors

#Plot histograms of proportion data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3)[c(2,1,3)], guide = F) +
  facet_grid(param ~ behav, scales = "fixed")



## Viz behavior over time
#Assign behaviors (via theta) to each time segment
theta.estim<- apply(theta.estim[,1:3], 1, function(x) x/sum(x)) %>% t()  #normalize probs for only first 3 behaviors being used
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)
names(theta.estim)<- c("id", "tseg","Exploratory","Resting","Transit")  #define behaviors
nobs<- data.frame(id = obs$id, tseg = obs$tseg, n = apply(obs[,11:16], 1, sum)) #calc obs per tseg using SL bins (more reliable than TA)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim2<- aug_behav_df(dat = dat_out[-1,] %>% mutate(date=1:nrow(dat_out[-1,])),
                            theta.estim = theta.estim, nobs = nobs)

#Change into long format
theta.estim.long<- theta.estim2 %>% gather(key, value, -id, -tseg, -time1, -date)
names(theta.estim.long)[5:6]<- c("behavior","prop")
theta.estim.long$behavior<- factor(theta.estim.long$behavior,
                                   levels = c("Resting","Exploratory","Transit"))

#generate long form of true behavior
true.behavior<- matrix(0, 5000, 3) %>% data.frame(., time1 = 1:5000)
names(true.behavior)[1:3]<- c("Resting","Exploratory","Transit")
tseg<- rep(1:50, each = 100)
tmp<- data.frame(behav = behav.full, tseg = tseg)
for(i in 1:length(behav)) {
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


dat2<- assign_behav(dat.list = dat.list, theta.estim2 = theta.estim2)
dat2$behav<- factor(dat2$behav, levels = c("Resting","Exploratory","Transit"))

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




##############################
### Overall Model Accuracy ###
##############################

true.b.coarse<- as.numeric(dat$behav_coarse[-1])
model.b<- as.numeric(dat2$behav[-1])

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)
# 95.8% accuracy when including all different behaviors together at coarse scale


## For 'Resting' behavior
true.b.coarse_rest<- which(true.b.coarse == 1)
model.b_rest<- which(model.b == 1)
(which(true.b.coarse_rest %in% model.b_rest) %>% length()) / length(true.b.coarse_rest)
# 98.7% accuracy for 'Resting' at coarse scale

## For 'Exploratory' behavior
true.b.coarse_exp<- which(true.b.coarse == 2)
model.b_exp<- which(model.b == 2)
(which(true.b.coarse_exp %in% model.b_exp) %>% length()) / length(true.b.coarse_exp)
# 89.7% accuracy for 'Exploratory' at coarse scale

## For 'Transit' behavior
true.b.coarse_transit<- which(true.b.coarse == 3)
model.b_transit<- which(model.b == 3)
(which(true.b.coarse_transit %in% model.b_transit) %>% length()) / length(true.b.coarse_transit)
# 98.4% accuracy for 'Transit' at coarse scale

