## CRW model

library(tidyverse)
library(circular)

source('Simulation Functions.R')


#################
#### Resting ####
#################

# length of walk
N<-1000

#make gamma distributed steps
steps<- rgamma(N, shape = 0.25, scale = 1)

# check out their distribution
ggplot(data.frame(steps), aes(steps)) +
  geom_density()

summary(steps)


# make clustered turning angles
theta <- rwrappedcauchy(N, mu=circular(pi), rho=0.8)
theta <- ifelse(theta > pi, theta-(2*pi), theta)
# check out their distribution
rose.diag(theta,bins=24)
hist(as.numeric(theta))
ggplot(data.frame(theta), aes(theta)) +
  geom_density()

# cumulative angle (absolute orientation)
Phi <- cumsum(theta)

# step length components
dX <- steps*cos(Phi)
dY <- steps*sin(Phi)

# actual X-Y values
X<-cumsum(dX)
Y<-cumsum(dY)
track<- data.frame(X,Y)

# plot that puppy
ggplot(data = track, aes(X,Y)) +
  geom_path() +
  geom_point(data = track[1,], aes(X,Y), color = "green", pch = 16, size = 3) +
  geom_point(data = track[N,], aes(X,Y), color = "red", pch = 17, size = 3) +
  coord_equal() +
  theme_bw()






###################
#### Exploring ####
###################

# length of walk
N<-1000

#make gamma distributed steps
steps<- rgamma(N, shape = 2, scale = 1)

# check out their distribution
ggplot(data.frame(steps), aes(steps)) +
  geom_density()

summary(steps)


# make clustered turning angles
theta <- rwrappedcauchy(N, mu=circular(pi), rho=0.2)
theta <- ifelse(theta > pi, theta-(2*pi), theta)
# check out their distribution
rose.diag(theta,bins=24)
hist(as.numeric(theta))
ggplot(data.frame(theta), aes(theta)) +
  geom_density()

# cumulative angle (absolute orientation)
Phi <- cumsum(theta)

# step length components
dX <- steps*cos(Phi)
dY <- steps*sin(Phi)

# actual X-Y values
X<-cumsum(dX)
Y<-cumsum(dY)
track<- data.frame(X,Y)

# plot that puppy
ggplot(data = track, aes(X,Y)) +
  geom_path() +
  geom_point(data = track[1,], aes(X,Y), color = "green", pch = 16, size = 3) +
  geom_point(data = track[N,], aes(X,Y), color = "red", pch = 17, size = 3) +
  coord_equal() +
  theme_bw()





###################
#### Dispersal ####
###################

# length of walk
N<-1000

#make gamma distributed steps
steps<- rgamma(N, shape = 10, scale = 1)

# check out their distribution
ggplot(data.frame(steps), aes(steps)) +
  geom_density()

summary(steps)


# make clustered turning angles
theta <- rwrappedcauchy(N, mu=circular(0), rho=0.8)
theta <- ifelse(theta > pi, theta-(2*pi), theta)
# check out their distribution
rose.diag(theta,bins=24)
hist(as.numeric(theta))
ggplot(data.frame(theta), aes(theta)) +
  geom_density()

# cumulative angle (absolute orientation)
Phi <- cumsum(theta)

# step length components
dX <- steps*cos(Phi)
dY <- steps*sin(Phi)

# actual X-Y values
X<-cumsum(dX)
Y<-cumsum(dY)
track<- data.frame(X,Y)

# plot that puppy
ggplot(data = track, aes(X,Y)) +
  geom_path() +
  geom_point(data = track[1,], aes(X,Y), color = "green", pch = 16, size = 3) +
  geom_point(data = track[N,], aes(X,Y), color = "red", pch = 17, size = 3) +
  coord_equal() +
  theme_bw()


####################################
####################################
#### Hard Clustering Simulation ####
####################################
####################################

### Simulate full track ###

#define behaviors and randomly sample 50 (for 50 time segments)
#weight probs so that behavior 1 (Resting) occurs 50%, behavior 2 (Exploratory) occurs 35%, and behavior 3 (Transit) occurs 15%

set.seed(1)

#create vector of behaviors
behav<- sample(c(1,2,3), 50, replace = TRUE, prob = c(0.5, 0.35, 0.15))
table(behav) #check freq

n=50
SL.params<- data.frame(shape=c(0.25, 2, 10), scale = c(1, 1, 1))
TA.params<- data.frame(mu=c(pi, pi, 0), rho = c(0.8, 0, 0.8))


#simulate track
#n=duration of each time segment (behav), behav is a vector of behaviors, SL.params and TA.params are DFs of the necessary params from which to generate distributions for SL and TA from gamma and wrapped cauchy distribs, and Z0 is the initial location

track<- CRW.sim(n=n, behav = behav, SL.params = SL.params, TA.params = TA.params, Z0=c(0,0))
track$true.behav<- factor(track$true.behav)
levels(track$true.behav)<- c("Resting","Exploratory","Transit")
true.brkpts<- which(diff(behav) != 0) * n


# plot that puppy
ggplot(data = track[-1,], aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=true.behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))


#compare distributions of SL and TA among behaviors
ggplot(track, aes(SL)) +
  geom_density(aes(fill=true.behav), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()

ggplot(track, aes(TA)) +
  geom_density(aes(fill=true.behav), alpha = 0.6, na.rm = T) +
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
ngibbs = 10000
dat.res<- behavior_segment(dat = behav.list, ngibbs = ngibbs)
#takes 3.5 min


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
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(param ~ behav, scales = "fixed")




## Viz behavior over time
#Assign behaviors (via theta) to each time segment
theta.estim<- apply(theta.estim[,1:3], 1, function(x) x/sum(x)) %>% t()  #normalize probs for only first 3 behaviors being used
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)
names(theta.estim)<- c("id", "tseg", "Resting","Exploratory","Transit")  #define behaviors
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
true.behavior<- matrix(0, 2500, 3) %>% data.frame(., time1 = 1:2500)
names(true.behavior)[1:3]<- c("Resting","Exploratory","Transit")
ind<- rep(behav, each = n)
for(i in 1:nrow(true.behavior)) {
  true.behavior[i, ind[i]]<- 1
}
true.behavior.long<- true.behavior %>% gather(key, value, -time1)
names(true.behavior.long)[2:3]<- c("behavior","prop")
true.behavior.long$behavior<- factor(true.behavior.long$behavior,
                                   levels = c("Resting","Exploratory","Transit"))

ggplot(theta.estim.long) +
  geom_line(aes(x=time1, y=prop, color = behavior), size = 1) +
  geom_line(data = true.behavior.long, aes(x=time1, y=prop, color=behavior), size = 0.5) +
  labs(x = "\nObservation", y = "Proportion of Behavior\n") +
  scale_color_viridis_d("Behavior", guide = F) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12)) +
  facet_wrap(~behavior)


dat2<- assign_behav(dat.list = dat.list, theta.estim2 = theta.estim2)
dat2$behav<- factor(dat2$behav, levels = c("Resting","Exploratory","Transit"))

ggplot() +
  geom_path(data = dat2, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat2[-1,], aes(x, y, fill=behav), size=2.5, pch=21,
             alpha=dat2$prop[-1]) +
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

true.b<- ind
model.b<- as.numeric(dat2$behav[-nrow(dat2)])

(which(true.b == model.b) %>% length()) / length(true.b)
# 97.7% accuracy when including all different behaviors together


## For 'Resting' behavior
true.b_rest<- which(true.b == 1)
model.b_rest<- which(model.b == 1)
(which(true.b_rest %in% model.b_rest) %>% length()) / length(true.b_rest)
# 98.9% accuracy for 'Resting'


## For 'Exploratory' behavior
true.b_exp<- which(true.b == 2)
model.b_exp<- which(model.b == 2)
(which(true.b_exp %in% model.b_exp) %>% length()) / length(true.b_exp)
# 96.4% accuracy for 'Exploratory'


## For 'Transit' behavior
true.b_transit<- which(true.b == 3)
model.b_transit<- which(model.b == 3)
(which(true.b_transit %in% model.b_transit) %>% length()) / length(true.b_transit)
# 98.0% accuracy for 'Exploratory'


