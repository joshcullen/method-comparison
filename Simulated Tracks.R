### Create simulated data with activity centers and multiple behaviors as defined by differences in step length and turning angle


library(waddle)
library(tidyverse)
library(viridis)

### My Simulation: Only for identifying ACs since not enough control over SL and TA ###
## Uses multistate BCRW model
set.seed(1)

#multistate BCRW w 5 phases
ns <- c(100,50,100,50,100)
Z.centers <- c(0,10-5i,10-5i,30,30)
attractions <- c(0.9,0.9,0.9,0.9,0.9)
rhos <- c(0.2,0.2,0.2,0.9,0.2)
# a <- c(5,1,5,1,5)
# b <- c(1,2,1,10,1)
BCRW.sim <- multiBCRW(rhos=rhos, attractions=attractions, Z.centers=Z.centers, ns=ns)
plot(BCRW.sim, col = viridis(n=3)[c(1,2,1,3,1)])
points(Z.centers, cex=3, col="orange", lwd=2)

#output is list of vectors for coords, behaviors, and duration of each behavior

# step length
S <- diff(BCRW.sim$Z)

nu.means<- matrix(NA, 5, 1)
for (i in 1:length(unique(BCRW.sim$Phase))) {
  ind<- which(BCRW.sim$Phase == i)
  
  nu.means[i,]<- mean(BCRW.sim$Z[ind] %>% diff() %>% Mod(), na.rm = T)
}

plot.ts(Mod(S))
abline(v = c(100,150,250,300), lwd=3, col="darkgrey", lty=3)
clip(1,100,0,20)
abline(h=nu.means[1], col="red", lwd=2)
clip(101,150,0,20)
abline(h=nu.means[2], col="red", lwd=2)
clip(151,250,0,20)
abline(h=nu.means[3], col="red", lwd=2)
clip(251,300,0,20)
abline(h=nu.means[4], col="red", lwd=2)
clip(300,400,0,20)
abline(h=nu.means[5], col="red", lwd=2)

# absolute orientation
Phi <- Arg(diff(BCRW.sim$Z))
# turning angles
Theta <- diff(Phi)
# rose diagrams
require(circular)
par(mfrow=c(1,5))
rose.diag(Theta[1:100], bins=8, main="Phase I", prop=2, col=viridis(n=3)[1])
rose.diag(Theta[101:150],bins=8, main="Phase II", prop=2, col=viridis(n=3)[2])
rose.diag(Theta[151:250], bins=8, main="Phase III", prop=2, col=viridis(n=3)[1])
rose.diag(Theta[251:300], bins=8, main="Phase IV", prop=2, col=viridis(n=3)[3])
rose.diag(Theta[301:400], bins=8, main="Phase V", prop=2, col=viridis(n=3)[1])
par(mfrow=c(1,1))


### Testing for modification of shape and scale params of Weibull distribution
# set.seed(1)
# test<- BCRW(n = 100, rho=0.9, attraction = 0.9, a=1, b=30, Z.center = 20)
# dist<- sqrt((test$X)^2 + (test$Y)^2)
# S<- diff(dist)
# plot(test)
# plot.ts(S); mean(S)
# 
# 
# #Evaluate modifications of Weibull distrib (for SL)
# W.5_1<- rweibull(500, 5, 1) %>% data.frame(value = ., type = "a=5,b=1")  #Resting/Encamped
# W.1_2<- rweibull(500, 1, 2) %>% data.frame(value = ., type = "a=1,b=2")  #ARS/Exploratory
# W.1_10<- rweibull(500, 1, 10) %>% data.frame(value = ., type = "a=1,b=10")  #Transit
# 
# W<- rbind(W.5_1, W.1_2, W.1_10)
# ggplot(W, aes(x=value)) +
#   geom_density(aes(fill = type), alpha = 0.6)







### My Simulation: Only for identifying behaviors ###
## Uses multistate CVM model; nu = mean speed; tau = tortuosity or time-scale of autocorr
set.seed(10)

#define behaviors and randomly sample 50 (for 50 time segments)
#weight probs so that behavior 1 (Resting) occurs 50%, behavior 2 (ARS) occurs 35%, and behavior 3 (Transit) occurs 15%
behav<- sample(c(1,2,3), 50, replace = TRUE, prob = c(0.5, 0.35, 0.15))
nu.vals<- c(1,2.5,10)
nus <- c(1,2.5,1,10,1)
taus <- c(1,1,1,10,1)
Ts <- c(100,50,100,50,100)
CVM.sim <- multiCVM(taus, nus, Ts)
plot(CVM.sim, col = viridis(n=3)[c(1,2,1,3,1)])

#calculate and plot mean speeds
S <- diff(CVM.sim$Z) %>% Mod()

nu.means<- matrix(NA, 5, 1)
for (i in 1:length(unique(CVM.sim$Phase))) {
  ind<- which(CVM.sim$Phase == i)
  
  nu.means[i,]<- mean(CVM.sim$Z[ind] %>% diff() %>% Mod(), na.rm = T)
}


plot.ts(S)
abline(v = c(100-1,150-2,250-3,300-4), lwd=3, col="darkgrey", lty=3)
clip(1,99,0,20)
abline(h=nu.means[1], col="red", lwd=2)
clip(100,148,0,20)
abline(h=nu.means[2], col="red", lwd=2)
clip(149,247,0,20)
abline(h=nu.means[3], col="red", lwd=2)
clip(248,296,0,20)
abline(h=nu.means[4], col="red", lwd=2)
clip(297,395,0,20)
abline(h=nu.means[5], col="red", lwd=2)




# absolute orientation
Phi <- diff(CVM.sim$Z) %>% Arg()
# turning angles
Theta <- diff(Phi)
# rose diagrams
require(circular)
par(mfrow=c(1,5))
rose.diag(Theta[1:100], bins=8, main="Phase I", prop=2, col=viridis(n=3)[1])
rose.diag(Theta[101:150],bins=8, main="Phase II", prop=2, col=viridis(n=3)[2])
rose.diag(Theta[151:250], bins=8, main="Phase III", prop=2, col=viridis(n=3)[1])
rose.diag(Theta[251:300], bins=8, main="Phase IV", prop=2, col=viridis(n=3)[3])
rose.diag(Theta[301:400], bins=8, main="Phase V", prop=2, col=viridis(n=3)[1])
par(mfrow=c(1,1))




## New method for CVM

set.seed(3)


#define behaviors and randomly sample 50 (for 50 time segments)
#weight probs so that behavior 1 (Resting) occurs 50%, behavior 2 (ARS) occurs 35%, and behavior 3 (Transit) occurs 15%

#create vector of behaviors
behav<- sample(c(1,2,3), 50, replace = TRUE, prob = c(0.5, 0.35, 0.15))
table(behav) #check freq

#create vector of speeds per behavior
nus<- rep(0, length(behav))
for (i in 1:length(behav)) {
  if (behav[i] == 1) {
    nus[i]<- 1  #Rest
  } else if (behav[i] == 2) {
    nus[i]<- 2.5  #ARS
  } else {
    nus[i]<- 10  #Transit
  }
}


#create vector of tortuosity per behavior
taus<- rep(0, length(behav))
for (i in 1:length(behav)) {
  if (behav[i] < 3) {
    taus[i]<- 1  #Rest/ARS
  } else {
    taus[i]<- 10  #Transit
  }
}

#create vector of behavior duration; set to equal for all segments
Ts <- rep(50, length(behav))

#Run CVM model
CVM.sim <- multiCVM(taus, nus, Ts)
plot(CVM.sim, col = viridis(n=3)[behav])

#calculate and plot mean speeds
S <- diff(CVM.sim$Z) %>% Mod()
S<- c(S, NA)

plot.ts(S, ylab = "Step Length")



# Relative turning angle
Phi <- diff(CVM.sim$Z) %>% Arg()
# turning angles
Phi<- c(NA, Phi)

plot.ts(Phi, ylab = "Turning Angle (rad)")


behav.aug<- rep(behav, each=49)  #assign behavior for all obs; not sure why each seg is missing 1 obs

sim_behav<- data.frame(id = 1, dt = 3600, dist = S, rel.angle = Phi, time1 = 1:2450, behav = behav.aug)
write.csv(sim_behav, "Sim track_behavior.csv", row.names = FALSE)
