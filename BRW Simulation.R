library(tidyverse)
library(circular)


#######################################
#######################################
#### Biased Random Walk Simulation ####
#######################################
#######################################

### Simulate full track ###

#randomly sample coordinates for 50 ACs

#create vector of coordinates
set.seed(1)
AC.x<- sample(50, 20, replace = FALSE)
set.seed(3)
AC.y<- sample(50, 20, replace = FALSE)
AC<- data.frame(x=AC.x, y=AC.y)

ggplot(AC, aes(x, y)) +
  geom_point(size = 2) +
  theme_bw() +
  coord_equal()


#simulate tracks as stationary process around each AC


BRW.sim=function (n = 50, a, b, Z.center = c(0,0), Z0 = c(3,4), rho) {
  #based on BCRW model described in Bailey et al. 2018 "Navigational efficiency in a biased and correlated random walk model of individual animal movement". Ecology. 99(1): 217-223.
  
  Z <- matrix(c(Z0, rep(NA, n*2 - 2)), n, 2, byrow = T)  #matrix of simulated locations
  ref.pt<- c(1,0)  #for reference vector
  
  #function to calc angle between 2 vectors
  angle <- function(x,y){
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    as.numeric(theta)
  }
  
  for (i in 2:n) {
    u<- as.matrix(c(Z.center-Z[i-1,]))  #vector formed by location and AC
    v<- as.matrix(ref.pt)  #reference vector for positive x-axis (representing 0 rad)
    
    #angle between loc-AC vector and the positive x-axis multiplied by the sign of the vector (x-component)
    if (u[1] < 0 & u[2] < 0) {
      ang.sign<- -1  #1st quadrant
    } else if (u[1] > 0 & u[2] < 0) {
      ang.sign<- -1  #2nd quadrant
    } else if (u[1] > 0 & u[2] > 0) {
      ang.sign<- 1  #3rd quadrant
    } else {
      ang.sign<- 1  #4th quadrant
    }
    
    #angle directed towards AC from location
    omega<- angle(t(u), v) * ang.sign
    if (is.na(omega))
      omega<- 0
    
    #generate navigation error (-pi < theta < pi) so that mean of cos(theta) is in [0,1]
    phi <- rnorm(1, 0, 0.5)
      
    #angle towards AC + error
    mu <- omega + phi
    if (mu < 0)
      mu<- 2 * pi + (omega + phi)
    theta <- rwrappedcauchy(1, mu = circular(mu), rho = rho) %>% as.numeric()
    
    # new step
    dist<- rgamma(1, a, b)  #step length
    dx <- dist * cos(theta)
    dy <- dist * sin(theta)
    
    # actual X-Y values
    x <- Z[i-1, 1] + dx
    y <- Z[i-1, 2] + dy
    Z[i,] <- c(x,y)
  }
  track <- data.frame(x = Z[,1], y = Z[,2])
  
  return(track)
}


set.seed(1)
test<- BRW.sim(n=500, a=2, b=1, Z.center = c(0,0), Z0 = c(3,4), rho = 0.3)

Z.center<- data.frame(x=0, y=0)

ggplot(data = test, aes(x, y)) +
  geom_path(size=0.5, color='gray75') +
  geom_point(size = 2) +
  geom_point(data = test[1,], aes(x,y), color = "green", pch = 1, size = 3, stroke = 1.5) +
  geom_point(data = test[nrow(test),], aes(x,y), color = "red", pch = 2, size = 3, stroke = 1.5) +
  geom_point(data = Z.center, aes(x,y), color = "gold", size = 4) +
  theme_bw() +
  coord_equal()







##############################


multiBRW.sim=function (Z0, Z.centers, n, nphases, a, b, rho, ...) {
  
  Z.list <- list()
  Z.centers<- as.matrix(Z.centers)
  
  for (i in 1:nphases) {
    if (i == 1) {
      Z0 = Z0
    } else {
      Z0 = Z.list[[i - 1]][n,] %>% as.numeric()
    }
    
    #for visiting each Z.center (AC) once depending on nphases
    if (nphases > nrow(Z.centers) & nphases %% nrow(Z.centers) == 0) {
      
      fold<- ceiling(nphases / nrow(Z.centers))  #number of times to sample
      ind<- c(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = FALSE),
              rep(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = TRUE), fold-1))
      
    } else if (nphases > nrow(Z.centers) & nphases %% nrow(Z.centers) != 0) {
      
      fold<- ceiling(nphases / nrow(Z.centers))  #number of times to sample
      tmp<- nphases %% nrow(Z.centers)  #remainder
      ind<- c(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = FALSE),
              rep(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = TRUE), fold-2),
              sample(x = nrow(Z.centers), size = tmp, replace = TRUE))
      
    } else {
      
    ind<- sample(x = nrow(Z.centers), size = nphases, replace = FALSE)
    
    }
    
    Z.list[[i]] <- BRW.sim(n = n, a = a, b = b, 
                           Z.center = Z.centers[ind[i],], 
                           Z0 = Z0, rho = rho)
  }
  
  track<- do.call(rbind.data.frame, Z.list)
  track$true.ac<- rep(ind, each = n)
  return(track)
}

set.seed(2)
track<- multiBRW.sim(n=250, a = 1, b = 1, nphases = 30, Z.center = AC, Z0 = c(3,4), rho = 0.3)
track$time1<- 1:nrow(track)


ggplot(data = track, aes(x, y)) +
  geom_path(size=0.5, color='gray75') +
  geom_point(size = 2, alpha = 0.3) +
  geom_point(data = AC, aes(x,y), color = "gold", size = 4) +
  geom_point(data = track[1,], aes(x,y), color = "green", pch = 1, size = 3, stroke = 1.5) +
  geom_point(data = track[nrow(track),], aes(x,y), color = "red", pch = 2, size = 3, stroke = 1.5) +
  theme_bw() +
  coord_equal()










#########################################


### Data Prep for Segmentation Model ###

library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(raster)
library(rgdal)
library(adehabitatLT)




### Create Grid to Discretize Space

# 2 unit res w 2 unit buffer on each side
grid<- raster(extent(min(track$x), max(track$x), min(track$y), max(track$y)) + 4)
res(grid)<- 2
grid[]<- 0
track$grid.cell<- cellFromXY(grid, track[,c("x","y")])


### Plot all points over grid

#Create grid cell borders
borders<- rasterToPolygons(grid, dissolve = F)
borders_f<- fortify(borders)

#Calc points per cell
tab<- table(cellFromXY(grid, track[,c("x","y")]))
grid[as.numeric(names(tab))] <- tab
grid_f<- as.data.frame(grid, xy = TRUE)
names(grid_f)[3]<- "count"


#plot points over grid
ggplot() +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = track, aes(x=x, y=y), color = "firebrick", size=0.5, alpha=0.5) +
  labs(x = "X", y = "Y") +
  theme_bw() +
  coord_equal()

#plot density surface of points in grid
ggplot() +
  geom_tile(data=grid_f, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "X", y = "Y") +
  theme_bw()






################################
#### Run Segmentation Model ####
################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_segmentation_model")


library(tidyverse)
library(progress)
library(furrr)
library(tictoc)
library(viridis)

source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')


## Prepare data

track$id<- 1
dat<- track
dat.list<- df.to.list(dat = dat)

#only select necessary cols and re-number grid cell IDs
dat.long<- map_dfr(dat.list, `[`) %>% dplyr::select(id, grid.cell, time1)  #create DF
names(dat.long)[2]<- "loc.id"
dat.long$loc.id<- dat.long$loc.id %>% factor()
levels(dat.long$loc.id)<- 1:length(unique(dat.long$loc.id))  #change from raw to modified cell ID
dat.long$loc.id<- dat.long$loc.id %>% as.character() %>% as.numeric()

#convert back to list
dat.list2<- df.to.list(dat.long)



identity<- unique(dat.long$id)
ngibbs = 10000

## Run Gibbs sampler

dat.res<- gibbs.time.seg(k = dat.list2[[1]], identity = identity, ngibbs = ngibbs)
###Takes 30 min to run for 1000 iterations


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


## Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- matrix(c(1,dat.res$breakpt[[ML]]), 1, length(dat.res$breakpt[[ML]])+1) %>% data.frame()
names(brkpts)[1]<- "id"


## Heatmaps

#since trouble w/ line of code in plot.heatmap(), including here with modified code
data<- dat.list2[[1]]

#re-define loc.id based only on those visited by this individual
uni.loc=unique(data$loc.id)
aux=data.frame(loc.id=uni.loc,loc.id1=1:length(uni.loc))
dat1=merge(data,aux,all=T)
dat1$loc.id=dat1$loc.id1
dat=dat1[order(dat1$time1),c('loc.id','time1')]

nloc<- length(uni.loc)
nobs<- nrow(data)
obs<- matrix(0, nobs, nloc)

for (i in 1:nrow(data)) {
  obs[i, dat$loc.id[i]]<- 1
}

obs<- data.frame(obs)
names(obs)<- 1:nloc
obs.long<- obs %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs.long$key<- as.numeric(obs.long$key)
obs.long$value<- factor(obs.long$value)
levels(obs.long$value)<- c("Absence","Presence")

ind=which(unique(data$id) == brkpts$id)
breakpt<- brkpts[ind,-1] %>% t() %>% data.frame()
names(breakpt)<- "breaks"


ggplot(obs.long, aes(x=time, y=key, fill=value)) +
    geom_tile() +
    scale_fill_viridis_d("") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = breakpt, aes(xintercept = breaks), color = viridis(n=9)[7], size = 0.15) +
    labs(x = "Observations", y = "Grid Cell") +
    theme_bw() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
          title = element_text(size = 20))


#Compare true and modeled breakpoints

true.brkpts<- which(diff(track$true.ac) != 0)
model.brkpts<- t(breakpt) %>% as.vector()
all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts), type = rep(c("True","Model"),
                                                                        c(length(true.brkpts),
                                                                          length(model.brkpts))))
ggplot(all.brkpts, aes(x=brks, y=type)) +
  geom_point(size=2) +
  theme_bw()


######################################
#### Assign Spatial Time Segments ####
######################################

breakpt<- brkpts[,-1] %>% as.numeric(.[1,])
breakpt1=c(0,breakpt,Inf)
n=length(breakpt1)
res=matrix(NA,nrow(dat.list[[1]]),1)
for (i in 2:n){
  ind=which(dat.list2[[1]]$time1>=breakpt1[i-1] & dat.list2[[1]]$time1<breakpt1[i])
  res[ind,]=i-1
}
dat_out<- dat.list[[1]]
dat_out$tseg<- as.vector(res)




###################################
#### Identify Activity Centers ####
###################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/activcenter_subset_locations")

set.seed(10)

#load libraries and read important functions
library('Rcpp')
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(progress)
library(rnaturalearth)
library(rnaturalearthdata)


sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')
source('helper functions.R')



#load data
dat<- dat_out
obs<- get.summary.stats_obs(dat)  #frequency of visitation in each location (column) for each time segment (row)

#geographical coordinates of locations
grid<- raster(extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y)) + 4)
res(grid)<- 2
grid[]<- 0

grid.cell.locs<- coordinates(grid) %>% data.frame()
names(grid.cell.locs)<- c("x", "y")
grid.cell.locs$grid.cell<- 1:length(grid)
grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]

#Define initial activity centers (top 50 by # of obs)
tmp<- colSums(obs[,-1]) %>% data.frame(grid.cell = colnames(obs[,-1]), nobs = .) %>%
  arrange(desc(nobs)) %>% slice(n=1:50) %>% dplyr::select(grid.cell)
tmp<- tmp$grid.cell %>% as.character() %>% as.numeric()
ind<- sample(tmp, size = 20, replace = F)

ac.coord.init<- grid.coord[ind,]

#top 50
ac.coord.init2<- grid.coord[tmp,]

#potential locations for activity centers (AC)
possib.ac=grid.coord #these don't have to be identical (i.e., we can define AC's on a coarser grid)

### Run Gibbs sampler

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=20
gamma1=0.1

#run gibbs sampler
options(warn=2)

res=gibbs.activity.center(dat=obs[,-1],grid.coord=grid.coord[,-3],n.ac=n.ac,
                          ac.coord.init=ac.coord.init[,-3],gamma1=gamma1,
                          possib.ac=possib.ac[,-3])

#plot output and look at frequency of AC visitation
plot(res$logl,type='l')
plot(res$phi,type='l')



##############################################
### Extract AC Coordinates and Assignments ###
##############################################

ac<- res$z[which.max(res$logl),]  #use ACs from iteration with max log likelihood
ac.coords<- matrix(NA, length(unique(ac)), 2)
colnames(ac.coords)<- c("x","y")
tmp<- res$coord[which.max(res$logl),]

for (i in 1:length(unique(ac))) {
  ac.coords[i,]<- round(c(tmp[i], tmp[i+length(unique(ac))]), 0)
}

ac.coords<- data.frame(ac.coords, ac=1:length(unique(ac)))

table(ac)


############################
### Add ACs to Dataframe ###
############################

tseg.length<- dat %>% group_by(id, tseg) %>% tally()
tseg.length<- tseg.length$n
ac.aug<- rep(ac, times = tseg.length)

dat$ac<- ac.aug

#Calculate number of obs per AC
dat %>% group_by(id) %>% dplyr::select(ac) %>% table()



## Map

ggplot() +
  geom_point(data = dat, aes(x, y), color="grey45", size=3, alpha = 0.5) +
  geom_point(data = ac.coords, aes(y, x), color = "blue", size = 4, pch = 1, stroke = 1) +
  geom_point(data = AC, aes(x, y), color = "gold", size = 4, pch = 1, stroke = 1) +
  scale_color_viridis_c("Activity Center") +
  labs(x = "X", y = "Y") +
  theme_bw() +
  coord_equal()
