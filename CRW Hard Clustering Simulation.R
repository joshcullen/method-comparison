####################################
####################################
#### Hard Clustering Simulation ####
####################################
####################################


library(tidyverse)
library(circular)

source('Simulation Functions.R')

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





write.csv(track, "CRW_HC_sim.csv", row.names = F)


