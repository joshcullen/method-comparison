###########################
### Hard-Clustering Sim ###
###########################

### Assess model accuracy
dat<- read.csv("Modeled HC Sim Tracks w Behav.csv", as.is = T)

## Overall
true.b<- dat$true.behav[-1] %>% factor(levels = c("Resting","Exploratory","Transit")) %>%
  as.numeric()
model.b<- dat$behav[-1] %>% factor(levels = c("Resting","Exploratory","Transit")) %>%
  as.numeric()

(which(true.b == model.b) %>% length()) / length(true.b)
# 97.7% accuracy when including all different behaviors together


## For 'Resting' behavior
true.b_rest<- which(true.b == 1)
model.b_rest<- which(model.b == 1)
(which(true.b_rest %in% model.b_rest) %>% length()) / length(true.b_rest)
# 97.9% accuracy for 'Resting'


## For 'Exploratory' behavior
true.b_exp<- which(true.b == 2)
model.b_exp<- which(model.b == 2)
(which(true.b_exp %in% model.b_exp) %>% length()) / length(true.b_exp)
# 97.3% accuracy for 'Exploratory'


## For 'Transit' behavior
true.b_transit<- which(true.b == 3)
model.b_transit<- which(model.b == 3)
(which(true.b_transit %in% model.b_transit) %>% length()) / length(true.b_transit)
# 98.3% accuracy for 'Exploratory'





############################
### Mixed-Membership Sim ###
############################

### Assess model accuracy

dat<- read.csv("Modeled MM Sim Tracks w Behav.csv", as.is = T)

## Overall
true.b.coarse<- dat$behav_coarse[-1] %>% factor(levels = c("Resting","Exploratory","Transit")) %>%
  as.numeric()
model.b<- dat$behav[-1] %>% factor(levels = c("Resting","Exploratory","Transit")) %>%
  as.numeric()

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)
# 92.7% accuracy when including all different behaviors together at coarse scale


## For 'Resting' behavior
true.b.coarse_rest<- which(true.b.coarse == 1)
model.b_rest<- which(model.b == 1)
(which(true.b.coarse_rest %in% model.b_rest) %>% length()) / length(true.b.coarse_rest)
# 99.4% accuracy for 'Resting' at coarse scale

## For 'Exploratory' behavior
true.b.coarse_exp<- which(true.b.coarse == 2)
model.b_exp<- which(model.b == 2)
(which(true.b.coarse_exp %in% model.b_exp) %>% length()) / length(true.b.coarse_exp)
# 78.9% accuracy for 'Exploratory' at coarse scale

## For 'Transit' behavior
true.b.coarse_transit<- which(true.b.coarse == 3)
model.b_transit<- which(model.b == 3)
(which(true.b.coarse_transit %in% model.b_transit) %>% length()) / length(true.b.coarse_transit)
# 98.7% accuracy for 'Transit' at coarse scale

