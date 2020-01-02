#Compare KDE against modeled ACs

library(tidyverse)
library(adehabitatHR)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ks)
library(raster)
library(sf)


dat<- read.csv("Snail Kite Gridded Data_AC.csv", as.is = TRUE)


###################
### Modeled ACs ###
###################

ac.coords<- read.csv("Activity Center Coordinates.csv", header = T, sep = ',')

## Map

#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- sf::st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N

# lakes
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
lakes10<- sf::st_transform(lakes10, crs = "+init=epsg:32617") %>%
  sf::st_crop(xmin = min(dat$x-20000), xmax = max(dat$x+20000), ymin = min(dat$y-20000),
              ymax = max(dat$y+20000))

nests<- dat %>% group_by(id) %>% dplyr::select(c(id, x, y)) %>% slice(n=1)


# ACs
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_point(data = dat, aes(x, y, fill = "Raw"), color = "grey55", shape = 21, size = 1,
             alpha = 0.2) +
  geom_point(data = nests, aes(x, y, fill = "Nests"), shape = 24, size = 2.5, alpha = 0.7) +
  geom_point(data = ac.coords, aes(x, y, fill = "ACs"), shape = 21, size = 3, alpha = 0.8) +
  labs(x="Longitude", y="Latitude") +
  scale_fill_manual("", values = c(viridis(n=3)[1],"red","grey55")) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,24,21)))) +
  theme_bw()




###########
### KDE ###
###########

## H plug-in
dat.coords<- dat[,c("x","y")]
Hpi.dat<- Hpi(x=dat.coords, nstage = 2)
kde.hpi<- kde(x=dat.coords, H=Hpi.dat, compute.cont = T)
kd_df <- expand.grid(x=kde.hpi$eval.points[[1]], y=kde.hpi$eval.points[[2]]) %>% 
  mutate(z = c(kde.hpi$estimate)/sum(kde.hpi$estimate))

ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "black", alpha = 0.5) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_tile(data=kd_df, aes(x, y, fill=z), alpha=0.75) +
  geom_point(data = ac.coords, aes(x, y), color = viridis(n=9)[7], size = 2, pch=1, stroke=0.5) +
  scale_fill_viridis_c("Proportional \nDensity \nEstimate") +
  theme_bw() +
  labs(x="Longitude", y="Latitude")



############################
### ACs vs KDE over Time ###
############################

dat2<- dat
dat2$date<- as.POSIXct(dat2$date, format = "%Y-%m-%d")
dat2<- dat2[order(dat2$date),]

# prop_ac<- matrix(0, length(unique(dat2$date)), 20) #empty matrix to fill proportions
# oo=1
# 
# for (i in unique(dat2$date)) {
#   ind<- which(dat2$date == i)
#   tmp<- prop.table(table(dat2[ind,]$ac))
#   
#   prop_ac[oo, as.numeric(names(tmp))]<- tmp
#   
#   oo=oo+1
# }
# 
# prop_ac<- cbind(date = unique(dat2$date), prop_ac) %>% data.frame()
# prop_ac$date<- as.POSIXct(prop_ac$date, origin = "1970-01-01 00:00:00")
# names(prop_ac)[2:21]<- as.character(1:20)
# prop_ac_long<- prop_ac %>% gather(key, value, -date)
# prop_ac_long$key<- factor(prop_ac_long$key, levels = 1:20)
# 
# 
# 
# ggplot() +
#   geom_area(data = prop_ac_long, aes(x=date, y=value, fill = key), color = "black", size = 0.25,
#             position = "fill") +
#   labs(x = "\nTime", y = "Proportion of Day Observed at AC\n") +
#   scale_fill_viridis_d("Behavior", direction = -1) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14),
#         panel.grid = element_blank())


ggplot() +
  geom_point(data = dat2, aes(x=date, y=ac)) +
  labs(x = "\nTime", y = "Activity Center\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        panel.grid = element_blank()) +
  scale_y_continuous(breaks = 1:20)
