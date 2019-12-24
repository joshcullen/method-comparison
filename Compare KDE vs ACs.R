#Compare KDE against modeled ACs

library(tidyverse)
library(adehabitatHR)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ks)
library(raster)
library(sf)


dat<- read.csv("Snail Kite Gridded Data_larger.csv", as.is = TRUE)


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
  geom_point(data = dat, aes(x, y, color = "Raw"), size = 1, alpha = 0.5) +
  geom_point(data = ac.coords, aes(x, y, color = "Model"), size = 3) +
  geom_point(data = nests, aes(x, y, color = "Nests"), shape = 17, size = 2) +
  labs(x="Longitude", y="Latitude") +
  scale_color_manual("", values = c(viridis(n=3)[1],"red","grey45")) +
  guides(color = guide_legend(override.aes = list(shape = c(16,17,16)))) +
  theme_bw()




###########
### KDE ###
###########

dat.spdf<- dat[,c("x","y")]
coordinates(dat.spdf)<- ~x+y
proj4string(dat.spdf)<- CRS("+init=epsg:32617")

# 5 km w 1 cell buffer
grid<- raster(extent(dat.spdf) + 10000)
res(grid)<- 500
proj4string(grid)<- CRS("+init=epsg:32617")
grid[]<- 0

grid.spxdf<- as(grid, 'SpatialPixels')

kde.href<- kernelUD(dat.spdf, h="href", kern = "bivnorm", grid = grid.spxdf)
image(kde.href)
points(ac.coords[,1:2], pch=16, cex=1)
