setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(mapdata)
library(sf)
#library(gisland)

world_sf <- st_as_sf(map("worldHires", fill = T , plot = F))

png("World_sites_sf.png", width = 17, height = 15, units = "cm", res = 300)

world_sf %>% 
  ggplot() +
  geom_sf(lwd = 0, color = NA, fill = "grey60")+
  coord_sf(expand = F)+
  # Nothern Hemisphere
  
  # Fram Strait
  geom_point(aes(x = -7, y = 78.9995), shape = 23, fill = "#0072B2", size = 2.25, stroke =0.75)+
  # L4 Channel
  geom_point(aes(x = -4.2167, y = 50.25), shape = 23, fill = "#D55E00", size = 2.25, stroke =0.75)+
  # Bedford Basin
  geom_point(aes(x = -63.64028, y = 44.69361), shape = 23, fill = "#F0E442", size = 2.25, stroke =0.75)+
  # BBMO
  geom_point(aes(x = 2.48, y = 41.4), shape = 23, fill = "#AA4499", size = 2.25, stroke =0.75)+
  # SPOTS
  geom_point(aes(x = -118.24, y = 33.33), shape = 23, fill = "#CC6677", size = 2.25, stroke =0.75)+
  
  # Southern Hemisphere
  
  # Yongala
  geom_point(aes(x = 147.618 , y = -19.3085), shape = 23, fill = "#009E73", size = 2.25, stroke =0.75)+
  # Rottnest
  geom_point(aes(x = 115.417 , y = -32), shape = 23, fill = "#88CCEE", size = 2.25, stroke =0.75)+
  # Maria Island
  geom_point(aes(x = 148.233 , y = -42.5967), shape = 23, fill = "#882255", size = 2.25, stroke =0.75)+
  
  theme_bw(20)+
  xlab(expression('Longitude'))+
  ylab(expression("Latitude"))
dev.off()

###### Station specific plotting

# Fram Strait
FSXlim <- c(-40, 40)
FSYlim <- c(70, 84)

FS <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = FSXlim, ylim = FSYlim))

png("FramStrait.png", width = 14, height = 14, units = "cm", res = 300)
FS %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = FSXlim, ylim = FSYlim)+ 
  geom_point(aes(x = -7, y = 78.9995), shape = 23, fill = "#0072B2", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# L4 Channel
L4Xlim <- c(-10, 10)
L4Ylim <- c(48, 58)

L4 <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = L4Xlim, ylim = L4Ylim))

png("L4EnglishChannel.png", width = 14, height = 14, units = "cm", res = 300)
L4 %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = L4Xlim, ylim = L4Ylim)+ 
  geom_point(aes(x = -4.2167, y = 50.25), shape = 23, fill = "#D55E00", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# Bedford Basin
BBXlim <- c(-67, -59)
BBYlim <- c(42.5, 48)

BB <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = BBXlim, ylim = BBYlim))

png("BedfordBasin.png", width = 14, height = 14, units = "cm", res = 300)
BB %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = BBXlim, ylim = BBYlim)+ 
  geom_point(aes(x = -63.64028, y = 44.69361), shape = 23, fill = "#F0E442", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# BBMO
BBMOXlim <- c(-2, 7)
BBMOYlim <- c(38, 44.5)

BBMO <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = BBMOXlim, ylim = BBMOYlim))

png("BBMO.png", width = 14, height = 14, units = "cm", res = 300)
BBMO %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = BBMOXlim, ylim = BBMOYlim)+ 
  geom_point(aes(x = 2.48, y = 41.4), shape = 23, fill = "#AA4499", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# SPOT
SPOTXlim <- c(-122, -113)
SPOTYlim <- c(31, 37)

SPOT <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = SPOTXlim, ylim = SPOTYlim))

png("SPOT.png", width = 14, height = 14, units = "cm", res = 300)
SPOT %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = SPOTXlim, ylim = SPOTYlim)+ 
  geom_point(aes(x = -118.24, y = 33.33), shape = 23, fill = "#CC6677", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# Yongala
YongXlim <- c(135, 155)
YongYlim <- c(-25, -10)

Yong <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = YongXlim, ylim = YongYlim))

png("Yong.png", width = 14, height = 14, units = "cm", res = 300)
Yong %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = YongXlim, ylim = YongYlim)+ 
  geom_point(aes(x = 147.618 , y = -19.3085), shape = 23, fill = "#009E73", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# Rottnest
RottXlim <- c(110, 130)
RottYlim <- c(-35, -20)

Rott <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = RottXlim, ylim = RottYlim))

png("Rott.png", width = 14, height = 14, units = "cm", res = 300)
Rott %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = RottXlim, ylim = RottYlim)+ 
  geom_point(aes(x = 115.417 , y = -32), shape = 23, fill = "#88CCEE", size = 7, stroke =1.1)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()

# Maria
MariaXlim <- c(144, 150)
MariaYlim <- c(-44, -40)

Maria <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = MariaXlim, ylim = MariaYlim))

png("Maria.png", width = 14, height = 14, units = "cm", res = 300)
Maria %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = MariaXlim, ylim = MariaYlim)+ 
  geom_point(aes(x = 148.233 , y = -42.5967), shape = 23, fill = "#882255", size = 4, stroke =1.25)+
  theme_bw(20)+
  xlab("")+
  ylab("")
dev.off()
######### plain ggplot example if sf package fails for some reason, I don't like this as much. It makes axis units tricky

worldx <- c(-180, 180)
worldy <- c(-90, 90)

world <- ggplot2::map_data("worldHires", xlim = worldx, ylim = worldy)

png("World_sites.png", width = 17, height = 15, units = "cm", res = 300)

ggplot(data = world, aes(x = long, y = lat, group = group)) +
  geom_polygon( fill = "grey50") +
  coord_quickmap(xlim = worldx, ylim = worldy, expand = FALSE) +
  
  # Nothern Hemisphere
  
  # Fram Strait
  geom_point(aes(x = -7, y = 78.9995), shape = 23, fill = "#0072B2", size = 2.1, stroke =0.85)+
  # L4 Channel
  geom_point(aes(x = -4.2167, y = 50.25), shape = 23, fill = "#D55E00", size = 2.1, stroke =0.85)+
  # Bedford Basin
  geom_point(aes(x = -63.64028, y = 44.69361), shape = 23, fill = "#F0E442", size = 2.1, stroke =0.85)+
  # BBMO
  geom_point(aes(x = -2.48, y = 41.4), shape = 23, fill = "#661100", size = 2.1, stroke =0.85)+
  # SPOTS
  geom_point(aes(x = -118.24, y = 33.33), shape = 23, fill = "#CC6677", size = 2.1, stroke =0.85)+
  
  # Southern Hemisphere
  
  # Yongala
  geom_point(aes(x = 147.618 , y = -19.3085), shape = 23, fill = "#009E73", size = 2.1, stroke =0.85)+
  # Rottnest
  geom_point(aes(x = 115.417 , y = -32), shape = 23, fill = "#88CCEE", size = 2.1, stroke =0.85)+
  # Maria Island
  geom_point(aes(x = 148.233 , y = -42.5967), shape = 23, fill = "#882255", size = 2.1, stroke =0.85)+
  
  theme_bw(20)+
  xlab(expression('Longitude'))+
  ylab(expression("Latitude"))

dev.off()

  