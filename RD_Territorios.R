
library(plyr)
library(dplyr)
library(raster)
library(sp)
library(rgdal)
library(rdrobust)
library(Hmisc)
library(data.table)


#Load data from rasters
setwd("~/Dropbox/Linea_Negra_R/Data/Dataframes/")
distance <- fread("distance_dataframe_protected_areas.csv") %>% plyr::rename(., c("layer" = "dist"))
defo <- fread("deforestation_dataframe.csv") %>% as.data.frame()
lights <- fread("rasters_lights_dataframe.csv")

#Aggregate deforestation (2001 - 2012)
defo$loss_sum <- rowSums(defo[, names(defo)[4:16]])
loss_sum <- dplyr::select(defo, c(ID, loss_sum))

#Merge data 
defo_dist <- distance %>%
  split(., .$buffer_id) %>% lapply(., function(x){
      merge(x, loss_sum, by = "ID") %>%
      merge(., lights, by = "ID") %>%
      mutate(., loss_sum = loss_sum * 100) %>%
      mutate(., dist_disc = ifelse(treatment == 1, 1, -1) * dist) %>%
      mutate(., dist_disc = dist_disc / 1000)
  })

# merge(., cov, by = "ID", all.x = T) %>%
# merge(., clump, by = "ID", all.x = T) %>%
# mutate(clumps = ifelse(is.na(clumps), 0, 1))

protect_areas[[5]] %>% spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
plot()
plot(list_polygons_clean_border[[5]], add = T, cex = 0.5, pch = 19, col = "red")

#RD - Deforestation
defo_dist_buffer <- 
  lapply(defo_dist, function(x){
    a <- rdrobust(x = x$dist_disc,
                  y = x$loss_sum,
                  c = 0,
                  vce = "nn",
                  all = T)
  })


#RD - Lights
lights_dist_buffer <- split(defo_dist, defo_dist$buffer_id) %>%
  lapply(., function(x){
    a <- rdrobust(x = x$dist_disc,
                  y = x$F182013.v4c_web.stable_lights.avg_vis,
                  c = 0,
                  vce = "nn",
                  all = T)
  })





#4. RD: Remove from RI all pixels that are not PNN

pnn <- defo_dist %>% filter(buffer_id == 3 & treatment == 1)
ri <- defo_dist %>% filter(buffer_id == 2 & treatment == 1)
ln <- defo_dist %>% filter(buffer_id == 1 & treatment == 1)

defo_dist_4 <- defo_dist %>% filter(buffer_id == 2) %>%
  filter(ID %in% pnn$ID) 


plot(territories[[1]])
plot(territories[[2]], add = T, border = "blue")
plot(territories[[3]], add = T, border = "red")
points(ri$x, ri$y)

# RD: Remove from PNN all pixels that are not RI []


# Saving history Commands' Backups
savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")



