#REFERENCE: https://github.com/ivanhigueram/deforestacion/blob/master/Cleaning/polygon_cleaning.R
library(raster)
library(rgdal)
library(rasterVis)
library(sp)
library(ggplot2)
library(maptools)
library(lattice)
library(plyr)
library(grid)
library(FNN)
library(dplyr)
library(broom)
library(stringr)
library(geosphere)
library(gtools)
library(magrittr)
library(rgeos)


                          #################################################
                          #               CLEANING POLYGONS               #
                          #       EFFECTIVE TREATMENTS AND CONTROLS       #
                         #      (CLEANING AND DISTANCE CALCULATIONS)      #
                          #################################################

######################################
# 1.  GETTING READY THE DATA
# #################################### 

#Function to clean holes from polygons
hole_free <- function(x){
  BCp <- slot(x, "polygons")
  holes <- lapply(BCp, function(y) sapply(slot(y, "Polygons"), slot, "hole")) 
  res <- lapply(1:length(BCp), function(i) slot(BCp[[i]], "Polygons")[!holes[[i]]]) 
  IDs <- row.names(x)
  SpatialPolygons(lapply(1:length(res), function(i) Polygons(res[[i]], ID=IDs[i])), proj4string = CRS(proj4string(x))) 
}

setwd("~/Dropbox/Linea_Negra_R/Data/")

# Loss-Year Brick Linea Negra 1km
res <- brick("~/Downloads/loss_year_brick_1km.tif")
#plot(res_ln)

# Resguardos indigenas
resguardos <- readOGR(dsn = "IGAC", layer = "Resguardos_Selected_LNegra") %>%
  spTransform(linea_negra, CRS=CRS("+init=epsg:3857"))
#plot(resguardos)

# Parques nacionales naturales
pnn_ln <- readOGR(dsn = "SIAC_ANLA", layer = "PNN_Selected_LNegra") %>%
  spTransform(CRS=CRS("+init=epsg:3857"))
# plot(pnn_ln)

# Polígono de la Línea Negra
linea_negra <- readOGR(dsn = "Linea_Negra", layer = "Linea_Negra_Polygon") %>%
  spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

linea_negra_proj <- linea_negra %>% spTransform(CRS=CRS("+init=epsg:3857"))

#Buffer linea negra
buffer_ln <- gBuffer(linea_negra_proj, width = 50000, byid = T)

# spTransform(CRS=CRS("+init=epsg:3857"))
# plot(linea_negra)

# Create a border colombia geometry (to clean shoreline)
  # Division buy departaments
    deptos_ln <- readOGR(dsn = "Colombia_Deptos", layer = "Colombia_Deptos") %>%
      spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
      spTransform(CRS=CRS("+init=epsg:3857"))
    # plot(deptos_ln)

# Reference: https://github.com/ivanhigueram/deforestacion/blob/master/Cleaning/colombia_raster.R
#Remove islands
#Remove municipalities that are out of continental land (Malpelo and Providencia) - No projection defined beforehand

deptos_ln <- deptos_ln[!(deptos_ln@data$DPTO_CCDGO == 88) , ] %>%
  gUnaryUnion(.)

# Remove holes
deptos_ln <- hole_free(x = deptos_ln)

# Proyect to meters
deptos_ln_proj <- spTransform(deptos_ln, CRS=CRS("+init=epsg:3857"))
deptos_ln <- list(deptos_ln, deptos_ln_proj)

# Densify lines using geosphere package
lineanegra_p <- linea_negra %>% as("SpatialLines") %>% as("SpatialPoints")
lineanegra_dens <- makePoly(linea_negra, interval=1000, sp = T)  %>% # The higuer the interval the lower the number of new points
  spTransform(CRS=CRS("+init=epsg:3857")) 

# lineanegra_dens_p <- lineanegra_dens %>% as("SpatialLines") %>% as("SpatialPoints")
# plot(linea_negra)
# plot(lineanegra_p, add = T)
# plot(lineanegra_dens_p, add = T, col = "red")

#Prepare data
resguardos_p <- resguardos %>% as("SpatialLines") %>% as("SpatialPoints")
resguardos_hole_free <- hole_free(resguardos)
resguardos_d <- unionSpatialPolygons(resguardos_hole_free, c(1:length(resguardos_hole_free@polygons))) %>%
  gUnaryUnion()

#Eliminate close to water (from colombia_raster.R)
deptos_ln_p <- deptos_ln[[2]] %>% as("SpatialLines") %>% as("SpatialPoints")

#Clean SpatialPoints (from polygons of Natural parks) -remove other treatments and get effective boundaries-
clean_treatments <- function(x, polygon, points_sp, points_border, shape){
  # print(x$ID)
  if(gIntersects(x, polygon)){
    #Remove inside points
    dif <- gDifference(x, polygon, drop_lower_td = T)
    if(!is.null(dif)){
      dif <- tidy(dif)[, 1:2] #Coordinates difference
      polygon2_coords <- tidy(x)[,1:2] #Coordinates polygon
      # Duplicated_coords is the non-intersecting points of the polygon2
      duplicated_coords <- merge(dif, polygon2_coords) 
      if(dim(duplicated_coords)[1] > 0){
        res <- SpatialPoints(duplicated_coords, proj4string = CRS("+init=epsg:3857"))
      } else {
        res <- SpatialPoints(polygon2_coords, proj4string = CRS("+init=epsg:3857"))
      }
      
    } else {
      return(0)
    }
    #Remove close cofounding treatments
    knn <- get.knnx(coordinates(points_sp), coordinates(res), k = 1, algorithm = "kd_tree") %>%
      data.frame(.)
    sp <- SpatialPointsDataFrame(res, knn, proj4string = CRS("+init=epsg:3857")) %>%
      .[!.@data$nn.dist < 1000, ]
    knn_border <- get.knnx(coordinates(points_border), coordinates(sp), k = 1, algorithm = "kd_tree") %>%
      data.frame(.)
    sp_border <- SpatialPointsDataFrame(sp, knn_border, proj4string = CRS("+init=epsg:3857")) %>%
      .[!.@data$nn.dist < 4000, ]
    # dif <- gDifference(shape, x) %>% as("SpatialLines") %>% as("SpatialPoints")
    # knn_final <- get.knnx(coordinates(dif), coordinates(sp_border), k = 1, algorithm = "kd_tree") %>%
    #   data.frame()
    # sp_final <- SpatialPointsDataFrame(sp_border, knn_final, proj4string = CRS("+init=epsg:3857")) %>%
    #   .[!.@data$nn.dist < 500, ]
    # 
  } else {
    # Remove close cofounding treatments
    points <- x %>% as("SpatialLines") %>% as("SpatialPoints")
    knn <- get.knnx(coordinates(points_sp), coordinates(points), k = 1, algorithm = "kd_tree") %>%
      data.frame(.)
    sp <- SpatialPointsDataFrame(points, knn, proj4string = CRS("+init=epsg:3857")) %>%
      .[!.@data$nn.dist < 1000, ]
    knn_border <- get.knnx(coordinates(points_border), coordinates(sp), k = 1, algorithm = "kd_tree") %>%
      data.frame(.)
    sp_border <- SpatialPointsDataFrame(sp, knn_border, proj4string = CRS("+init=epsg:3857")) %>%
      .[!.@data$nn.dist < 4000, ]
    # dif <- gDifference(shape, x) %>% as("SpatialLines") %>% as("SpatialPoints")
    # knn_final <- get.knnx(coordinates(dif), coordinates(sp_border), k = 1, algorithm = "kd_tree") %>%
    #   data.frame()
    # sp_final <- SpatialPointsDataFrame(sp_border, knn_final, proj4string = CRS("+init=epsg:3857")) %>%
    #   .[!.@data$nn.dist < 500, ]
  }
  
}

#Clean polygon 
clean <- clean_treatments(x = lineanegra_dens, polygon = resguardos_hole_free, points_sp = resguardos_p, 
                          points_border = deptos_ln_p)


#Calculate distance



  # Clear all plots
    dev.off() 

  # Saving history Commands' Backups
  savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")











