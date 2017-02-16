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
library(pbapply)



                          #################################################
                          #                                               #
                          #   CLEANING POLYGONS AND DISTANCE CALCULATIONS #
                          #       EFFECTIVE TREATMENTS AND CONTROLS       #
                          #                                               #
                          #################################################

######################################
#                                    #
# 1.  GETTING THE DATA READY         #
#                                    #
###################################### 

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
#plot(res)

# Create a border colombia geometry (to clean shoreline)
# Division by departaments
colombia <- readOGR(dsn = "Colombia_Deptos", layer = "Colombia_Deptos") %>%
  spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  spTransform(CRS=CRS("+init=epsg:3857")) %>%
  .[!(.@data$DPTO_CCDGO == 88) , ] %>%
  gUnaryUnion(.) %>%
  hole_free()
# plot(deptos_ln)

# Polígono de la Línea Negra
linea_negra <- readOGR(dsn = "Linea_Negra", layer = "Linea_Negra_Polygon") %>%
  spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

linea_negra_proj <- linea_negra %>% spTransform(CRS=CRS("+init=epsg:3857"))


# Parques nacionales naturales
pnn <- readOGR(dsn = "WDPA/Protected_Areas/", layer = "WDPA_Jan2017_COL-shapefile-polygons") %>% 
  spTransform(CRS=CRS("+init=epsg:3857")) %>%
  .[!.@data$STATUS_YR > 2012 & !.@data$GIS_AREA < 1 & .@data$MANG_AUTH == "Parques Nacionales Naturales de Colombia" , ] %>%
  gIntersection(., colombia)

pnn_ln <- gIntersection(pnn, linea_negra_proj)

# Resguardos indigenas
resguardos <- readOGR(dsn = "IGAC", layer = "Resguardos_Selected_LNegra") %>%
  spTransform(linea_negra, CRS=CRS("+init=epsg:3857"))
resguardos@data$AREA_KM2 <- gArea(resguardos, byid = T)/1e6
resguardos <- resguardos %>% .[!.@data$AREA_KM2 < 1, ]

resguardos_ln <- gIntersection(resguardos, linea_negra_proj)

# Densify lines using geosphere package
lineanegra_p <- linea_negra %>% as("SpatialLines") %>% as("SpatialPoints")  # _p for points
lineanegra_dens <- makePoly(linea_negra, interval=1000, sp = T)  %>% # The higuer the interval the lower the number of new points
  spTransform(CRS=CRS("+init=epsg:3857")) 

#Plot the densified geometry
# lineanegra_dens_p <- lineanegra_dens %>% as("SpatialLines") %>% as("SpatialPoints")
# plot(linea_negra)
# plot(lineanegra_p, add = T)
# plot(lineanegra_dens_p, add = T, col = "red")


#Buffer linea negra
territories <-  list(linea_negra_proj, resguardos_ln, pnn_ln) %>%
  lapply(., function(x){
    spTransform(x, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  })

territories_buffer <- territories %>%
  spTransform(CRS=CRS("+init=epsg:3857"))
  lapply(., gBuffer, width = 50000, byid = T)


########################################
#                                      #
# 2.   CLEANING POLYGONS               #
#                                      #
########################################


#Prepare data (union of poly and points)
resguardos_p <- resguardos %>% as("SpatialLines") %>% as("SpatialPoints")
pnn_ln_p <- pnn_ln %>% as("SpatialLines") %>% as("SpatialPoints")

resguardos_hole_free <- hole_free(resguardos) %>% gBuffer(0.0001, byid = F)
pnn_ln_hole_free <- hole_free(pnn_ln) %>% gBuffer(0.0001, byid = F)

resg_pnn_p <- rbind(resguardos_p, pnn_ln_p)
resg_pnn <- raster::union(resguardos_hole_free, pnn_ln_hole_free) # raster is the package and union is the name of the specific function within the package. :: helps to access the exact function from that specific package
colombia_p <- colombia %>% as("SpatialLines") %>% as("SpatialPoints")
#Eliminate close to water (from colombia_raster.R)


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
clean <- clean_treatments(x = lineanegra_dens, polygon = resg_pnn, points_sp = resguardos_p, 
                          points_border = colombia_p)

#It worked?
plot(linea_negra_proj)
plot(clean, add = T, col = "red", pch = 19)

#Clean shoreline of resguardos and PNN

#Clean SpatialPoints (from polygons of Natural parks) -remove points ONLY near to national frontiers and border points-
clean_treatments_border <- function(x, points_border){
  sp <- x %>% as("SpatialLines") %>% as("SpatialPoints")
  knn_border <- get.knnx(coordinates(points_border), coordinates(sp), k = 1, algorithm = "kd_tree") %>%
    data.frame(.)
  sp_final <- SpatialPointsDataFrame(sp, knn_border, proj4string = CRS("+init=epsg:3857")) %>%
    .[!.@data$nn.dist < 5000, ] 
  
}

list_poly <- list(pnn_ln, resguardos_ln)
list_polygons_clean_border <- lapply(list_poly, clean_treatments_border, points_border = colombia_p)


#Resguardos witout shore
plot(list_poly[[1]])
plot(list_polygons_clean_border[[1]], add = T, col = "red", pch = 20)


#PNN witout shore
plot(list_poly[[2]])
plot(list_polygons_clean_border[[2]], add = T, col = "red", pch = 20)


#################################################
#                                               #
# 3.          CALCULATE DISTANCES               #
#                                               #
#################################################


#Project all geometries to WGS84 (easier than to reproject raster)
clean_frontiers <- c(clean,list_polygons_clean_border) %>%
  lapply(., function(x){
  if(typeof(x) == "S4" & length(x) > 0){
    sp <- spTransform(x, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    return(sp)
  } else { 
    return(0)
  }
})

territories_buffer <- territories_buffer %>%
  lapply(., function(x){
    if(typeof(x) == "S4" & length(x) > 0){
      sp <- spTransform(x, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      return(sp)
    } else { 
      return(0)
    }
  })


#Function to crop/mask raster and calculate distance using all cores from CPU.
calculate_distances_parallel <- function(buffer, points){
  if(length(points)  > 2 & typeof(points) == "S4"){
    crop(res[[1]], buffer) %>%
      mask(buffer) %>%
      clusterR(.,distanceFromPoints, args = list(xy = points)) %>%
      mask(buffer) %>%
      resample(res[[1]])
  } else
    return(0)
}

beginCluster()
system.time(mask_distance <- mapply(calculate_distances_parallel,
                                    buffer = territories_buffer, 
                                    points = clean_frontiers))
endCluster()


#################################################
#                                               #
# 4.       IDENTIFY CONTROL/TREATMENT           #
#                     CELLS                     #
#                                               #
#################################################


#Identify cells inside national parks and buffers and their identifier


cells_territories <- lapply(territories, function(x){
  cellFromPolygon(res[[1]], x)
})


#################################################
#                                               #
# 5.          EXTRACT DISTANCE DATA FROM        #
#                       RASTERS                 #
#                                               #
#################################################

#1. Extract distance as data frame per buffer (list element)
list_dataframes <- pblapply(mask_distance, as.data.frame, xy = T, na.rm = T)
names(list_dataframes) <- c(1, 2, 3)

#2. Extract row names (id cells) and define treatments
list_dataframes <- lapply(list_dataframes, function(x){
  x$ID <- row.names(x); x
})

# FAIL - Create treatments for all lists
# list_dataframes <- mapply(function(x, y){
#   x %>%
#   mutate(., treatment = ifelse(ID %in% unlist(y), 1, 0)) %>%
#     data.frame()
# }, x = list_dataframes, y = cells_territories, SIMPLIFY = T)
# 
# 

list_dataframes[[1]]$treatment <- ifelse(list_dataframes[[1]]$ID %in% unlist(cells_territories[[1]]), 1, 0)
list_dataframes[[2]]$treatment <- ifelse(list_dataframes[[2]]$ID %in% unlist(cells_territories[[2]]), 1, 0)
list_dataframes[[3]]$treatment <- ifelse(list_dataframes[[3]]$ID %in% unlist(cells_territories[[3]]), 1, 0)

#3. Append all elements of the list (1: LN, 2: Resg, 3: PNN) 
distance_dataframe <- do.call(rbind, list_dataframes)
distance_dataframe$buffer_id <- rep(names(list_dataframes), sapply(list_dataframes, nrow)) #identify cells from buffers


#4. Export
write.csv(distance_dataframe, "distance_dataframe_territories", row.names = FALSE)


  # Saving history Commands' Backups
  savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")











