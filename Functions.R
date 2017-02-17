################################################################################
#                     FUNCTIONS FOR GEOMETRIES PROCESSING                      #
###############################################################################


# FUNCIÓN PARA CREAR EL DATA SLOT

DataSlot <- function(poly){
  (poly.df <- data.frame(ID=1:length(poly))) 
  rownames(poly.df)
  
  # Extract polygon ID's
  (poly_id <- sapply(slot(poly, "polygons"), function(x) slot(x, "ID")))
  
  # Create dataframe with correct rownames
  (poly.df <- data.frame( ID=1:length(poly), row.names = poly_id))    
  
  # Try coersion again and check class
  poly <- SpatialPolygonsDataFrame(poly, poly.df)
}


##################################################################


# Función que crea un hueco en un polígono con la forma de otro
# Reference: http://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe


AddHoleToPolygon <-function(poly,hole){
  # invert the coordinates for Polygons to flag it as a hole
  coordsHole <-  hole@polygons[[1]]@Polygons[[1]]@coords
  newHole <- Polygon(coordsHole,hole=TRUE)
  
  # punch the hole in the main poly
  listPol <- poly@polygons[[1]]@Polygons
  listPol[[length(listPol)+1]] <- newHole
  punch <- Polygons(listPol,poly@polygons[[1]]@ID)
  
  # make the polygon a SpatialPolygonsDataFrame as the entry
  new <- SpatialPolygons(list(punch),proj4string=poly@proj4string)
  new <- SpatialPolygonsDataFrame(new,data=as(poly,"data.frame"))
  
  return(new)
}


###################################################################################

# HOLE FREE FUNCTION

hole_free <- function(x){
  BCp <- slot(x, "polygons")
  holes <- lapply(BCp, function(y) sapply(slot(y, "Polygons"), slot, "hole")) 
  res <- lapply(1:length(BCp), function(i) slot(BCp[[i]], "Polygons")[!holes[[i]]]) 
  IDs <- row.names(x)
  SpatialPolygons(lapply(1:length(res), function(i) Polygons(res[[i]], ID=IDs[i])), proj4string = CRS(proj4string(x))) 
}

###################################################################################

# SHAPE TO RASTER FUNCTION

shp2raster <- function(shp, mask.raster, label, value, transform = FALSE, proj.from = NA,
                       proj.to = NA, map = TRUE) {
  require(raster, rgdal)
  
  # use transform==TRUE if the polygon is not in the same coordinate system as
  # the output raster, setting proj.from & proj.to to the appropriate
  # projections
  if (transform == TRUE) {
    proj4string(shp) <- proj.from
    shp <- spTransform(shp, proj.to)
  }
  
  # convert the shapefile to a raster based on a standardised background
  # raster
  r <- rasterize(shp, mask.raster)
  # set the cells associated with the shapfile to the specified value
  r[!is.na(r)] <- value
  # merge the new raster with the mask raster and export to the working
  # directory as a tif file
  r <- mask(merge(r, mask.raster), mask.raster, filename = label, format = "GTiff",
            overwrite = T)
  
  # plot map of new raster
  if (map == TRUE) {
    plot(r, main = label, axes = F, box = F)
  }
  
  names(r) <- label
  return(r)
}


################################################################################
#                     FUNCTIONS FOR CLEANING GEOMETRIES                        #
###############################################################################

#Clean SpatialPoints (from polygons of Natural parks) -remove other treatments and get effective boundaries-
#Clean SpatialPoints (from polygons of Natural parks) -remove other treatments and get effective boundaries-
clean_treatments <- function(x, polygon, points_sp, points_border, shape){
  print(x$ID)
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
      .[!.@data$nn.dist < 1000, ]
    # dif <- gDifference(shape, x) %>% as("SpatialLines") %>% as("SpatialPoints")
    # knn_final <- get.knnx(coordinates(dif), coordinates(sp_border), k = 1, algorithm = "kd_tree") %>%
    #   data.frame()
    # sp_final <- SpatialPointsDataFrame(sp_border, knn_final, proj4string = CRS("+init=epsg:3857")) %>%
    #   .[!.@data$nn.dist < 500, ]
    
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
      .[!.@data$nn.dist < 1000, ]
    # dif <- gDifference(shape, x) %>% as("SpatialLines") %>% as("SpatialPoints")
    # knn_final <- get.knnx(coordinates(dif), coordinates(sp_border), k = 1, algorithm = "kd_tree") %>%
    #   data.frame()
    # sp_final <- SpatialPointsDataFrame(sp_border, knn_final, proj4string = CRS("+init=epsg:3857")) %>%
    #   .[!.@data$nn.dist < 500, ]
  }
  
}
#Clean SpatialPoints (from polygons of Natural parks) -remove points ONLY near to national frontiers and border points-
clean_treatments_border <- function(x, points_border){
  sp <- x %>% as("SpatialLines") %>% as("SpatialPoints")
  knn_border <- get.knnx(coordinates(points_border), coordinates(sp), k = 1, algorithm = "kd_tree") %>%
    data.frame(.)
  sp_final <- SpatialPointsDataFrame(sp, knn_border, proj4string = CRS("+init=epsg:3857")) %>%
    .[!.@data$nn.dist < 5000, ] 
  
}


################################################################################
#                           FUNCTIONS FOR RASTERS                              #
###############################################################################


#Function to crop/mask raster and calculate distance using all cores from CPU.
calculate_distances_parallel <- function(points){
  if(length(points)  > 2 & typeof(points) == "S4"){
    crop(res[[1]], buffer_areas) %>%
      mask(buffer_areas) %>%
      clusterR(.,distanceFromPoints, args = list(xy = points)) %>%
      mask(buffer_areas) %>%
      resample(res[[1]])
  } else
    return(0)
}


###################################################################################


# PROCESING RASTERS FUNCTION: crop, mask, make calculations and extract raster data

processing_rasters <- function(layer.list, ext, shape){
  layer.list %>%
    lapply(setExtent, ext) %>%
    lapply(crop, shape) %>%
    stack() %>% 
    mask(shape)
}



################################################################################
#                             TABLES MANAGMENT                                #
###############################################################################

rd_to_df <- function(list, dataframe){
  rd <- lapply(list, "[", "tabl3.str") %>%
    lapply(as.data.frame) %>%
    lapply( "[", 3 , ) %>%
    ldply() %>% mutate(N_l = unlist(lapply(list, "[", "N_l"))) %>%
    mutate(N_r = unlist(lapply(list, "[", "N_r"))) %>%
    mutate(N = N_l + N_r) %>%
    mutate(bws = unlist(lapply(list, function(x) x$bws[1, 1])))
  
  defo_mean <- mapply(function(x, y){
    y %>%
      filter(abs(dist_disc) <= x$bws[1, 1] & treatment == 0) %>% 
      dplyr::summarize(mean = mean(loss_sum))
  }, x = list , y = dataframe, SIMPLIFY = F) %>% unlist()
  
  df <- rd %>% cbind(., defo_mean) %>% t() %>% 
   as.data.frame() 
  names(df) <- NULL
  names(df) <- c("LineaNegra", "Doble_LN_PNN", "Doble_LN_Resg", "Triple")
  
  # %>% dplyr::rename(LineaNegra = Linea_Negra.mean,
  #                                     Doble_LN_PNN = Proteccion_Doble_LN_&_PNN.mean, Doble_LN_Resg = Proteccion_Doble_LN_&_Resg.mean,
  #                                     Triple = Triple_Proteccion_LN.mean)
  row.names(df) <- c("Territorio","Tratamiento", "StdErr", "Z", "p", "CI_l", "CI_u", "N_left","N_right", "N", "bws", "Media control")
  return(df)
}

rd_to_df_2 <- function(list, dataframe){
  rd <- lapply(list, "[", "tabl3.str") %>%
    lapply(as.data.frame) %>%
    lapply( "[", 3 , ) %>%
    ldply() %>% mutate(N_l = unlist(lapply(list, "[", "N_l"))) %>%
    mutate(N_r = unlist(lapply(list, "[", "N_r"))) %>%
    mutate(N = N_l + N_r) %>%
    mutate(bws = unlist(lapply(list, function(x) x$bws[1, 1])))
  
  lights_mean <- mapply(function(x, y){
    y %>%
      filter(abs(dist_disc) <= x$bws[1, 1] & treatment == 0) %>% 
      dplyr::summarize(mean = mean(F182013.v4c_web.stable_lights.avg_vis))
  }, x = list , y = dataframe, SIMPLIFY = F) %>% unlist()
  
  df <- rd %>% cbind(., lights_mean) %>% t() %>% 
    as.data.frame() 
  names(df) <- NULL
  names(df) <- c("LineaNegra", "Doble_LN_PNN", "Doble_LN_Resg", "Triple")
  
  # %>% dplyr::rename(LineaNegra = Linea_Negra.mean,
  #                                     Doble_LN_PNN = Proteccion_Doble_LN_&_PNN.mean, Doble_LN_Resg = Proteccion_Doble_LN_&_Resg.mean,
  #                                     Triple = Triple_Proteccion_LN.mean)
  row.names(df) <- c("Territorio","Tratamiento", "StdErr", "Z", "p", "CI_l", "CI_u", "N_left","N_right", "N", "bws", "Media control")
  return(df)
}


