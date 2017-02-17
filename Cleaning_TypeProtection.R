library(raster)
library(rgdal)
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

setwd("~/GitHub/linea_negra/")
source("Functions.R")

#Directories
#data <-  setwd("~/Dropbox/Linea_Negra_R/Data/")
data <- setwd("~/Dropbox/BANREP/Linea_Negra_R/Data/")


#################################################
#                                               #
#             CLEANING POLYGONS                 #
#           (TYPES OF PROTECTION)               #
#                                               #
#################################################  

# Open basic geometries
# ----------------------

# Create a border colombia geometry (to clean shoreline)
# Division by departaments
colombia <- readOGR(dsn = "Colombia_Deptos", layer = "Colombia_Deptos") %>%
  spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  spTransform(CRS=CRS("+init=epsg:3857")) %>%
  .[!(.@data$DPTO_CCDGO == 88) , ] %>%
  gUnaryUnion(.) %>%
  hole_free()
# plot(colombia)

# Resguardos de Magdalena, La Guajira y Cesar
# ------------------------------------------------
resguardos <- readOGR(dsn = "IGAC", layer = "Resguardos_Selected_LNegra") %>%
  spTransform(CRS=CRS("+init=epsg:3857"))

resguardos_ln <- raster::intersect(resguardos, linea_negra_proj) %>%
  spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
resguardos_ln_proj <- resguardos_ln %>% spTransform(CRS=CRS("+init=epsg:3857"))

# PNN dentro de la Línea Negra
# ---------------------------------
pnn <- readOGR(dsn = "WDPA/Protected_Areas/", layer = "WDPA_Jan2017_COL-shapefile-polygons", stringsAsFactors = F) %>% 
  spTransform(CRS=CRS("+init=epsg:3857")) %>%
  .[!.@data$STATUS_YR > 2012 & !.@data$GIS_AREA < 1 & .@data$MANG_AUTH == "Parques Nacionales Naturales de Colombia" , ]

pnn_ln <- raster::intersect(pnn, linea_negra_proj)

# Polígono Línea Negra
# ---------------------------------
linea_negra <- readOGR(dsn = "Linea_Negra", layer = "Linea_Negra_Polygon") %>%
  spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

linea_negra_proj <- linea_negra %>% spTransform(CRS=CRS("+init=epsg:3857"))

# Densify lines using geosphere package
lineanegra_p <- linea_negra %>% as("SpatialLines") %>% as("SpatialPoints")  # _p for points
lineanegra_dens <- makePoly(linea_negra, interval=1000, sp = T)  %>% # The higuer the interval the lower the number of new points
  spTransform(CRS=CRS("+init=epsg:3857")) 


# Open protection levels
# ----------------------

#Set directories
# setwd("~/Dropbox/Linea_Negra_R/Data/")
setwd("~/Dropbox/BANREP/Linea_Negra_R/Data/")

protect_areas <- list.files("Zonas_Protegidas_LN") %>%
  .[str_detect(., ".shp")] %>% strsplit(., ".", fixed = T) %>%
  unlist() %>% .[!str_detect(., "shp")] %>% .[c(1:2, 5)] %>%
  lapply(., function(x){
    readOGR(dsn = "Zonas_Protegidas_LN" , layer = x) %>%
      spTransform(CRS=CRS("+init=epsg:3857"))
  }) %>%
  c(., lineanegra_dens)

names <- list.files("Zonas_Protegidas_LN") %>%
  .[str_detect(., ".shp")] %>% strsplit(., ".", fixed = T) %>%
  unlist() %>% .[!str_detect(., "shp")]  %>% .[c(1:2, 5)] %>%
  c(., "Linea_Negra")

#Graphs
plot(protect_areas[[3]], add = T, border = "red")
plot(protect_areas[[1]], add = T, border = "blue")
plot(protect_areas[[2]], add = T, border = "orange")
plot(protect_areas[[4]], add = T, border = "yellow")

#Prepare data
protect_areas_p <- lapply(protect_areas, function(x){
  as(x, "SpatialLines") %>% as("SpatialPoints")
})

buffer_areas <- gBuffer(protect_areas[[4]], width = 50000) %>%
  spTransform(.,CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


# protect_areas_buffer <- protect_areas %>%
#   lapply(., gBuffer, width = 50000, byid = T) %>%
#   lapply(spTransform, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Visualize types of protection
plot(protect_areas_buffer[[3]], add = T, border = "red")
plot(protect_areas_buffer[[1]], add = T, border = "blue")
plot(protect_areas_buffer[[2]], add = T, border = "orange")
plot(protect_areas_buffer[[5]], add = T, border = "yellow")

#Clean geometries
#Load functions
setwd("~/GitHub/linea_negra/")
source("Functions.R")

#Prepare data for cleaning (intersections and unions)

#National parks and indigenous to points
resguardos_p <- resguardos %>% as("SpatialLines") %>% as("SpatialPoints")
pnn_ln_p <- pnn_ln %>% as("SpatialLines") %>% as("SpatialPoints")
resg_pnn_p <- rbind(resguardos_p, pnn_ln_p)

#Colombia points
colombia_p <- colombia %>% as("SpatialLines") %>% as("SpatialPoints")

#National parks and indigenous (intersection)
resguardos_hole_free <- hole_free(resguardos) %>% gBuffer(0.0001, byid = F)
pnn_ln_hole_free <- hole_free(pnn_ln) %>% gBuffer(0.0001, byid = F)
resg_pnn <- raster::union(resguardos_hole_free, pnn_ln_hole_free) # raster is the package and union is the name of the specific function within the package. :: helps to access the exact function from that specific package


############################### WARNING ######################################
# CLEANING ALL CASES LISTED IN PROTECT_AREAS IN BATCH WILL LEAD TO ERRORS   #
# SINCE ALL GEOMETRIES HAVE COMMON BOUNDARIES. THE BATCH CLEANING GENERATES #
# UNREPRESENTATIVE GEOMETRIES WITHOUT BOUNDARIES                            #
#############################################################################

#Clean by case
polygons_contiguous <- list(resguardos_ln_proj, pnn_ln, protect_areas[[1]], resg_pnn)
points_contiguous <- list(resguardos_p, pnn_ln_p, protect_areas_p[[1]], resg_pnn_p)

list_polygons_clean_border <- mapply(function(x, y, z){
  clean_treatments(
    x = x,
    polygon = y,
    points_sp = z,
    points_border = colombia_p)
  }, x = protect_areas, y = polygons_contiguous, z = points_contiguous)

#Does it work?
plot(protect_areas[[1]])
plot(list_polygons_clean_border[[1]], add = T, cex = 0.1, pch = 19, col = "red")

plot(protect_areas[[4]])
plot(list_polygons_clean_border[[4]], add = T, cex = 1, pch = 19, col = "red")


######################## DISTANCE CALCULATIONS ###########################

################################ WARNING ##################################
# We need a raster as a template. We will use deforestation raster masked #
# to the buffer to mantain the size and ID's of rasters. For that we need #
# to have all the geometries in WGS84 (re-project is hard!)               #
###########################################################################


# Loss-Year Brick Linea Negra 1km
res <- brick("loss_year_brick_1km.tif")
#plot(res)

#Project all geometries to WGS84 (easier than to reproject raster)
clean_frontiers <- list_polygons_clean_border %>%
  lapply(., function(x){
    if(typeof(x) == "S4" & length(x) > 0){
      sp <- spTransform(x, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      return(sp)
    } else { 
      return(0)
    }
  })


#I only use the buffer of LN since all the geometries are inside. 
# protect_areas_buffer <- protect_areas %>%
#   lapply(., function(x){
#     if(typeof(x) == "S4" & length(x) > 0){
#       sp <- gBuffer(x, width = 50000) %>%
#       spTransform(x, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#       return(sp)
#     } else { 
#       return(0)
#     }
#   })



beginCluster()
system.time(mask_distance <- mapply(calculate_distances_parallel,
                                    points = clean_frontiers))
endCluster()

#Identify treatment/control cells
cells_territories <- protect_areas %>%
  lapply(spTransform, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  lapply(., function(x){
    cellFromPolygon(res[[1]], x)
  })

lapply(cells_territories, function(x){
  unlist(x) %>% length()
})


#################################################
#                                               #
#         EXTRACT DISTANCE DATA FROM            #
#                       RASTERS                 #
#                                               #
#################################################

#1. Extract distance as data frame per buffer (list element)
list_dataframes <- pblapply(mask_distance, as.data.frame, xy = T)
names(list_dataframes) <- names

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
list_dataframes[[4]]$treatment <- ifelse(list_dataframes[[4]]$ID %in% unlist(cells_territories[[4]]), 1, 0)

#3. Append all elements of the list (1: LN, 2: Resg, 3: PNN) 
distance_dataframe <- do.call(rbind, list_dataframes)
distance_dataframe$buffer_id <- rep(names(list_dataframes), sapply(list_dataframes, nrow)) #identify cells from buffers
distance_dataframe <- distance_dataframe[complete.cases(distance_dataframe[, c("layer")]), ]

#4. Export
setwd(str_c(data, "/", "Dataframes"))
write.csv(distance_dataframe, "distance_dataframe_protected_areas.csv", row.names = FALSE)


