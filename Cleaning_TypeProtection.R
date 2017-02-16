#################################################
#                                               #
#             CLEANING POLYGONS                 #
#           (TYPES OF PROTECTION)               #
#                                               #
#################################################  


# Open protection levels
# ----------------------
setwd("~/Dropbox/Linea_Negra_R/Data/")
protect_areas <- list.files("Zonas_Protegidas_LN") %>%
  .[str_detect(., ".shp")] %>% strsplit(., ".", fixed = T) %>%
  unlist() %>% .[!str_detect(., "shp")] %>%
  lapply(., function(x){
    readOGR(dsn = "Zonas_Protegidas_LN" , layer = x) %>%
      spTransform(CRS=CRS("+init=epsg:3857"))
  })

names <- list.files("Zonas_Protegidas_LN") %>%
  .[str_detect(., ".shp")] %>% strsplit(., ".", fixed = T) %>%
  unlist() %>% .[!str_detect(., "shp")] 

#Graphs
plot(protect_areas[[4]])
plot(protect_areas[[3]], add = T, border = "red")
plot(protect_areas[[1]], add = T, border = "blue")
plot(protect_areas[[2]], add = T, border = "orange")
plot(protect_areas[[5]], add = T, border = "yellow")


#Prepare data
protect_areas_p <- lapply(protect_areas, function(x){
  as(x, "SpatialLines") %>% as("SpatialPoints")
})

protect_areas_buffer <- protect_areas %>%
  lapply(., gBuffer, width = 50000, byid = T) %>%
  lapply(spTransform, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(protect_areas_buffer[[4]])
plot(protect_areas_buffer[[3]], add = T, border = "red")
plot(protect_areas_buffer[[1]], add = T, border = "blue")
plot(protect_areas_buffer[[2]], add = T, border = "orange")
plot(protect_areas_buffer[[5]], add = T, border = "yellow")

#Clean geometries (remove shorelines)
list_polygons_clean_border <- lapply(protect_areas, clean_treatments_border, points_border = colombia_p) %>%
  lapply(spTransform, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


#Does it work?
plot(protect_areas_p[[4]], cex = 0.1, pch = 19)
plot(list_polygons_clean_border[[4]], add = T, cex = 0.1, pch = 19, col = "red")

plot(protect_areas_p[[1]], cex = 0.1, pch = 19)
plot(list_polygons_clean_border[[1]], add = T, cex = 0.1, pch = 19, col = "red")

#Calculate distances

beginCluster()
system.time(mask_distance <- mapply(calculate_distances_parallel,
                                    buffer = protect_areas_buffer, 
                                    points = list_polygons_clean_border))
endCluster()

#Identify treatment/control cells
cells_territories <- protect_areas %>%
  lapply(spTransform, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  lapply(., function(x){
    cellFromPolygon(res_mask_ln, x)
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
list_dataframes[[5]]$treatment <- ifelse(list_dataframes[[5]]$ID %in% unlist(cells_territories[[5]]), 1, 0)

#3. Append all elements of the list (1: LN, 2: Resg, 3: PNN) 
distance_dataframe <- do.call(rbind, list_dataframes)
distance_dataframe$buffer_id <- rep(names(list_dataframes), sapply(list_dataframes, nrow)) #identify cells from buffers

#4. Export
setwd("~/Dropbox/Linea_Negra_R/Data/Dataframes/")
write.csv(distance_dataframe, "distance_dataframe_protected_areas.csv", row.names = FALSE)


