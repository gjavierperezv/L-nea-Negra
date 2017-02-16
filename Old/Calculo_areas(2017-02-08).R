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
library(maptools)

# #################################################################
# FUNCIÓN PARA CREAR EL DATA SLOT
# -------------------------------
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
# Refference: http://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe
# ----------------------------------------------------------------
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
# PROCESING RASTERS FUNCTION: crop, mask, make calculations and extract raster data
# ---------------------------------------------------------------------------------
processing_rasters <- function(layer.list, ext, shape){
  layer.list %>%
    lapply(setExtent, ext) %>%
    lapply(crop, shape) %>%
    stack() %>% 
    mask(shape)
}
###################################################################################

setwd("~/Dropbox/Linea_Negra_R/Data/")

###############################
# 1. POLÍGONO DE LA LÍNEA NEGRA
###############################
  
  # 1.1. Línea Negra
  # ----------------
  linea_negra <- readOGR(dsn = "Linea_Negra", layer = "Linea_Negra_Polygon") %>%
    spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  linea_negra_proj <- linea_negra %>% spTransform(CRS=CRS("+init=epsg:3857"))
  #plot(linea_negra_proj)
  linea_negra_proj@data$AREA_KM2_ln <- gArea(linea_negra_proj, byid = T)/1e6
  # Otra alternativa:
  # rgeos::gArea(linea_negra_proj)/1e6
  sum(linea_negra_proj@data$AREA_KM2_ln)
  plot(linea_negra_proj)

  # 1.2. Buffer Línea Negra (50 km2)
  # --------------------------------
  linea_negra_buffer <- gBuffer(linea_negra_proj, width = 50000, byid = T)
  dev.off() # Clear all plots
  plot(linea_negra_buffer)
  plot(linea_negra_proj, add = T)

#########################################
# 2. POLÍGONO DE LOS RESGUARDOS INDÍGENAS
#########################################
  
  # 2.1. Resguardos de Magdalena, La Guajira y Cesar
  # ------------------------------------------------
  resguardos <- readOGR(dsn = "IGAC", layer = "Resguardos_Selected_LNegra") 
  resguardos <-  spTransform(resguardos, CRS=CRS("+init=epsg:3857"))
  resguardos@data$AREA_KM2_resg <- gArea(resguardos, byid = T)/1e6
  dev.off() # Clear all plots
  plot(linea_negra_proj)
  plot(resguardos, add = T)
  
  # 2.2. Resguardos de la Línea Negra
  # ---------------------------------
  resguardos_ln <- raster::intersect(resguardos, linea_negra_proj) %>%
    spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  resguardos_ln_proj <- resguardos_ln %>% spTransform(CRS=CRS("+init=epsg:3857"))
  
  resguardos_ln_proj@data$AREA_KM2_resln <- gArea(resguardos_ln_proj, byid = T)/1e6
  plot(linea_negra_proj)
  plot(resguardos_ln_proj, add = T, border = "red")
  sum(resguardos_ln_proj@data$AREA_KM2_resln)

#######################################################
# 3. POLÍGONO DE LOS PARQUES NACIONALES NATURALES (PNN)
#######################################################
  
  # 3.1. PNN de toda Colombia
  # -------------------------
  pnn <- readOGR(dsn = "WDPA/Protected_Areas/", layer = "WDPA_Jan2017_COL-shapefile-polygons") %>% 
    spTransform(CRS=CRS("+init=epsg:3857")) %>%
    .[!.@data$STATUS_YR > 2012 & !.@data$GIS_AREA < 1 & .@data$MANG_AUTH == "Parques Nacionales Naturales de Colombia" , ]
  #plot(pnn)
  
  # 3.2. PNN dentro de la Línea Negra
  # ---------------------------------
  pnn_ln <- raster::intersect(pnn, linea_negra_proj)
  dev.off() # Clear all plots
  pnn_ln@data$AREA_KM2_pnn <- gArea(pnn_ln, byid = T)/1e6
  plot(linea_negra_proj)
  plot(pnn_ln, add = T, border = "red")
  sum(pnn_ln@data$AREA_KM2_pnn)

########################
# 4. TIPOS DE PROTECCIÓN
########################
  
  # 4.1. SIN PROTECCIÓN ((Fuera de la Línea Negra (dentro de un Buffer de 50 km2))
  # ------------------------------------------------------------------------------
    colombia <- readOGR(dsn = "Colombia_Deptos", layer = "Colombia_Deptos")
    colombia <- spTransform(colombia, CRS=CRS("+init=epsg:3857"))
    plot(colombia)
    slotNames(colombia)
    
    # Disolving
      croquis_colombia <- gUnaryUnion(colombia, id = colombia@data$DEPTOCCDGO)
      plot(croquis_colombia) 
      plot(linea_negra_buffer, add = T)
    
    # Intersect
       buffer_continental <- raster::intersect(croquis_colombia, linea_negra_buffer)
       plot(croquis_colombia)
       plot(buffer_continental, add = T, border = "red")
       plot(buffer_continental)
    
     # Creando el polígono 
       buffer_continental@data$AREA_KM2_buffercontin <- gArea(buffer_continental, byid = T)/1e6
       sum(buffer_continental@data$AREA_KM2_buffercontin)
      
     # Calculando el buffer continantal de la Linea Negra
      buffer_continental <- spTransform(buffer_continental, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
     
     # Calculando el polígono entre la Línea Negra y el buffer continental
       sin_proteccion <- AddHoleToPolygon(buffer_continental[1,], linea_negra[1,]) 
       plot(sin_proteccion, col = "blue", border = "red")
    
    # Calculando el área del polígono sin protección
       sin_proteccion <- DataSlot(sin_proteccion) %>%
         spTransform(sin_proteccion, CRS=CRS("+init=epsg:3857"))
       sin_proteccion@data$AREA_KM2_sinprotecc <- gArea(sin_proteccion, byid = T)/1e6
       sum(sin_proteccion@data$AREA_KM2_sinproteccion)
       
    # Exportando el polígono en formato shape
       writeOGR(sin_proteccion, dsn = 'Zonas_Protegidas_LN', layer = 'Sin_Proteccion_LN', driver = "ESRI Shapefile", overwrite_layer = T)
  
  # 4.2. PROTECCIÓN SENCILLA (LN - PNN - RESGUARDOS)
  # -----------------------------------------------
     
      # Se crea la Unión entre PNN y Resguardos
      resg_union_pnn <- gUnion(resguardos_ln_proj, pnn_ln, byid = T)
      plot(resg_union_pnn, col = "blue", border = "red")
      
      # Creando la zona de protección individual
      proteccion_sencilla <- gDifference(linea_negra_proj, resg_union_pnn, byid = F, drop_lower_td = T)
      plot(proteccion_sencilla, col = "blue", border = "red")
      
      # Creando el Data Slot para luego poder covertirlo a SpatialPOlygonDataFrame y así poder exportarlo a shape
      proteccion_sencilla <- DataSlot(proteccion_sencilla)
      proteccion_sencilla@data$AREA_KM2_protecsencilla <- gArea(proteccion_sencilla, byid = T)/1e6
      sum(proteccion_sencilla@data$AREA_KM2_protecsencilla)
      
      # Exportando el polígono en formato shape
      #proteccion_sencilla <- SpatialPolygonsDataFrame(proteccion_sencilla)
      writeOGR(proteccion_sencilla, dsn = 'Zonas_Protegidas_LN', layer = 'Proteccion_Sencilla_LN', driver = "ESRI Shapefile", overwrite_layer = T)

  # 4.3. DOBLE PROTECCIÓN (LINEA NEGRA+RESGUARDOS = RESGUARDOS - PNN)
  # -----------------------------------------------------------------
  # Opción (polígono correcto pero sin el slot data)
  doble_resg_menos_pnn <- gDifference(resguardos_ln_proj, pnn_ln)
  plot(linea_negra_proj)
  plot(doble_resg_menos_pnn, add = T,  col = "blue", border = "red")
  
    # Creando el Data Slot
      doble_resg_menos_pnn <- DataSlot(doble_resg_menos_pnn)
      slotNames(doble_resg_menos_pnn)
      class(doble_resg_menos_pnn)
      
    # Calculando el área
      doble_resg_menos_pnn@data$AREA_KM2_resmenospnn <- gArea(doble_resg_menos_pnn, byid = T)/1e6
      sum(doble_resg_menos_pnn@data$AREA_KM2_resmenospnn)
    
    # Graficando  
      plot(linea_negra_proj)
      plot(resguardos_ln_proj, add = T, border = "green")
      plot(resg_menos_pnn, add = T, border = "blue")
      plot(pnn_ln, add = T, border = "red")
      
  # 4.4. DOBLE PROTECCIÓN (LINEA NEGRA+PNN = PNN - RESGUARDOS)
  # -----------------------------------------------------------------
    doble_pnn_menos_resg <- gDifference(pnn_ln, resguardos_ln_proj, byid = F, drop_lower_td = T)
    plot(linea_negra_proj)
    plot(resguardos_ln_proj, add = T, border = "green")
    plot(pnn_ln, add = T, border = "red")
    plot(doble_pnn_menos_resg, add = T, col = "blue", border = "red")
    
     # Creando el Data Slot
      doble_pnn_menos_resg <- DataSlot(doble_pnn_menos_resg)
      
    # Calculando las áreas
    doble_pnn_menos_resg@data$AREA_KM2_pnnmenosresg <- gArea(doble_pnn_menos_resg, byid = T)/1e6
    sum(doble_pnn_menos_resg@data$AREA_KM2_pnnmenosresg)

  # 4.5. TRIPLE PROTECCIÓN (INTERSECCIÓN ENTRE PNN Y RESGUARDOS)
  # ------------------------------------------------------------
  triple_resg_in_pnn <- raster::intersect(resguardos_ln_proj, pnn_ln)
  triple_resg_in_pnn@data$AREA_KM2_resinpnn <- gArea(triple_resg_in_pnn, byid = T)/1e6
  triple_resg_in_pnn <- gUnaryUnion(triple_resg_in_pnn)
  dev.off() # Clear all plots
  plot(linea_negra_proj)
  #plot(resguardos_ln_proj, add = T)
  #plot(pnn_ln, add = T)
  plot(triple_resg_in_pnn, add = T, col = "blue", border = "red")
  sum(triple_resg_in_pnn@data$AREA_KM2_resinpnn)

###################################################################  
# 5. CÁLCULO E LAS ÁREAS DE AFECTACIÓN POR LOS DIFERENTES FENÓMENOS
###################################################################
  
    # Gráficos de las diferentes áreas de protección
    # ##############################################
      # 1. Sin protección
      plot(sin_proteccion, col = "blue", border = "red")  
      plot(linea_negra_proj, add = T)
      # 2. Protección sencilla
      plot(linea_negra_proj)
      plot(proteccion_sencilla, add = T, col = "blue", border = "red")
      # 3. Protección doble (LINEA NEGRA + RESGUARDOS (Resguardos - PNN))
      plot(linea_negra_proj)
      plot(doble_resg_menos_pnn, add = T, col = "blue", border = "red")
      # 4. Protección doble (LINEA NEGRA + PNN (PNN - Resguardos))
      plot(linea_negra_proj)
      plot(doble_pnn_menos_resg, add = T, col = "blue", border = "red")
      #5. Triple protección (LINEA NEGRA + RESGUARDOS + PNN)
      plot(linea_negra_proj)
      plot(triple_resg_in_pnn, add = T, col = "blue", border = "red")
    
  # ------------------
  # 5.1. DEFORESTACIÓN
  # ------------------
    # Leyendo la información de deforestación (Raster 1 km2)
    res <- brick("~/Dropbox/Linea_Negra_R/Data/WDPA/loss_year_brick_1km.tif/")
    plot(res)
    
    # Calculando el área de deforestación acumulada entre 2001 y 2012
      defo_agregada <- stackApply(res[[c(2:13)]], 1, fun = sum,  filename = "loss_year_brick_1km.tif",
                                format = "GTiff",
                                options = "INTERLEAVE=BAND", 
                                progress = "text", overwrite = T)
      
    # Creando una lista y el cáclculo  
    
      
#################################################################################################      
     
          
      # Calculando el área deforestada
      list <- c(sin_proteccion, proteccion_sencilla, doble_resg_menos_pnn, doble_pnn_menos_resg, triple_resg_in_pnn)
      defo <- lapply(list, function(x){
        sp <- spTransform(x, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        raster::extract(res,
                        sp,
                        fun = sum, na.rm = T, df = T) %>%
          mutate(area_pixel = length(unlist(cellFromPolygon(res[[1]], sp))))
      }) %>% ldply() 
      
      defo_tot <- defo %>%
        mutate(., loss_sum = rowSums(defo[, c(4:length(names(defo)) - 1 )])) %>%
        dplyr::select(., c(ID, loss_sum, area_pixel)) %>% mutate(.,loss_sum = (loss_sum * 100) / 14) %>%
        mutate(defo_rate = loss_sum/area_pixel)
      
      xtable(defo_tot)
  # ----------------------------------
  # 5.2. POBLAMIENTO (LUCES NOCTURNAS)
  # ----------------------------------
  
  
  # ---------
  # 5.3. VÍAS
  # ---------

  
  
  # -------------------------------
  # 5.4. CONFLICTO DE USO DEL SUELO
  # -------------------------------
 
  
  
   
  # ------------------------
  # 5.5. AFECTACIÓN RELATIVA
  # ------------------------
  
  
  
  dev.off() # Clear all plots
  savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")
    
    
    




