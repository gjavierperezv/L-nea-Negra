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
library(spatstat)
library(xtable)

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
# HOLE FREE FUNCTION
# ------------------
hole_free <- function(x){
  BCp <- slot(x, "polygons")
  holes <- lapply(BCp, function(y) sapply(slot(y, "Polygons"), slot, "hole")) 
  res <- lapply(1:length(BCp), function(i) slot(BCp[[i]], "Polygons")[!holes[[i]]]) 
  IDs <- row.names(x)
  SpatialPolygons(lapply(1:length(res), function(i) Polygons(res[[i]], ID=IDs[i])), proj4string = CRS(proj4string(x))) 
}

###################################################################################
# SHAPE TO RASTER FUNCTION
# ------------------------
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
      
     # Calculando el buffer continental de la Linea Negra
       buffer_continental <- spTransform(buffer_continental, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
     
     # Calculando el polígono entre la Línea Negra y el buffer continental
       sin_proteccion <- AddHoleToPolygon(buffer_continental[1,], linea_negra[1,]) 
       plot(sin_proteccion, col = "blue", border = "red")
    
       # Quitándole la línea negra
       sin_proteccion <- gDifference(buffer_continental, linea_negra, byid = F, drop_lower_td = T) %>%
         spTransform(., CRS=CRS("+init=epsg:3857")) %>% # Proyectando la capa (metros)
         gDifference(., pnn, byid = F, drop_lower_td = T)  # Quitándole los Parques Nacionales 
         # Leyendo la información de los Resguardos Indígenas
         resguardos_colombia_proj <- readOGR(dsn = "IGAC", layer = "Resguardos_Indigenas") %>% 
            spTransform(., CRS=CRS("+init=epsg:3857"))
            
       sin_proteccion <- gDifference(sin_proteccion, resguardos_colombia_proj, byid = F, drop_lower_td = T)  # Quitándole los Resguardos Indígenas  
        # Mapa del polígono
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
    
    # Exportando el shape file
    writeOGR(doble_resg_menos_pnn, dsn = 'Zonas_Protegidas_LN', layer = 'Proteccion_Doble_LN_&_Resg', driver = "ESRI Shapefile", overwrite_layer = T)

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

    # Exportando el shape file
    writeOGR(doble_pnn_menos_resg, dsn = 'Zonas_Protegidas_LN', layer = 'Proteccion_Doble_LN_&_PNN', driver = "ESRI Shapefile", overwrite_layer = T)
    
  # 4.5. TRIPLE PROTECCIÓN (INTERSECCIÓN ENTRE PNN Y RESGUARDOS)
  # ------------------------------------------------------------
  triple_resg_in_pnn <- raster::intersect(resguardos_ln_proj, pnn_ln)
  triple_resg_in_pnn <- gUnaryUnion(triple_resg_in_pnn)
  triple_resg_in_pnn <- DataSlot(triple_resg_in_pnn)
  triple_resg_in_pnn@data$AREA_KM2_resinpnn <- gArea(triple_resg_in_pnn, byid = T)/1e6
  dev.off() # Clear all plots
  plot(linea_negra_proj)
  #plot(resguardos_ln_proj, add = T)
  #plot(pnn_ln, add = T)
  plot(triple_resg_in_pnn, add = T, col = "blue", border = "red")
  sum(triple_resg_in_pnn@data$AREA_KM2_resinpnn)

  # Exportando el polígono en formato shape
  writeOGR(triple_resg_in_pnn, dsn = 'Zonas_Protegidas_LN', layer = 'Triple_Proteccion_LN', driver = "ESRI Shapefile", overwrite_layer = T)
  
#############################################################################################  
# 5. VARIABLES DE RESULTADO: CÁLCULO DE LAS ÁREAS DE AFECTACIÓN POR LOS DIFERENTES FENÓMENOS
#############################################################################################
  
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
    res_df <- as.data.frame(res, xy = T)
    res_df$ID <- row.names(res_df)
    setwd("~/Dropbox/Linea_Negra_R/Data/Dataframes/")
    write.csv(res_df, "deforestation_dataframe.csv", row.names = FALSE)
    
    # Calculando el área de deforestación acumulada entre 2001 y 2012
      defo_agregada <- stackApply(res[[c(2:13)]], 1, fun = sum,  filename = "loss_year_brick_1km.tif",
                                format = "GTiff",
                                options = "INTERLEAVE=BAND", 
                                progress = "text", overwrite = T)
      
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
        dplyr::select(., c(ID, loss_sum, area_pixel)) %>% 
        mutate(.,loss_sum = (loss_sum * 100) / 14) %>%
        mutate(defo_rate = loss_sum/area_pixel)
      
      xtable(defo_tot)
  # ----------------------------------
  # 5.2. POBLAMIENTO (LUCES NOCTURNAS)
  # ----------------------------------
    # Leyendo la infoamción de Luces Nocturnas
      setwd("~/Dropbox/Linea_Negra_R/Data/NOAA/AVER_AND_CLOUDFREE/Web_Stable_Average_Visible/")
      list_raster <- list.files() %>%
        lapply(raster)
      rasters_extent <- extent(list_raster[[1]]) #We need to put all rasters into the same extent (all have the same resolution). [] extracts a list, [[]] extracts elements within the list
      
      setwd("~/Dropbox/Linea_Negra_R/Data/")
        colombia <- readOGR(dsn = "Colombia_Deptos", layer = "Colombia_Deptos") %>%
        spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
        #spTransform(CRS=CRS("+init=epsg:3857")) %>%
        .[!(.@data$DPTO_CCDGO == 88) , ] %>%
        gUnaryUnion(.) %>%
        hole_free() 
        
      rasters_lights <- processing_rasters(list_raster, rasters_extent, colombia) %>%
        resample(res[[1]]) #Set both resolutions equal 
      plot(rasters_lights)
      
      #Extract light raster data 
      rasters_lights_df <- as.data.frame(rasters_lights, xy = T)
      rasters_lights_df$ID <- row.names(rasters_lights_df)
      setwd("~/Dropbox/Linea_Negra_R/Data/Dataframes/")
      write.csv(rasters_lights_df ,"rasters_lights_dataframe.csv", row.names = FALSE)
      
  
    # Calculando el área con influencia de Luces Nocturnas  
      list <- c(sin_proteccion, proteccion_sencilla, doble_resg_menos_pnn, doble_pnn_menos_resg, triple_resg_in_pnn)
      lights <- lapply(list, function(x){
        sp <- spTransform(x, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        raster::extract(rasters_lights,
                        sp,
                        fun = sum, na.rm = T, df = T) %>%
          mutate(area_pixel = length(unlist(cellFromPolygon(rasters_lights[[1]], sp))))
      }) %>% ldply() 
      
      lights_tot <- lights %>%
        mutate(., lights_sum = rowSums(lights[, c(1:length(names(lights)) - 1 )])) %>%
        dplyr::select(., c(ID, lights_sum, area_pixel)) %>% 
        mutate(.,lights_sum = (lights_sum * 100) / 14) %>%
        mutate(lights_rate = lights_sum/area_pixel)
      
  # ---------
  # 5.3. VÍAS
  # ---------
      # Leyendo las vías del DANE-MGN a 2005 (todas las vías)
      vias_mgn <- readOGR(dsn = "DANE_MGN", layer = "VIAS_TODAS_RANGO_LN_2005") %>%
        spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      vias_mgn_proj <- spTransform(vias_mgn, CRS=CRS("+init=epsg:3857"))
      plot(vias_mgn)
      plot(linea_negra, add = T, border = "red")
  
      # Leyendo las vías del IGAC-SIGOT 
        # Vias Tipo 1
        vias_sigot_1 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_1") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_1)
        # Vias Tipo 2
        vias_sigot_2 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_2") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_2)
        # Vias Tipo 3
        vias_sigot_3 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_3") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_3)
        # Vias Tipo 4
        vias_sigot_4 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_4") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_4)
        # Vias Tipo 5
        vias_sigot_5 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_5") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_5)
        # Vias Tipo 6
        vias_sigot_6 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_6") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_6)
        # Vias Tipo 7
        vias_sigot_7 <- readOGR(dsn = "IGAC/Vias_SIGOT_PorTipo", layer = "Vias_SIGOT_Tipo_7") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(vias_sigot_7)
        
      # Haciendo Merge de las vías del DANE-MGN y del IGAC-SIGOT
        # vias_dane_igac <- gUnion(vias_mgn, vias_sigot_4, byid = T)  # Se bloquea R Studio

      # extract all polygons and convert to a raster where cells
      # associated with vias_mgn have a value of 1 and everything else has a value of 0
      # vias_mgn_raster <- shp2raster(shp = vias_mgn, mask.raster = rasters_lights[[1]], label = "Vias Rasterizadas", transform = FALSE, value = 1)
      
      # Rasterizando la capa de vías
        vias_mgn_raster <- rasterize(vias_mgn, rasters_lights[[1]], vias_mgn@data$DPTO_NAREA)
        writeRaster(vias_mgn_raster, filename = "vias_mgn_raster_1km.tif", format = "GTiff", progress = "text", overwrite = T)  
        
        plot(vias_mgn_raster)
        plot(linea_negra, add = T, border = "red")
      
      # Calculando los kiómetros de vías en cada una de las zonas de análisis
      # Referencia: http://gis.stackexchange.com/questions/138861/calculating-road-density-in-r-using-kernel-density
        
        # OPCIÓN 1
        # --------
        # Convert SpatialLines to psp object using maptools library
        #psp_vias_mgn_proj <- as.psp(vias_mgn_proj)
        #plot(psp_vias_mgn_proj)
        # Pixellate with resolution
        #px_vias_mgn_proj <- pixellate(psp_vias_mgn_proj)
        #plot(px_vias_mgn_proj)
        # This can be converted to raster as desired
        #rLength <- raster(px_vias_mgn_proj)
        # Values:
        round(as.matrix(rLength),3)
      
        plot(pixellate(psp_vias_mgn_proj, what = "l"))
        
    # Otras opciones:
    # Referencia: http://gis.stackexchange.com/questions/119993/convert-line-shapefile-to-raster-value-total-length-of-lines-within-cell
       
        # OPCIÓN 2
        #---------
        # Need to convert to a line segment pattern object with maptools
        roadsPSP <- as.psp(as(vias_mgn_proj, 'SpatialLines'))
        # Calculate lengths per cell
        roadLengthIM <- pixellate.psp(roadsPSP)
        # Convert pixel image to raster in km
        roadLength <- raster(roadLengthIM / 1000, crs=projection(vias_mgn_proj))
        writeRaster(roadLength, filename = "road_density_km_km2.tif", format = "GTiff", progress = "text", overwrite = T)
        round(as.matrix(roadLength),3)
        sum(lengths.psp(roadsPSP))
        sum(pixellate(roadsPSP))
        
        # Plot
        spplot(roadLength, scales = list(draw=TRUE), xlab="x", ylab="y", 
               col.regions=colorRampPalette(brewer.pal(9, "YlOrRd")), 
               sp.layout=list("sp.lines", vias_mgn), 
               par.settings=list(fontsize=list(text=15)))

  # -------------------------------
  # 5.4. CONFLICTO DE USO DEL SUELO
  # -------------------------------
    # Leyendo los polígonos de Conflicto de uso del suelo
        setwd("~/Dropbox/Linea_Negra_R/Data/")
        
        # Conflictos por subutilización
        # -----------------------------
          conflicto_sub <- readOGR(dsn = "IGAC", layer = "Conflictos_Subutlizacion") %>%
            spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          plot(conflicto_sub)
          plot(linea_negra, add = T, border = "red")
    
          # Rasterizando la capa de conflicto por Subutilización
          conflicto_sub_raster <- rasterize(conflicto_sub, rasters_lights[[1]], conflicto_sub@data$AREA)
          writeRaster(conflicto_sub_raster, filename = "conflicto_sub_raster_1km.tif", format = "GTiff", progress = "text", overwrite = T)  
          plot(conflicto_sub_raster)
          
        # Conflictos por Sobretilización
        # -------------------------------
          conflicto_sobre <- readOGR(dsn = "IGAC", layer = "Conflictos_Sobreutlizacion") %>%
            spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          plot(conflicto_sobre)
          plot(linea_negra, add = T, border = "red")
        
          # Rasterizando la capa de conflicto por Sobreutilización
          conflicto_sobre_raster <- rasterize(conflicto_sobre, rasters_lights[[1]], conflicto_sobre@data$AREA)
          writeRaster(conflicto_sobre_raster, filename = "conflicto_sobre_raster_1km.tif", format = "GTiff", progress = "text", overwrite = T)  
          plot(conflicto_sobre_raster)
          
  # ------------------------
  # 5.5. AFECTACIÓN RELATIVA
  # ------------------------
      # Leyendo los polígonos del ìndice de Afectación Relativa (IRA)
        setwd("~/Dropbox/Linea_Negra_R/Data/")
      
      # Ìndice de Afectación Relativa (IRA) - TODOS
      # -------------------------------------------
        ira <- readOGR(dsn = "SIAC_ANLA", layer = "IndiceRelativoAfectacion_IRA_2010") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(ira)
        plot(linea_negra, add = T, border = "red")
      
        # Rasterizando la capa de IRA
        ira_raster <- rasterize(ira, rasters_lights[[1]], ira@data$afectacion)
        writeRaster(ira_raster, filename = "ira_raster_1km.tif", format = "GTiff", progress = "text", overwrite = T)  
        plot(ira_raster)   
      
      # Ìndice de Afectación Relativa (IRA) - ALTO/MUY ALTO
      # ---------------------------------------------------
        ira_alto_muyalto <- readOGR(dsn = "SIAC_ANLA", layer = "IndiceRelativoAfectacion_IRA_2010_Alto_MuyAlto") %>%
          spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        plot(ira_alto_muyalto)
        plot(linea_negra, add = T, border = "red")
        
        # Rasterizando la capa de IRA
        ira_alto_muyalto_raster <- rasterize(ira_alto_muyalto, rasters_lights[[1]])
        writeRaster(ira_alto_muyalto_raster, filename = "ira_alto_muyalto_raster_1km.tif", format = "GTiff", progress = "text", overwrite = T)  
        plot(ira_alto_muyalto_raster)   
        
#############################################################################################  
# 6. VARIABLES DE CONTROL: CÁLCULO DE LAS ÁREAS DE AFECTACIÓN POR LOS DIFERENTES FENÓMENOS
#############################################################################################
    
  # SOIL QUALITY
  # ------------
    setwd("~/Dropbox/Linea_Negra_R/Data/FAO/Soil_Quality")
    list_raster <- list.files()
    rasters <- lapply(list_raster, raster)
    rasters_extent <- extent(rasters[[1]])
    list(rasters)
    
    rasters_soil <- processing_rasters(rasters, rasters_extent, colombia) %>%
      raster::resample(rasters_lights, method = "ngb")
    plot(rasters_soil)
            
   # ELEVATION
   # ---------
    setwd("~/Dropbox/Linea_Negra_R/Data/USGS")
    rasters_extent <- extent(rasters_lights)
    elevation <- raster("altura_mean_30arc.tif") %>%
      crop(colombia) %>%
      setExtent(rasters_extent) %>%
      mask(colombia)
    plot(elevation)
    
    slope <- terrain(elevation, opt = "slope")
    aspect <- terrain(elevation, opt = "aspect")
    hills <- hillShade(slope, aspect, angle = 40, 0)
    roughness <- terrain(elevation, opt = "roughness")
    flowdir <- terrain(elevation, opt = "flowdir")
    names(hills) <- "hill"
    #Plot
    plot(slope)
    plot(aspect)
    plot(hills)
    plot(roughness)
    plot(flowdir)
    
    
        
      
    
      
  dev.off() # Clear all plots
  savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")
  save.image(file = "calculo_areas.RData")   
  load("calculo_areas.RData") 
    




