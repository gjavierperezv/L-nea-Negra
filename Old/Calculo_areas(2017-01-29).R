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
     
     # Create a dataframe and display default rownames
       (buffer_continental.df <- data.frame(ID=1:length(buffer_continental))) 
       rownames(buffer_continental.df)
      
        # Extract polygon ID's
       (buffer_continental_id <- sapply(slot(buffer_continental, "polygons"), function(x) slot(x, "ID")))
       
       # Create dataframe with correct rownames
       (buffer_continental.df <- data.frame( ID=1:length(buffer_continental), row.names = buffer_continental_id))    
       
       # Try coersion again and check class
       buffer_continental <- SpatialPolygonsDataFrame(buffer_continental, buffer_continental.df)
       class(buffer_continental) 
       slotNames(buffer_continental)
       
       # Now we can add a column
       buffer_continental@data$AREA_KM2_buffercontin <- gArea(buffer_continental, byid = T)/1e6
       sum(buffer_continental@data$AREA_KM2_buffercontin)
       
       
       
       
       
       
       ########################
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
       ####################
       
    buffer_continental2 <- DataSlot(croquis_colombia)
    slotNames(buffer_continental2) 
    
    
    
    
    
    
  # 4.2. PROTECCIÓN SENCILLA (LN - PNN - RESGUARDOS)
  # -----------------------------------------------
      # Función que crea un hueco en un polígono con la forma de otro
      # Refference: http://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe
      
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
      
      # Se crea la Unión entre PNN y Resguardos
      #resg_union_pnn <- gUnion(resguardos_ln_proj, pnn_ln, byid = T)
      #plot(resg_union_pnn, col = "blue", border = "red")
      
      # Se aplica la función para crear los huecos en el polígono de la Línea Negra
        temp15 <- AddHoleToPolygon(linea_negra_proj[1,], resguardos_ln_proj[15,]) 
        plot(temp15, col = "blue", border = "red")
        temp15@data$AREA_KM2_t15 <- gArea(temp15, byid = T)/1e6
        sum(temp15@data$AREA_KM2_temp15)
       
        temp16 <- AddHoleToPolygon(temp15[1,], resguardos_ln_proj[16,]) 
        plot(temp16, col = "blue", border = "red")
        temp16@data$AREA_KM2_temp16 <- gArea(temp16, byid = T)/1e6
        sum(temp16@data$AREA_KM2_temp16)
        
        
        
        
        ####################
        temp_1 <- AddHoleToPolygon(linea_negra_proj[1,], resguardos_ln_proj[1,])
        for (i in 2:length(resguardos_ln_proj)){
          for(j in 1:length(resguardos_ln_proj)){
            temp_i <- AddHoleToPolygon(temp_j[1,], resguardos_ln_proj[i,])
            i=i+1
            j=j-1
          }
        }
        #################
        plot(temp_2)
        dev.off() # Clear all plots
      
        
        
        
        
        
        # Posible comandos: unionSpatialPolygons, 
        lps <- getSpPPolygonsLabptSlots(resguardos_ln_proj)
        IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
        diss_resg <- unionSpatialPolygons(resguardos_ln_proj, IDOneBin)
        plot(diss_resg, col = "blue", border = "red")
        diss_resg
        slotNames(diss_resg)
            # Create a dataframe and display default rownames
            (diss_resg.df <- data.frame(ID=1:length(diss_resg))) 
            rownames(diss_resg.df)
            
            # Try to coerce to SpatialPolygonsDataFrame (will throw error)
            #resg_menos_pnn <- SpatialPolygonsDataFrame(resg_menos_pnn, resg_menos_pnn.df) 
            
            # Extract polygon ID's
            (diss_resg_id <- sapply(slot(diss_resg, "polygons"), function(x) slot(x, "ID")))
            
            # Create dataframe with correct rownames
            (diss_resg.df <- data.frame( ID=1:length(diss_resg), row.names = diss_resg_id))    
            
            # Try coersion again and check class
            diss_resg <- SpatialPolygonsDataFrame(diss_resg, diss_resg.df)
            class(diss_resg) 
            slotNames(diss_resg)
            
            # Now we can add a column
            diss_resg@data$AREA_KM2_dissresg <- gArea(diss_resg, byid = T)/1e6
            sum(diss_resg@data$AREA_KM2_dissresg)
  
        plot(diss_resg)
        diss_resg_lnhole <- AddHoleToPolygon(linea_negra_proj[1,], diss_resg[1,]) 
        diss_resg_lnhole@data$AREA_KM2_dissresglnhole <- gArea(diss_resg_lnhole, byid = T)/1e6
        sum(diss_resg_lnhole@data$AREA_KM2_dissresglnhole)
        
        plot(diss_resg_lnhole, col = "blue", border = "red")
        
 
      
  # 4.3. DOBLE PROTECCIÓN (LINEA NEGRA+RESGUARDOS = RESGUARDOS - PNN)
  # -----------------------------------------------------------------
  # Opción (polígono correcto pero sin el slot data)
  resg_menos_pnn <- gDifference(resguardos_ln_proj, pnn_ln)
  plot(linea_negra_proj)
  plot(resg_menos_pnn, add = T,  border = "red")
  
    # Construyendo el slot data
    # -------------------------
      # Create a dataframe and display default rownames
      (res_menos_pnn.df <- data.frame(ID=1:length(resg_menos_pnn))) 
      rownames(res_menos_pnn.df)
      
      # Try to coerce to SpatialPolygonsDataFrame (will throw error)
      #resg_menos_pnn <- SpatialPolygonsDataFrame(resg_menos_pnn, resg_menos_pnn.df) 
      
      # Extract polygon ID's
      (resg_menos_pnn_id <- sapply(slot(resg_menos_pnn, "polygons"), function(x) slot(x, "ID")))
      
      # Create dataframe with correct rownames
      (resg_menos_pnn.df <- data.frame( ID=1:length(resg_menos_pnn), row.names = resg_menos_pnn_id))    
      
      # Try coersion again and check class
      resg_menos_pnn <- SpatialPolygonsDataFrame(resg_menos_pnn, resg_menos_pnn.df)
      class(resg_menos_pnn) 
      slotNames(resg_menos_pnn)
      
      # Now we can add a column
      resg_menos_pnn@data$AREA_KM2_resmenospnn <- gArea(resg_menos_pnn, byid = T)/1e6
      sum(resg_menos_pnn@data$AREA_KM2_resmenospnn)
      
    plot(linea_negra_proj)
    plot(resguardos_ln_proj, add = T, border = "green")
    plot(resg_menos_pnn, add = T, border = "blue")
    plot(pnn_ln, add = T, border = "red")
    
  # 4.4. DOBLE PROTECCIÓN (LINEA NEGRA+PNN = PNN - RESGUARDOS)
  # -----------------------------------------------------------------
    pnn_menos_resg <- gDifference(pnn_ln, resguardos_ln_proj)
    plot(linea_negra_proj)
    plot(resguardos_ln_proj, add = T, border = "green")
    plot(pnn_ln, add = T, border = "red")
    plot(pnn_menos_resg, add = T, border = "blue")
    
      # Construyendo el slot data
      # Create a dataframe and display default rownames
      (pnn_menos_resg.df <- data.frame(ID=1:length(pnn_menos_resg))) 
      rownames(pnn_menos_resg.df)
      
      # Try to coerce to SpatialPolygonsDataFrame (will throw error)
      #resg_menos_pnn <- SpatialPolygonsDataFrame(resg_menos_pnn, resg_menos_pnn.df) 
      
      # Extract polygon ID's
      (pnn_menos_resg_id <- sapply(slot(pnn_menos_resg, "polygons"), function(x) slot(x, "ID")))
      
      # Create dataframe with correct rownames
      (pnn_menos_resg.df <- data.frame( ID=1:length(pnn_menos_resg), row.names = pnn_menos_resg_id))    
      
      # Try coersion again and check class
      pnn_menos_resg <- SpatialPolygonsDataFrame(pnn_menos_resg, pnn_menos_resg.df)
      class(pnn_menos_resg) 
      slotNames(pnn_menos_resg)
      
      # Now we can add a column
      pnn_menos_resg@data$AREA_KM2_pnnmenosresg <- gArea(pnn_menos_resg, byid = T)/1e6
      sum(pnn_menos_resg@data$AREA_KM2_pnnmenosresg)

  # 4.5. TRIPLE PROTECCIÓN (INTERSECCIÓN ENTRE PNN Y RESGUARDOS)
  # ------------------------------------------------------------
  resg_in_pnn <- raster::intersect(resguardos_ln_proj, pnn_ln)
  resg_in_pnn@data$AREA_KM2_resinpnn <- gArea(resg_in_pnn, byid = T)/1e6
  dev.off() # Clear all plots
  plot(linea_negra_proj)
  #plot(resguardos_ln_proj, add = T)
  #plot(pnn_ln, add = T)
  plot(resg_in_pnn, add = T, border = "red")
  sum(resg_in_pnn@data$AREA_KM2_resinpnn)

###################################################################  
# 5. CÁLCULO E LAS ÁREAS DE AFECTACIÓN POR LOS DIFERENTES FENÓMENOS
###################################################################

  # 5.1. DEFORESTACIÓN
  # ------------------
  
    # 5.1.1. SIN PROTECCIÓN (Fuera de la Línea Negra (dentro de un Buffer de 50 km2))
    # -------------------------------------------------------------------------------
    
    
    # 5.1.2. Dentro de la Línea Negra
    # -------------------------------
    
    
    # 5.1.3. Dentro de los PNN
    # ------------------------
  
  
    # 6.4. Dentro de los Resguardos
    # -----------------------------
  
  
    # 6.5. Dentro de los Resguardos que están al interior de los PNN
    # --------------------------------------------------------------
  
  
  
    # 6.6. Dentro de los PNN que están al interior de los Resguardos
    # --------------------------------------------------------------

  # 5.2. POBLAMIENTO (LUCES NOCTURNAS)
  # ----------------------------------
  
  
  # 5.3. VÍAS
  # ---------

  
  
  # 5.4. CONFLICTO DE USO DEL SUELO
  # -------------------------------
  
  
  # 5.5. AFECTACIÓN RELATIVA
  # ------------------------
  
  
  
  dev.off() # Clear all plots
  savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")
    
    
    




