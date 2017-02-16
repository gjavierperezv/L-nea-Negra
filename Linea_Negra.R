library(raster)
library(rgdal)
library(rasterVis)
library(sp)
library(ggplot2)
library(maptools)
library(lattice)
library(plyr)
library(grid)
library(dplyr)
library(stringr)
library(gtools)
library(magrittr)
library(rgeos)


########################################################################################################
#######################################    GETTING THE DATA READY     ##################################
########################################################################################################

# HOLE FREE FUNCTION
# ------------------
  hole_free <- function(x){
    BCp <- slot(x, "polygons")
    holes <- lapply(BCp, function(y) sapply(slot(y, "Polygons"), slot, "hole")) 
    res <- lapply(1:length(BCp), function(i) slot(BCp[[i]], "Polygons")[!holes[[i]]]) 
    IDs <- row.names(x)
    SpatialPolygons(lapply(1:length(res), function(i) Polygons(res[[i]], ID=IDs[i])), proj4string = CRS(proj4string(x))) 
  }

# PROCESING RASTERS FUNCTION: crop, mask, make calculations and extract raster data
# ---------------------------------------------------------------------------------
  processing_rasters <- function(layer.list, ext, shape){
    layer.list %>%
      lapply(setExtent, ext) %>%
      lapply(crop, shape) %>%
      stack() %>% 
      mask(shape)
  }
  
  
# 1. PREPARING POLYGONS IN THE AREA OF INFLUENCE
# -----------------------------------------------

    # 1.1. DEPARTAMENTS POLYGON: LA GUAJIRA, MAGDALENA, CESAR
    # -------------------------------------------------------
    setwd("~/Dropbox/Linea_Negra_R/Data/Colombia_Deptos/")
    deptos_ln <- readOGR(".", "Guaj_Magd_Cesa_Deptos")  #The dsn parameter is the directory, but if it is the current working directory, we can use the "." option
    deptos_ln_proj <- spTransform(deptos_ln, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    #plot(deptos_ln_proj)
    
    deptos_ln_croquis <- gUnaryUnion(deptos_ln_proj)
    #plot(deptos_ln_croquis)
  
    # 1.2. MUNICIPALITIES POLYGON: LA GUAJIRA, MAGDALENA, CESAR
    # ---------------------------------------------------------
    setwd("~/Dropbox/Linea_Negra_R/Data/Colombia_Municipios/")
    mpios_ln <- readOGR(".", "Guaj_Magd_Cesa_Mpios")    #The dsn parameter is the directory, but if it is the current working directory, we can use the "." option
    mpios_ln_proj <- spTransform(mpios_ln, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    #plot(mpios_ln_proj)  

  # 1.3. LINEA NEGRA POLYGON
  # ------------------------
    setwd("~/Dropbox/Linea_Negra_R/Data/Linea_Negra/")
    linea_negra <- readOGR(dsn = ".", layer = "Linea_Negra_Polygon") %>%
      spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    linea_negra_proj <- linea_negra %>% spTransform(CRS=CRS("+init=epsg:3857"))
    # plot(linea_negra_proj)
  
  # 1.4. COLOMBIA
  # -------------
    setwd("~/Dropbox/Linea_Negra_R/Data/Colombia_Deptos/")
    colombia <- readOGR(dsn = ".", layer = "Colombia_Deptos") %>%
      spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
      #spTransform(CRS=CRS("+init=epsg:3857")) %>%
      .[!(.@data$DPTO_CCDGO == 88) , ] %>%
      gUnaryUnion(.) %>%
      hole_free()
     #plot(colombia)

# 2. NIGHTLIGHT DATA
#-------------------
# Open .tif files as a raster (the raster package allow to read these files in the disk and not in the memory, this improves the efficiency of functions in R)

setwd("~/Dropbox/Linea_Negra_R/Data/NOAA/AVER_AND_CLOUDFREE/Web_Stable_Average_Visible/")
list_raster <- list.files() %>%
  lapply(raster)
list(list_raster)
rasters_extent <- extent(list_raster[[1]]) #We need to put all rasters into the same extent (all have the same resolution). [] extracts a list, [[]] extracts elements within the list
rasters_lights <- processing_rasters(list_raster, rasters_extent, colombia)
plot(rasters_lights)

# 4. DEFORESTATION DATA
#----------------------
  # 4.1. COLOMBIA
  # -------------

  # 4.1.1. Open rasters using raster package and crop to Colombias's shape 
  # ----------------------------------------------------------------------
    # Oppening Colombia's Shape File  (NOTE:OPEN THIS FILE BEFORE DE RASTERS)
      # setwd("~/Dropbox/Linea_Negra_R/Data/Colombia_Municipios")
      # colombia_mpios_proj <- 
      #   readOGR(dsn = ".", layer="pladiv_polygon") %>%  #The dsn parameter is the directory, but if it is the current working directory, we can use the "." option
      #   spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      #plot(colombia_mpios_proj)
      
      # Oppening Deforestation raster files
      setwd("~/Dropbox/Linea_Negra_R/Data/WDPA")
      files <- list.files() %>%
        str_detect("lossyear")  #List of files with the word lossyear
      
      # Applying the crop function to the raster files
      loss_year <- lapply(list.files()[files], raster) %>%
        lapply(crop, colombia) #lapply returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X.
      head(loss_year)
      
  # 4.1.2. Mosaic images to create a RasterLayer
  # --------------------------------------------
    # Raster file at 30 meters
    # -------------------------------
      loss_year$fun <- mean
      loss_year$na.rm <- TRUE
      system.time(loss_year <- do.call(mosaic, loss_year)) # (15 min)  do.call constructs and executes a function call from a name or a function and a list of arguments to be passed to it.
      
      writeRaster(loss_year, filename = "loss_year_mosaic_30m.tif", format = "GTiff", progress = "text", overwrite = T)  
      removeTmpFiles(0.1)    
    
      # Map
      setwd("~/Dropbox/Linea_Negra_R/Data/WDPA")
      loss_year_mosaic_30m <- raster("loss_year_mosaic_30m.tif")
      #plot(loss_year_mosaic_30m)
     
    # Raster file at 1 kilometer
    # --------------------------
    # Refference: https://github.com/ivanhigueram/deforestacion/blob/master/Extracting/raster_loss.R
      
      #Layerize create big files in .grd binary format which does not have any compression. For that reason we change 
      # the tmpdir directory to a new one with more space (1 Tb LaCie EHD). 
      
      #dir.create("Temp")
      #rasterOptions(tmpdir = "~/Dropbox/Linea_Negra_R/Data/WDPA/Temp")
        system.time(
        loss_year_brick <- layerize(loss_year_mosaic_30m, filename = "loss_year_brick.tif",  # For layerize: https://www.rdocumentation.org/packages/raster/versions/2.5-8/topics/layerize
                                    overwrite = T,
                                    format = "GTiff",
                                    options = "INTERLEAVE=BAND", 
                                    progress = "text"))
       #plot(loss_year_brick)
     
      # Now, we want to create a grid of 1 km2 (using the night light data as reference grid)
      setwd("~/Dropbox/Linea_Negra_R/Data/WDPA")
      loss_year_brick <-  brick("loss_year_brick.tif") 
      fact <- round(dim(loss_year_brick)[1:2] / dim(rasters_lights[[1]])[1:2]) 
      
      agg <- aggregate(loss_year_brick, fact, filename = "loss_year_brick_agg.tif",
                       overwrite = T,
                       format = "GTiff",
                       options = "INTERLEAVE = BAND",
                       progress = "text" ) #Proportion of deforestation in the 1 km grid 
      
       # Lectura de archivo ya existente
      #agg <- blick("loss_year_brick_agg.tif")
      
      agg_mask <- mask(agg, colombia)
      rasters_lights <- setExtent(rasters_lights, agg_mask)
      
      res <- resample(agg_mask, rasters_lights, filename = "loss_year_brick_1km.tif",  # Resample: transfers values between non matching Raster* objects (in terms of origin and resolution).
                      format = "GTiff",
                      options = "INTERLEAVE=BAND", 
                      progress = "text", 
                      overwrite = T)
      plot(res)
      # Lectura de archivo ya existente
      res <- brick("loss_year_brick_1km.tif")
      plot(res)
      #Verification that all cells are the same
      identical(coordinates(res), coordinates(rasters_lights))
      

  # 4.2. "LINEA NEGRA": Area of influence (La Guajira, Magdalena, Cesar)
  # --------------------------------------------------------------------
      #res_ln <- mask(res, mpios_ln_proj, filename = "loss_year_brick_ln_1k.tif",
      #             overwrite = T)
      #plot(res_ln)      
      
# 5. SOIL QUALITY
#----------------  
 
  # 5.1. SOIL DATA FROM FAO
  # -----------------------
    #Download data
    #setwd("/Volumes/LaCie/Datos") 
    #url <- "http://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/"
    #files <- str_c("sq", c(1:7), ".asc")
    #filenames_list <- as.list(files)
    
    #l_ply(filenames_list, download,
    #      baseurl = "http://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/",
    #      folder = "soil_quality"
    #)

  # 5.2. PROCESSING DATA
  # --------------------
  setwd("~/Dropbox/Linea_Negra_R/Data/FAO/Soil_Quality")
  list_raster <- list.files()
  rasters <- lapply(list_raster, raster)
  rasters_extent <- extent(rasters[[1]])
  list(rasters)
  
  #setwd("~/Dropbox/Linea_Negra_R/Data/Colombia_Municipios/")
  #mpios_ln <- readOGR(".", "Guaj_Magd_Cesa_Mpios")    #The dsn parameter is the directory, but if it is the current working directory, we can use the "." option
  #mpios_ln_proj <- spTransform(mpios_ln, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  #plot(mpios_ln_proj)
  
  rasters_soil <- processing_rasters(rasters, rasters_extent, colombia) %>%
    raster::resample(rasters_lights, method = "ngb")
  # plot(rasters_soil)
  
# 6. ELEVATION DATA FROM  USGS
# ----------------------------
  
  # 6.1. DOWNLOADING DATA
  # ---------------------
  #download.file(     
  #  url = "http://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Global_tiles_GMTED/300darcsec/mea/W090/10S090W_20101117_gmted_mea300.tif" ,
  #  destfile = "altura_mean_30arc.tif", mode="wb")
  
  # 6.2. PROCESSING DATA
  # --------------------
  setwd("~/Dropbox/Linea_Negra_R/Data/USGS")
  rasters_extent <- extent(rasters_lights)
  elevation <- raster("altura_mean_30arc.tif") %>%
    crop(colombia) %>%
    setExtent(rasters_extent) %>%
    mask(colombia)
  plot(elevation)

# 7. SLOPE AND ASPECTS
# --------------------
  slope <- terrain(elevation, opt = "slope")
  aspect <- terrain(elevation, opt = "aspect")
  hills <- hillShade(slope, aspect, angle = 40, 0)
  roughness <- terrain(elevation, opt = "roughness")
  flowdir <- terrain(elevation, opt = "flowdir")
  names(hills) <- "hill"
  #Plots
  plot(slope)
  plot(aspect)
  plot(hills)
  plot(roughness)
  plot(flowdir)
  
# 9. RASTER CLUMPS
# ----------------
  # Refference: https://github.com/ivanhigueram/deforestacion/blob/master/Extracting/raster_clumps.R
  
  # 9.1 Get deforestation raster for reference
  # -----------------------------------------
    setwd("~/Dropbox/Linea_Negra_R/Data/WDPA")
    res_ln <- brick("loss_year_brick_ln_1k.tif")
    #plot(res_ln)
  
  # 9.2 GET ADMINISTATIVE GIS DATA
  # -------------------------------
    # setwd("~/Dropbox/Linea_Negra_R/Data/Colombia_Municipios/")
    # mpios_ln_proj <- 
    #   readOGR(".", "Guaj_Magd_Cesa_Mpios")  %>%  #The dsn parameter is the directory, but if it is the current working directory, we can use the "." option
    #   spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    #plot(mpios_ln_proj)  
 
  # 9.3. GETTING NIGHTLIGHT DATA
  # ----------------------------
   
    # Open .tif files as a raster (the raster package allow to read these files in the disk and not in the memory, this improves the efficiency of functions in R)
    setwd("~/Dropbox/Linea_Negra_R/Data/NOAA/AVER_AND_CLOUDFREE/Web_Stable_Average_Visible/")
    list_raster <- list.files() %>%
      lapply(raster)
    list(list_raster)
    rasters_extent <- extent(list_raster[[1]]) #We need to put all rasters into the same extent (all have the same resolution). [] extracts a list, [[]] extracts elements within the list
    rasters_lights <- processing_rasters(list_raster, rasters_extent, colombia)
    list(rasters_lights)
    #plot(rasters_lights)
    
  # 9.4. CLUMPS TO IDENTIFY CITIES (QUEEN) FOR 1993 (THIS YEAR was chosen BECAUSE THE LEGISLATION STARTED IN AUGUST 1995)
  # ---------------------------------------------------------------------------------------------------------------------
    clump_lights <- clump(rasters_lights[[2]], directions = 8) %>%  # Detect clumps (patches) of connected cells. Each clump gets a unique ID. NA and zero are used as background values (i.e. these values are used to separate clumps). You can use queen's or rook's case, using the directions argument.
      resample(., res[[1]])     # Resample transfers values between non matching Raster* objects (in terms of origin and resolution)
    clump_lights[clump_lights > 1] <- 1
    list(clump_lights)
    #plot(clump_lights)
    
  # 9.5. CLUMPS TO POLYGONS
  # -----------------------
    p1 <- rasterToPolygons(clump_lights, dissolve = T)  #Cells with NA are not converted. Dissolve: If TRUE, polygons with the same attribute value will be dissolved into multi-polygon regions. This option requires the rgeos package.
    setwd("~/Dropbox/Linea_Negra_R/Data")
    writeOGR(p1, "polygon_clump_layer_1993.shp",
             layer = "clumps", 
             driver = "ESRI Shapefile", 
             overwrite_layer = T)
    clumps <- readOGR(".", "polygon_clump_layer_1993")
    #plot(clumps)
    #list(p1)
    #plot(p1)   
  
# 10. PROCESS DATA
# -----------------
    setwd("~/Dropbox/Linea_Negra_R/Data/USGS")
    rasters_extent <- extent(rasters_lights)
    elevation <- raster("altura_mean_30arc.tif")  %>%
      crop(colombia) %>%
      setExtent(rasters_extent) %>%
      mask(colombia) %>%
      resample(., res[[1]])
    #list(elevation)  
    #plot(elevation)
    
   
  # 10.1. SLOPE AND ASPECTS
  # ------------------------
  slope <- terrain(elevation, opt = "slope")  #Compute slope, aspect and other terrain characteristics from a raster with elevation data. The elevation data should be in map units (typically meter) for projected (planar) raster data. 
  roughness  <- terrain(elevation, opt = "roughness")
  tri <- terrain(elevation, opt = "TRI")  #Terrain Ruggedness Index (TRI)
  #plot(slope)
  #plot(roughness)  
  #plot(tri)   
  
  # 10.2. GET CLIMATE
  # -----------------
  # FUENTE DE LOS DATOS: http://worldclim.org/current
  # The current file was previously procesed by https://github.com/ivanhigueram/
  
  # Precipitation
  # -------------
  setwd("~/Dropbox/Linea_Negra_R/Data/World_Climate/")
  # list_files <- list.files() %>%
  #   str_detect("prec")
  # 
  # prec <- lapply(list.files()[list_files], raster) %>%
  #   brick() %>%
  #   lapply(crop, colombia_municipios) %>%
  #   lapply(mask, colombia_municipios) %>%
  #   resample(res[[1]], filename = "prec_1km.tif",
  #            format = "GTiff",
  #            options = "INTERLEAVE=BAND", 
  #            progress = "text", overwrite = T)
  prec <- brick("prec_1km.tif")
  plot(prec)
  
   # 10.3. SOIL QUALITY
  # ------------------
  # NOTA: El archivo sq_1km_ln_tif.tif no lo guardo en la subcarpeta Soil_Quality dentro de FAO sino en FAO directamente.
  # La razón es que en un procedimiento anterior se hace un listado de los rasters en esa subcarpeta (sq1.asc a sq7.asc),
  # y si q_1km_ln_tif.tif se encunetra allí, lo tama tambien para definir el extent y genera confusión y el siguiente
  # error:  Error in compareRaster(x) : different extent. Por tanto se cambia a una ubicación más arriba, dentro de FAO.
  
  setwd("~/Dropbox/Linea_Negra_R/Data/FAO/Soil_Quality")
  list_files <- list.files()
  
  sq <- lapply(list_files, raster) %>%
    brick() %>%
    crop(colombia) %>%
    mask(colombia) %>%
    resample(res[[1]], filename = "sq_1km_ln_tif",
             format = "GTiff",
             options = "INTERLEAVE = BAND",
             progress = "text", overwrite = T)
  #list(sq)
  plot(sq)
  sq <- brick("sq_1km_ln_tif.tif")
  #plot(sq)

  # 10.4. ROADS
  # -----------
  #setwd("~/Dropbox/Linea_Negra_R/Data/IGAC")
  #roads_deptos_ln <- readOGR(dsn = ".", layer="Vias_Pavimentadas_Deptos_ln") %>%
    # spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) to mercator
  #spTransform(CRS=CRS("+init=epsg:3857"))
  #plot(roads_deptos_ln)
  
  #buffer_roads <- gBuffer(roads_deptos_ln, 50000, byid = T, id = row.names(roads_deptos_ln)) %>%
    #spTransform(CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  #plot(buffer_roads)
  # cells_roads <- cellFromPolygon(res[[1]], buffer_roads)
 
  #cells_roads <- unlist(cells_roads)
  #write.csv(cells_roads, "roads.csv")
  
  # 10.5. EXTRACT DATA TO data.frame
  # --------------------------------
  dataframes_extract <- list(elevation, slope, roughness, tri, clump_lights, sq, prec) %>%  # Falta inlcuir precipitacion prec
    lapply(as.data.frame, na.rm  = T) %>%
    lapply(function(x){
      x$ID <- row.names(x); x
    })
  
    # 10.5.1. MERGE
    # -------------
    merge_rasters_dataframes <- Reduce(function(...) merge(..., by="ID", all = T), dataframes_extract) %>%
      #mutate(roads_deptos_ln = ifelse(ID %in% unlist(cells_roads), 1, 0)) %>%
      mutate(clumps_1 = ifelse(is.na(clumps), 0, 1)) %>%
      mutate(prec = select(., starts_with("prec")) %>% rowMeans(na.rm = TRUE))
    
    # 10.5.2 Export .csv with clumps and ID
    # -------------------------------------
    setwd("~/Dropbox/Linea_Negra_R/Data/")
    write.csv(merge_rasters_dataframes, "geographic_covariates.csv")
  
  
  dev.off() # Clear all plots
  savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")








