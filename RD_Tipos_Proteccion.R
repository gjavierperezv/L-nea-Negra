library(plyr)
library(dplyr)
library(raster)
library(sp)
library(stringr)
library(rgdal)
library(rdrobust)
library(Hmisc)
library(stargazer)
library(data.table)

setwd("~/GitHub/linea_negra/")
source("Functions.R")

#Directories
# Javier
#data <- "~/Dropbox/Linea_Negra_R/Data/Dataframes/"
# Iván 
data <- "~/Dropbox/BANREP/Linea_Negra_R/Data/"
setwd(data)

#Load data from rasters
distance <- fread("Dataframes/distance_dataframe_protected_areas.csv") %>% plyr::rename(., c("layer" = "dist"))
defo <- fread("Dataframes/deforestation_dataframe.csv") %>% as.data.frame() %>% .[complete.cases(.[, c(3:17)]), ]
lights <- fread("Dataframes/rasters_lights_dataframe.csv") %>% as.data.frame() %>% .[complete.cases(.[, c(3:37)]), ]

#Aggregate deforestation (2001 - 2012)
defo$loss_sum <- rowSums(defo[, names(defo)[4:16]])
loss_sum <- dplyr::select(defo, c(ID, loss_sum))

#Merge data  (1: LN, 2: Resg, 3: PNN) 
defo_dist <- distance %>%
  merge(., loss_sum, by.x = "ID", by.y = "ID") %>%
  merge(., lights, by = "ID") %>%
    mutate(., loss_sum = loss_sum * 100) %>%
    mutate(., dist_disc = ifelse(treatment == 1, 1, -1) * dist) %>%
    mutate(., dist_disc = dist_disc / 1000)
    # merge(., cov, by = "ID", all.x = T) %>%
    # merge(., clump, by = "ID", all.x = T) %>%
    # mutate(clumps = ifelse(is.na(clumps), 0, 1))


#RD - Deforestation
defo_dist_buffer <- split(defo_dist, defo_dist$buffer_id) %>%
  lapply(., function(x){
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



#RD results to tables (df) and export to txt

#Deforestation
setwd(str_c(data, "Tables"))
df_optimal <- rd_to_df(defo_dist_buffer, split(defo_dist, defo_dist$buffer_id))
stargazer(df_optimal[2:dim(df_optimal)[1], ], type = "text", summary = FALSE, rownames = TRUE, out = "rd_type_deforestation.txt")

#Lights
df_optimal <- rd_to_df_2(lights_dist_buffer, split(defo_dist, defo_dist$buffer_id))
stargazer(df_optimal[2:dim(df_optimal)[1], ], type = "text", summary = FALSE, rownames = TRUE, out = "rd_type_lights.txt")

#Graphs!

defo_dist_all <- split(defo_dist, defo_dist$buffer_id) #Collapse all dataframes into one list and remove "all"

l <- lapply(defo_dist_all, function(x){
  mutate(x, bin = cut(x$dist_disc, breaks = c(-50:50), include.lowest = T)) %>%
    group_by(bin) %>%
    summarise(meanbin = mean(loss_sum), sdbin = sd(loss_sum)) %>%
    .[complete.cases(.),] %>%
    as.data.frame() %>%
    mutate(treatment = ifelse(as.numeric(row.names(.)) > 50, 1, 0), bins = row.names(.)) %>%
    mutate(bins = mapvalues(.$bins, from = c(1:100), to = c(-50:49)))
})



#Individual graphs for all territories (natural parks + territories)
setwd(str_c(data, "Graphs"))
mapply(function(x, type){
  g <- ggplot(x, aes(y = (meanbin), x = as.numeric(bins), colour = as.factor(treatment))) 
  g <- g + stat_smooth(method = "auto") 
  g <- g + geom_point(colour = "black", size = 1)
  g <- g + labs(x = "Distancia (Km.)", y = "Deforestación (Ha x Km2)")
  # g <- g + scale_x_continuous(limits = c(-20, 20))
  # g <- g + scale_y_continuous(limits = c(0, 0.3))
  g <- g + ggtitle(str_c("Discontinuidad\n", "para", type, sep = " "))
  g <- g + guides(colour = FALSE)
  g <- g + theme_bw()
  g
  ggsave(str_c("RDggplot_", type, "strategy2",".pdf"), width=30, height=20, units="cm")
}, x = l, type = c("LineaNegra", "Doble_LN_PNN", "Doble_LN_Resg", "Triple"))




# Saving history Commands' Backups
savehistory("~/Dropbox/Linea_Negra_R/Programs/Buckups.Rhistory")




