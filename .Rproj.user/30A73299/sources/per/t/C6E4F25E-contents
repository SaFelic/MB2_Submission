##############################################################################################
### Script for classification of satellite data with subsequent Land Cover Change Assessment
### and further analysis of MODIS NDVI Raster Time Series (RTS) 
### including statistical RTS analysis and breakpoint detection using bfastmonitor
##############################################################################################

##############################################################################################
### Main web sources:
### https://geocompr.github.io/post/2019/ggplot2-inset-maps/
### https://www.earthdatascience.org/courses/earth-analytics/lidar-raster-data-r/classify-raster/
### https://github.com/wegmann/R_scripts/blob/master/supervised_classification.R
##############################################################################################

##############################################################################################
### Written by: Sarah Schneider
##############################################################################################

# function checking if package exists and if necessary installation of missing packages (adapted from Wegmann, s. web soruces)
loadandinstall <- function(mypkg) {if (!is.element(mypkg, installed.packages()[,1])){install.packages(mypkg)}; library(mypkg, character.only=TRUE)  }

loadandinstall("raster")
loadandinstall("sf")
loadandinstall("randomForest")
loadandinstall("raster")
loadandinstall("rgdal")
loadandinstall("maptools")
loadandinstall("RStoolbox")
loadandinstall("raster")
loadandinstall("cowplot")
loadandinstall("rasterVis")
loadandinstall("dplyr")
loadandinstall("ggplot2")
loadandinstall("lubridate")
loadandinstall("gridExtra")

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 1:                             ###
###                      LAND COVER CLASSIFICATION                      ###
###                                                                     ###
###########################################################################
###########################################################################

#setwd("path/to/downloaded/data/folder")

# importing Raster and Training Data
ras <- list.files(path=".", pattern="clipped", full.names=TRUE)
ras <- lapply(ras, brick)

vec <- list.files(path=".", pattern="gpkg$", full.names=TRUE)
vec <- lapply(vec, readOGR)

# Random Forest Classification 

class_2001 <- superClass(ras[[1]], vec[[1]],
                         model = "rf", trainPartition = 0.7, responseCol = "class")

class_2020 <- superClass(ras[[2]], vec[[2]],
                         model = "rf", trainPartition = 0.7, responseCol = "class")

# Validation results
class_2001$validation$performance # Accuracy Assessment of 2001 Classification

class_2020$validation$performance # Accuracy Assessment of 2020 Classification

# Validation results visualized
cm_2001 <- as.data.frame(class_2001$validation$performance$table)
ggplot(data=cm_2001, aes(Prediction, Reference)) + geom_point(aes(size = Freq))

cm_2020 <- as.data.frame(class_2020$validation$performance$table)
ggplot(data=cm_2020, aes(Prediction, Reference)) + geom_point(aes(size = Freq))

# Classification maps
p_2001 <- ggplot() + ggR(class_2001$map, geom_raster = T, ggLayer = T) +
          scale_fill_manual(values = c("blue", "darkgreen", "grey", "red", "white"), name = "Land Cover", na.value = NA)+
          coord_sf(crs=st_crs(ras[[1]]))+
          labs(title = "Land Cover Classification Brazil, 2001",
               x = "Longitude", 
               y = "Latitude")

p_2020 <- ggplot() + ggR(class_2020$map, geom_raster = T, ggLayer = T) +
          scale_fill_manual(values = c("blue", "darkgreen", "grey", "red", "white"), name = "Land Cover", na.value = NA)+
          coord_sf(crs=st_crs(ras[[1]]))+
          labs(title = "Land Cover Classification Brazil, 2020",
               x = "Longitude", 
               y = "Latitude")

pdf("Classification_Maps_2001_2020.pdf", width = 20, height = 8) 
grid.arrange(p_2001, p_2020, ncol=2)
dev.off()

############################################################################
############################################################################
###                                                                      ###
###                              SECTION 2:                              ###
###               CHANGE DETECTION WITH POSTCLASSIFICATION               ###
###                                                                      ###
############################################################################
############################################################################

# Change Map calculation
class_2001_rast <- class_2001$map*10 
class_2020_rast <- class_2020$map

class_2001_2020 <- class_2001_rast+class_2020_rast

# Visualization
ggplot() + ggR(class_2001_2020, geom_raster = T, ggLayer = T, forceCat = T) +
  coord_sf(crs=st_crs(ras[[1]]))

# Reclassifying raster
reclass_df <- c(11, 1, # all Unchanged classes are assigned to class 1
                22, 1,
                33, 1,
                44, 1,
                21, 2, # all Vegetation loss classes are assigned to class 2
                23, 2,
                24, 2,
                12, 3, # all Vegetation gain classes are assigned to class 3
                32, 3,
                42, 3,
                13, 4, # all Rest classes are assigned to class 4
                14, 4,
                31, 4, 
                34, 4,
                41, 4,
                43, 4)

reclass_m <- matrix(reclass_df,
                    ncol = 2,
                    byrow = TRUE)

reclass_2001_2020 <- reclassify(x=class_2001_2020, rcl=reclass_m, include.lowest=TRUE)

# displaying the distribution of pixels per Land Cover Class
barplot(reclass_2001_2020,
        main = "Number of pixels per class") 

# renaming numerical classes 
reclass_2001_2020_df <- as.data.frame(reclass_2001_2020, xy = TRUE, na.rm = TRUE) 
reclass_2001_2020_df['names'] <- NA # new column for names of LC Change classes

reclass_2001_2020_df <- reclass_2001_2020_df %>% # renaming numerical classes in new "names" column
  mutate(names = case_when(
    layer == 1 ~ "Unchanged",
    layer == 2 ~ "Vegetation Loss",
    layer == 3 ~ "Vegetation Growth",
    layer == 4 ~ "Rest"
  ))

# first Visualization of of reclassified Land Cover Change Map
ggplot() + 
  geom_raster(data = reclass_2001_2020_df,
              aes(x = x, y = y, fill = names)) +
  labs(title = "Land Cover Change Map in Brazil",
       subtitle = "2001 - 2020",
       x = "Longitude",
       y = "Latitude",
       fill = "Change Class") +
  coord_sf(crs=st_crs(ras[[1]]))+
  theme_minimal()

# calculating changed area 
change_area <- freq(reclass_2001_2020$layer, useNA = "no")
change_area[,"count"] * 30^2 * 1e-06

ca <- as.data.frame(change_area)
ca["nms"] <- NA

ca <- ca %>% # renaming numerical classes in new "names" column for subsequent plotting
  mutate(nms = case_when(
    value == 1 ~ "Unchanged",
    value == 2 ~ "Vegetation Loss",
    value == 3 ~ "Vegetation Growth",
    value == 4 ~ "Rest"
  ))

# plotting calculated Change Area 
coul <- RColorBrewer::brewer.pal(4, "Pastel1") 

ggplot(data=ca, aes(y=nms, x=(count/1000))) +
  geom_bar(stat="identity", fill=coul)+
  geom_text(aes(label=(round(count/1000))),hjust = -0.25)+
  labs(y="Change Class",
       x="Area (km²)")+
  ggtitle("Area of Land Cover Change in km²")


# plotting final reclassified change map with inset map
Area <- read_sf("Area.shp") # loading area for inset map bounding box
ins_map <- read_sf("bra_adm_ibge_2020_shp/bra_admbnda_adm1_ibge_2020.shp") #loading Brazil shapefile
ins_map_bb = st_as_sfc(st_bbox(Area)) # bounding box of Area

ins_map <- st_set_crs(ins_map, "EPSG:4326") # setting crs 
ins_map <- st_transform(ins_map, crs = "EPSG:4326") # and reprojecting loaded shapefiles
ins_map_bb <- st_transform(ins_map_bb, crs = "EPSG:4326")

# plot of inset map
inset_map = ggplot() + 
  geom_sf(data = ins_map, fill = "white") + 
  geom_sf(data = ins_map_bb, fill = NA, color = "red", size = 0.5) +
  theme_void()+
  coord_sf(crs=st_crs(ras[[1]]))
inset_map

# plot of Change Map
reclass_map <- ggplot(reclass_2001_2020_df) +  
  geom_raster(aes(x=x, y=y, fill= names))+ 
  scale_fill_manual(values = c("bisque4", "lavenderblush2", "chartreuse4", "coral3"))+
  labs(title = "2001 - 2020 Land Cover Change in Novo Progresso, Brazil",
       fill="Change Class", 
       x = "Longitude", 
       y = "Latitude")+
  coord_sf(crs=st_crs(ras[[1]]))
reclass_map

# plot of both maps
final_map = ggdraw() +
  draw_plot(reclass_map) +
  draw_plot(inset_map, x = 0.65, y = 0.65, width = 0.3, height = 0.3)
ggsave("LCC_Map.pdf")
