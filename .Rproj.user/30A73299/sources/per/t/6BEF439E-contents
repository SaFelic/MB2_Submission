### ##########################################################################################
### Script for Analysis of MODIS NDVI Raster Time Series (RTS) 
### including statistical RTS analysis and breakpoint detection using bfastmonitor
##############################################################################################

##############################################################################################
### Main web sources:
### https://cran.r-project.org/web/packages/MODISTools/vignettes/modistools-vignette.html
### http://www.loicdutrieux.net/bfastSpatial/
### https://github.com/khufkens/MODISTools_bfast_tutorial/blob/master/tutorial.md
### https://github.com/wegmann/R_scripts/blob/master/supervised_classification.R
##############################################################################################

##############################################################################################
### Written by: Sarah Schneider
##############################################################################################

# function checking if package exists and if necessary installation of missing packages (adapted from Wegmann, s. web soruces)
loadandinstall <- function(mypkg) {if (!is.element(mypkg, installed.packages()[,1])){install.packages(mypkg)}; library(mypkg, character.only=TRUE)  }

loadandinstall("MODISTools")
loadandinstall("raster")
loadandinstall("grDevices")
loadandinstall("rasterVis")
loadandinstall("viridis")
loadandinstall("bfast")
loadandinstall("bfastSpatial")

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 3:                             ###
###                   MODIS NDVI TIME SERIES ANALYSIS                   ###
###                                                                     ###
###########################################################################
###########################################################################

# listing products and bands to extract NDVI and Quality data
mod_products <- mt_products()
head(mod_products)

mod_bands <- mt_bands("MOD13Q1")
head(mod_bands)

# Downloading NDVI Data for Site with 20 km Buffer
brazil_NDVI <- mt_subset(product = "MOD13Q1",
                         lat = -7.03110502,
                         lon = -55.66010780,
                         band = "250m_16_days_NDVI",
                         start = "2001-01-01",
                         end = "2020-12-31",
                         km_lr = 20,
                         km_ab = 20,
                         site_name = "Brazil", 
                         progress = T)

# Downloading the QA Data for Site with 20 km Buffer
brazil_QA <- mt_subset(product = "MOD13Q1",
                       lat = -7.03110502,
                       lon = -55.66010780,
                       band = "250m_16_days_pixel_reliability",
                       start = "2001-01-01",
                       end = "2020-12-31",
                       km_lr = 20,
                       km_ab = 20,
                       site_name = "Brazil", 
                       progress = T)

# convert subsets to raster
NDVI_r <- mt_to_raster(df = brazil_NDVI)
QA_r <- mt_to_raster(df = brazil_QA)

# creating a mask to remove unreliable pixels from NDVI data and apply mask to data
msk <- QA_r
msk[(QA_r < 0 | QA_r > 1)] <- NA 

NDVI_masked <- mask(NDVI_r, msk, maskvalue=NA, updatevalue=NA)

# converting MODIS dates to date object
dates <- as.Date(names(NDVI_masked), "X%Y.%m.%d")

## first statistical calculations of NDVI raster brick

# 1 Calculating monthly averaged values over the years

# create monthly indices to prepare it for stackApply, which takes the means for all the monthly NDVI over the years
mths <- as.numeric(format(dates, format = "%m")) # extract month from dates object

# stackApply -> calculating monthly means of all NDVI values over the years
MonthNDVI<- stackApply(NDVI_masked, mths, fun = mean)

names(MonthNDVI) <- month.abb # renaming month numbers to abbreviated month names

cols <- colorRampPalette(c("red", "yellow", "darkgreen"))
levelplot(MonthNDVI, col.regions=cols, main = "Average monthly NDVI values") # raster plot of averaged monthly NDVI values
bwplot(MonthNDVI) # violin plot of averaged monthly NDVI values

# 2 Calculating yearly averaged values
yrs <- as.numeric(format(dates, format = "%Y")) # extract years from dates object

# stackApply -> calculating yearly means of all NDVI values 
YearlyNDVI<- stackApply(NDVI_masked, yrs, fun = mean)

levelplot(YearlyNDVI, col.regions=cols, main = "Average yearly NDVI values") # raster plot of averaged yearly NDVI values
bwplot(YearlyNDVI) # violin plot of averaged yearly NDVI values

## applying bfastmonitor to cleaned Raster Brick to find out when changes happened (breakpoints) and how strong (magnitudes of breaks)
# helper function from the bfastmonitor documentation to be used with the calc() function
fn_bfm <- function(x, timestamps = dates) {
  ndvi <- bfastts(x, timestamps, type = c("16-day"))
  ndvi <- window(ndvi)/10000
  bfm <- bfastmonitor(data = ndvi, start = c(2005, 1), history = c("ROC")) # detecting all changes starting from 2005
  return(c(breakpoint = bfm$breakpoint, magnitude = bfm$magnitude))
}

NDVI_bfm <- raster::calc(NDVI_masked, fun=fn_bfm)
raster::plot(NDVI_bfm) # time of break (layer.1) and magnitude of change (layer.2)

# plotting time of break and magnitude of change separately
# extract breakpoint raster
NDVI_breakpoints <- raster(NDVI_bfm, 1)
months <- changeMonth(NDVI_breakpoints)

cols <- viridis::inferno(12)
breaks <- c(1:12)

pdf("changeMonths.pdf", width=15,height=8)
levelplot(months,at=breaks, col.regions=cols, main="Months of Changes in NDVI")
dev.off()

# extract magnitude raster
NDVI_magnitude <- raster(NDVI_bfm, 2)

# plotting only breakpoint pixels
magnitude_breakpoints <- NDVI_magnitude
magnitude_breakpoints[is.na(NDVI_breakpoints)] <- NA

pdf("changeMagnitudeBP.pdf", width=8,height=8)
rasterVis::levelplot(magnitude_breakpoints, main="Magnitude: breakpoints")
dev.off()
