library(raster)
library(sf)


# EPSG:3412
inputcrs_ice <- CRS('+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=km +no_defs')

file.directory <- "data/env covariates/Nihashi & Ohshima 2015/"
ice.prod <- readBin(paste0(file.directory, "2003-2010_ave.icepro.data"), 
                    what = 'double', n = 1264 * 1328, endian = "little", size = 4)
ice.prod <- matrix(ice.prod,ncol = 1264,nrow = 1328, byrow = T) # imax=1264,jmax=1328
ice.prod <- raster(ice.prod,  crs= inputcrs_ice, xmn = -3950, xmx = 3950, ymn = -3950, ymx = 4350)

ice.prod.bin <- ice.prod
ice.prod.bin[ice.prod.bin<4] <- NA
ice.prod.bin[!is.na(ice.prod.bin)] <- 1

ice.prod.pol <- rasterToPolygons(ice.prod.bin, fun = function(x) {x==1}, dissolve = T)
ice.prod.pol <- st_as_sf(ice.prod.pol)

st_crs(ice.prod.pol) <- inputcrs_ice
st_write(obj = ice.prod.pol, dsn = paste0(file.directory, "coastal polynya climatology"), 
         layer = "coastal polynya climatology", driver = 'ESRI Shapefile')
