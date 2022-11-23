library(raster)
library(terra)
library(sf)
library(concaveman)
library(ncdf4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### DATA EXTENT -----------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sf_use_s2(F)

# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 
# EPSG:3412
inputcrs_ice <- CRS('+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=km +no_defs')

grat.40s             <- st_read("data/map data/ne_50m_graticules_5.shp")
grat.40s             <- st_transform(grat.40s[grat.40s$direction=="S" & grat.40s$degrees==40,][,1], south_pole_equal_area.proj)
grat.40s             <- rbind(grat.40s, grat.40s[1,])
grat.40s             <- concaveman(grat.40s) # works
grat.40s             <- st_buffer(grat.40s,dist = 0)
grat.40s             <- as_Spatial(st_zm(grat.40s))
background_raster40  <- raster(grat.40s,res=10)
values(background_raster40) <- 1:ncell(background_raster40)
# background_raster40  <- rast(background_raster40)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SEA ICE CLIMATOLOGY 1978 - 2018 ---------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file.directory <- "data/env covariates/NSIDC-0192 Sea Ice Trends and Climatologies from SMMR and SSMI-SSMIS, Version 3/"
ice.files  <- list.files(file.directory)
ice.files  <- ice.files[!grepl("2018",ice.files)]
pers.files <- ice.files[!grepl("mean",ice.files)]
ice.files  <- ice.files[grepl("mean",ice.files)]

ice <- pers <- vector(mode="list")
for(i in 1:length(ice.files)){
  ic <- readBin(paste0(file.directory, ice.files[i]), what = 'integer', n=104912, size = 1, endian = "little")
  ic[ic < 0]<-NA
  ic <- matrix(ic,ncol = 316,nrow = 332, byrow = T)
  #ice <- raster(ice,  crs= inputcrs_ice, xmn = -3950, xmx = 3950, ymn = -4150, ymx = 4150)
  ic <- raster(ic,  crs= inputcrs_ice, xmn = -3950, xmx = 3950, ymn = -3950, ymx = 4350)
  
  ice[[i]] <- ic
  
  ic <- readBin(paste0(file.directory, pers.files[i]), what = 'integer', n=104912, size = 1, endian = "little")
  ic[ic < 0]<-NA
  ic <- matrix(ic,ncol = 316,nrow = 332, byrow = T)
  ic <- raster(ic,  crs= inputcrs_ice, xmn = -3950, xmx = 3950, ymn = -3950, ymx = 4350)
  
  pers[[i]] <- ic
  
}

#NSIDC_ice[NSIDC_ice<15] <-0
NSIDC_ice_conc <- projectRaster((mean(stack(ice))), background_raster40)
NSIDC_ice_conc[is.na(NSIDC_ice_conc)] <-0
names(NSIDC_ice_conc) <- 'NSIDC_ice_conc'

NSIDC_ice_pers <- projectRaster((mean(stack(pers))), background_raster40)
NSIDC_ice_pers[is.na(NSIDC_ice_pers)] <-0
names(NSIDC_ice_pers) <- 'NSIDC_ice_pers'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SEA ICE BREAKUP TIMING 1978 - 2020 ------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# retreat and advance as days since 15 February

NSIDC_ice_retreat <- raster("data/env covariates/sea ice retreat climatology.tif")
NSIDC_ice_retreat <- projectRaster((NSIDC_ice_retreat),background_raster40)
names(NSIDC_ice_retreat) <- 'NSIDC_ice_retreat'

NSIDC_ice_advance <- raster("data/env covariates/sea ice advance climatology.tif")
NSIDC_ice_advance <- projectRaster((NSIDC_ice_advance),background_raster40)
names(NSIDC_ice_advance) <- 'NSIDC_ice_advance'

NSIDC_ice_duration <- raster("data/env covariates/sea ice duration climatology.tif")
NSIDC_ice_duration <- projectRaster((NSIDC_ice_duration),background_raster40)
names(NSIDC_ice_duration) <- 'NSIDC_ice_duration'

NSIDC_ice_edge <- raster("data/env covariates/sea ice edge duration climatology.tif")
NSIDC_ice_edge <- projectRaster((NSIDC_ice_edge),background_raster40)
names(NSIDC_ice_edge) <- 'NSIDC_ice_edge'
NSIDC_ice_edge[is.na(NSIDC_ice_edge)] <- 0

NSIDC_spring_edge <- raster("data/env covariates/sea ice edge duration climatology Nov-Jan.tif")
NSIDC_spring_edge <- crop(projectRaster((NSIDC_spring_edge), background_raster40), background_raster40)
NSIDC_spring_edge[is.na(NSIDC_spring_edge)] <- 0
NSIDC_spring_edge <- mask(NSIDC_spring_edge, NSIDC_ice_duration)
NSIDC_spring_edge[is.na(NSIDC_spring_edge)] <- 0
names(NSIDC_spring_edge) <- 'NSIDC_spring_edge'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### BATHYMETRY ------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rerun = F
if(rerun){
  gebco1 <- raster('data/env covariates/gebco_2020_n0.0_s-90.0_w-180.0_e-90.0.tif')
  gebco2 <- raster('data/env covariates/gebco_2020_n0.0_s-90.0_w-90.0_e0.0.tif')
  gebco3 <- raster('data/env covariates/gebco_2020_n0.0_s-90.0_w0.0_e90.0.tif')
  gebco4 <- raster('data/env covariates/gebco_2020_n0.0_s-90.0_w90.0_e180.0.tif')
  gebco1 <- aggregate(crop(gebco1,extent(-180,180,-90,-30)), factor=100)
  gebco2 <- aggregate(crop(gebco2,extent(-180,180,-90,-30)), factor=100)
  gebco3 <- aggregate(crop(gebco3,extent(-180,180,-90,-30)), factor=100)
  gebco4 <- aggregate(crop(gebco4,extent(-180,180,-90,-30)), factor=100)
  gebco  <- merge(gebco1,gebco2)
  gebco  <- merge(gebco,gebco3)
  gebco  <- merge(gebco,gebco4)
  writeRaster(gebco,"data/gebco 30S and below.tif", format="GTiff", overwrite=TRUE)
}

rerun = F
if(rerun){
  gebco <- raster('data/env covariates/gebco 30S and below.tif')
  gebco <- resample(gebco, raster(crs=4326, resolution=0.01, xmn=-180,xmx=180,ymn=-90,ymx=30))
  gebco <- project(rast(gebco), rast(background_raster40))
  
  ibsco <- raster("data/env covariates/IBCSO_v2_ice-surface.tif")
  ibsco <- projectRaster((ibsco), background_raster40)
  
  bath <- ibsco
  bath[is.na(bath)] <- raster(gebco)[is.na(bath)]
  
  writeRaster((bath),"data/env covariates/IBSCO2022-GEBCO merge.tif", format="GTiff", overwrite=TRUE)
} else {bath <- raster("data/env covariates/IBSCO2022-GEBCO merge.tif")}

shelf                 <- bath
shelf[shelf <= -1000] <- 1
shelf[shelf != 1]     <- NA
bound                 <- terra::boundaries(shelf)
bound[bound==0]       <- NA
d                     <- terra::distance(bound)
dis_1000              <- -d
dis_1000[shelf == 1]  <- d[shelf == 1]
names(dis_1000)       <- 'dis_1000'

bath[bath>0]          <- NA
names(bath)           <- 'bath'

slope                 <- terrain(bath,opt = "slope")
names(slope)          <- 'slope'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### coastline ---------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

land                  <- st_read("data/map data/ne_10m_land.shp")
land                  <- st_transform(st_crop(land, extent(-180,180,-90,0)), south_pole_equal_area.proj)
coast                 <- mask(background_raster40, land)
coast[!is.na(coast)]  <- 1
dis_coastline         <- terra::distance(coast)
names(dis_coastline)  <- "dis_coastline"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### distance to coastal polynyas ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pol_ice.pol <- st_read("data/env covariates/coastal polynyas/coastal polynya climatology.shp")

pol_ice.pol         <- st_transform(pol_ice.pol, south_pole_equal_area.proj)
dis_ice.pol         <- mask(background_raster40, pol_ice.pol)
dis_ice.pol[!is.na(dis_ice.pol)]  <- 1
dis_ice.pol         <- terra::distance(dis_ice.pol)
names(dis_ice.pol)  <- "dis_ice.pol"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### mixed layer depth -------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# from Pellichero et al. 2017
# ship data: 1906-2012; argo data: 2002-2014; seal data: 2004-2014?

rerun = F
if(rerun){
  nc4 <- nc_open('data/env covariates/Pellichero et al. 2017/CLIM_MLD_Pellichero_etal2016.nc')
  temp <- ncvar_get(nc4,"ML_TEMPERATURE")
  sal <- ncvar_get(nc4,"ML_SALINITY")
  pressure <- ncvar_get(nc4,"ML_MAX_PRESSURE")
  lon <- ncvar_get(nc4,"lon")
  lat <- ncvar_get(nc4,"lat")
  month <- ncvar_get(nc4,'time')
  nc_close(nc4)
  
  
  lon2 <- lon
  lon2[lon2>180] <- lon2[lon2>180]-360
  dat <- data.frame(lon      = rep(lon2, length(lat)),
                    lat      = sort(rep(lat, length(lon))),
                    pressure = c(pressure),
                    temp     = c(temp),
                    sal      = c(sal))
  coordinates(dat) <- dat[,1:2]
  proj4string(dat) <- 4326
  
  empty.mld <- raster(crs =proj.lonlat, vals=NULL, res=0.5,
                      xmn=-180, xmx=180, ymn=-80, ymx=-40)
  
  dat_ras <- rasterize(dat, empty.mld, fun=mean, field="pressure")
  writeRaster(dat_ras,"data/env covariates/Mixed Layer Depth Jan-Mar climatology.tif", format="GTiff", overwrite=TRUE)
  
  dat_ras <- rasterize(dat, empty.mld, fun=mean, field="temp")
  writeRaster(dat_ras,"data/env covariates/Mixed Layer Temperature Jan-Mar climatology.tif", format="GTiff", overwrite=TRUE)
  
  dat_ras <- rasterize(dat, empty.mld, fun=mean, field="sal")
  writeRaster(dat_ras,"data/env covariates/Mixed Layer Salinity Jan-Mar climatology.tif", format="GTiff", overwrite=TRUE)
}

Pellichero_ml_depth <- rast('data/env covariates/Mixed Layer Depth Jan-Mar climatology.tif')
Pellichero_ml_depth <- projectRaster(raster(Pellichero_ml_depth), raster(background_raster40))
# Pellichero_ml_depth <- project(Pellichero_ml_depth,background_raster40, method = "bilinear")
names(Pellichero_ml_depth) <- 'Pellichero_ml_depth'

Pellichero_ml_temp <- rast('data/env covariates/Mixed Layer Temperature Jan-Mar climatology.tif')
Pellichero_ml_temp <- projectRaster(raster(Pellichero_ml_temp),raster(background_raster40))
names(Pellichero_ml_temp) <- 'Pellichero_ml_temp'

Pellichero_ml_sal <- rast('data/env covariates/Mixed Layer Salinity Jan-Mar climatology.tif')
Pellichero_ml_sal <- projectRaster(raster(Pellichero_ml_sal),raster(background_raster40))
names(Pellichero_ml_sal) <- 'Pellichero_ml_sal'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Chla, phenology and bloom stats --------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1998 - 2019
OCCCI_bd <- raster("data/env covariates/Chla bd mean climatology 98_18.tif")
OCCCI_bd <- projectRaster((OCCCI_bd),background_raster40)
names(OCCCI_bd) <- 'OCCCI_bd'

OCCCI_bt <- raster("data/env covariates/Chla bt mean climatology 98_18.tif")
OCCCI_bt <- projectRaster((OCCCI_bt),background_raster40)
names(OCCCI_bt) <- 'OCCCI_bt'

OCCCI_bi <- raster("data/env covariates/Chla bi mean climatology 98_18.tif")
OCCCI_bi <- projectRaster((OCCCI_bi),background_raster40)
names(OCCCI_bi) <- 'OCCCI_bi'

OCCCI_be <- raster("data/env covariates/Chla be mean climatology 98_18.tif")
OCCCI_be <- projectRaster((OCCCI_be),background_raster40)
names(OCCCI_be) <- 'OCCCI_be'

OCCCI_bm <- raster("data/env covariates/Chla bm mean climatology 98_18.tif")
OCCCI_bm <- projectRaster((OCCCI_bm),background_raster40)
names(OCCCI_bm) <- 'OCCCI_bm'

OCCCI_bs <- raster("data/env covariates/Chla bs mean climatology 98_18.tif")
OCCCI_bs <- projectRaster((OCCCI_bs),background_raster40)
names(OCCCI_bs) <- 'OCCCI_bs'

OCCCI_bv <- raster("data/env covariates/Chla bv mean climatology 98_18.tif")
OCCCI_bv <- projectRaster((OCCCI_bv),background_raster40)
names(OCCCI_bv) <- 'OCCCI_bv'

OCCCI_ba <- raster("data/env covariates/Chla ba mean climatology 98_18.tif")
OCCCI_ba <- projectRaster((OCCCI_ba),background_raster40)
names(OCCCI_ba) <- 'OCCCI_ba'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### WOA2018 ------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WOA_temp_0 <- raster("data/env covariates/WOA 18 temperature 0m Dec-Mar climatology.tif")
WOA_temp_0 <- projectRaster(WOA_temp_0,background_raster40)
names(WOA_temp_0) <- 'WOA_temp_0'

WOA_temp_200 <- raster("data/env covariates/WOA 18 temperature 200m Dec-Mar climatology.tif")
WOA_temp_200 <- projectRaster(WOA_temp_200,background_raster40)
names(WOA_temp_200) <- 'WOA_temp_200'

WOA_sal_0 <- raster("data/env covariates/WOA 18 salinity 0m Dec-Mar climatology.tif")
WOA_sal_0 <- projectRaster(WOA_sal_0,background_raster40)
names(WOA_sal_0) <- 'WOA_sal_0'

WOA_sal_200 <- raster("data/env covariates/WOA 18 salinity 200m Dec-Mar climatology.tif")
WOA_sal_200 <- projectRaster(WOA_sal_200,background_raster40)
names(WOA_sal_200) <- 'WOA_sal_200'

WOA_ox_0 <- raster("data/env covariates/WOA 18 dissloved oxygen 0m Dec-Mar climatology.tif")
WOA_ox_0 <- projectRaster(WOA_ox_0,background_raster40)
names(WOA_ox_0) <- 'WOA_ox_0'

WOA_ox_200 <- raster("data/env covariates/WOA 18 dissloved oxygen 200m Dec-Mar climatology.tif")
WOA_ox_200 <- projectRaster(WOA_ox_200,background_raster40)
names(WOA_ox_200) <- 'WOA_ox_200'

WOA_si_0 <- raster("data/env covariates/WOA 18 silicate 0m Dec-Mar climatology.tif")
WOA_si_0 <- projectRaster(WOA_si_0,background_raster40)
names(WOA_si_0) <- 'WOA_si_0'

WOA_si_200 <- raster("data/env covariates/WOA 18 silicate 200m Dec-Mar climatology.tif")
WOA_si_200 <- projectRaster(WOA_si_200,background_raster40)
names(WOA_si_200) <- 'WOA_si_200'

WOA_ni_0 <- raster("data/env covariates/WOA 18 nitrate 0m Dec-Mar climatology.tif")
WOA_ni_0 <- projectRaster(WOA_ni_0,background_raster40)
names(WOA_ni_0) <- 'WOA_ni_0'

WOA_ni_200 <- raster("data/env covariates/WOA 18 nitrate 200m Dec-Mar climatology.tif")
WOA_ni_200 <- projectRaster(WOA_ni_200,background_raster40)
names(WOA_ni_200) <- 'WOA_ni_200'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## stack all predictors  --------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

env      <- stack((WOA_temp_0), (WOA_temp_200),
                  (WOA_sal_0), (WOA_sal_200),
                  (Pellichero_ml_depth), (Pellichero_ml_temp), (Pellichero_ml_sal),
                  (OCCCI_bd), (OCCCI_bi), (OCCCI_be), (OCCCI_bt), (OCCCI_bm), (OCCCI_bs), (OCCCI_bv), (OCCCI_ba),
                  (NSIDC_ice_retreat), (NSIDC_ice_advance), (NSIDC_ice_duration), 
                  (NSIDC_ice_edge), (NSIDC_spring_edge), (NSIDC_ice_conc), (NSIDC_ice_pers),
                  bath, dis_1000, dis_coastline, dis_ice.pol,
                  (WOA_ox_0), (WOA_ox_200),
                  (WOA_si_0), (WOA_si_200),
                  (WOA_ni_0), (WOA_ni_200))

env[["OCCCI_bm"]]       <- log(env[["OCCCI_bm"]])
env[["OCCCI_bv"]]       <- log(env[["OCCCI_bv"]])
env[["OCCCI_bs"]]       <- log(env[["OCCCI_bs"]])
env[["OCCCI_ba"]]       <- log(env[["OCCCI_ba"]])

e1 <- env[["OCCCI_bi"]]
values(e1)[values(e1)>182 & !is.na(values(e1))] <- values(e1)[values(e1)>182 & !is.na(values(e1))] - 365
env[["OCCCI_bi"]] <- e1

e1 <- env[["OCCCI_bt"]]
values(e1)[values(e1)>182 & !is.na(values(e1))] <- values(e1)[values(e1)>182 & !is.na(values(e1))] - 365
env[["OCCCI_bt"]] <- e1

e1 <- env[["OCCCI_be"]]
values(e1)[values(e1)>182 & !is.na(values(e1))] <- values(e1)[values(e1)>182 & !is.na(values(e1))] - 365
env[["OCCCI_be"]] <- e1

env[["NSIDC_ice_duration"]][is.na(env[["NSIDC_ice_duration"]])] <- 0
env[["NSIDC_ice_retreat" ]][is.na(env[["NSIDC_ice_retreat" ]])] <- 46     # 15 Feb
env[["NSIDC_ice_advance" ]][is.na(env[["NSIDC_ice_advance" ]])] <- 365+45 # 14 Feb

iceNA               <- rast("data/env covariates/sea ice NA.tif")
iceNA               <- terra::project(iceNA, background_raster40)
# iceNA[iceNA==1]     <- 0
# iceNA[is.na(iceNA)] <- 1
# iceNA[iceNA==0]     <- NA
iceNA               <- resample((iceNA), rast(env[["NSIDC_ice_duration"]]))
iceNA               <- raster(iceNA)

env[["NSIDC_ice_duration"]][iceNA==1] <- NA
env[["NSIDC_ice_edge"]]    [iceNA==1] <- NA
env[["NSIDC_spring_edge"]] [iceNA==1] <- NA
env[["NSIDC_ice_retreat" ]][iceNA==1] <- NA
env[["NSIDC_ice_advance" ]][iceNA==1] <- NA
env[["NSIDC_ice_conc"]]    [iceNA==1] <- NA
env[["NSIDC_ice_pers"]]    [iceNA==1] <- NA


saveRDS(env,"data/Environmental covariate stack.rds")


