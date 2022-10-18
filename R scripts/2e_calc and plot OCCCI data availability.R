#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages

library(sf)
library(sp)
library(raster)
library(terra)
library(concaveman)
library(biomod2)
library(stringr)
library(RColorBrewer)
library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define directories

base_dir       <- getwd()
gbm_output_dir <- paste0(base_dir, "/data/gbm")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load mapping data and create polygon for pseudo abscences  -------

sf_use_s2(F)

# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 

ice_shelf          <- st_read("data/map data/Ice_shelf.shp")
ice_shelf          <- st_transform(ice_shelf, south_pole_equal_area.proj)
ice_shelf          <- st_buffer(ice_shelf,dist = 0)

grat.50S           <- st_read("data/map data/ne_50m_graticules_5.shp")
grat.50S           <- st_transform(grat.50S[grat.50S$direction=="S" & grat.50S$degrees==50,][,1], south_pole_equal_area.proj)
grat.50S           <- rbind(grat.50S, grat.50S[1,])
grat.50S           <- concaveman(grat.50S) # works
grat.50S           <- st_buffer(grat.50S,dist = 0)
circumpolar        <- st_zm(grat.50S, drop = T,what = "ZM")

land               <- st_read("data/map data/ne_10m_land.shp")
land               <- st_transform(st_crop(land, extent(-180,180,-90,-40)), south_pole_equal_area.proj)
land               <- st_intersection(land, circumpolar)

background_polygon <- st_difference((grat.50S),st_geometry(ice_shelf))
background_polygon <- st_difference(background_polygon,st_geometry(land))
background_raster  <- raster(as_Spatial(background_polygon),res=10)
values(background_raster) <- 1:ncell(background_raster)


species                 <- "E.crystallorophias"
ec_mean <- raster(paste0("data/",species,"_circumpolar_GBM_MEDIAN_3000_ensemble_TSS_weighted_Dec-Mar.tif"))
ec_sd   <- raster(paste0("data/",species,"_circumpolar_GBM_SD_3000_ensemble_Dec-Mar.tif"))
ec_dat  <- readRDS(paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
ec_eval <- readRDS(paste0("data/",species,"_circumpolar_GBM_100_ensemble model evaluation_Dec-Mar.rds"))


icer <- crop(icer, ec_mean)
ec_mean[is.na(icer)]<-NA
ec_sd[is.na(icer)]<-NA

species                 <- "E.superba" 
es_mean <- raster(paste0("data/",species,"_circumpolar_GBM_MEDIAN_3000_ensemble_TSS_weighted_Dec-Mar.tif"))
es_sd   <- raster(paste0("data/",species,"_circumpolar_GBM_SD_3000_ensemble_Dec-Mar.tif"))
es_dat  <- readRDS(paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
es_eval <- readRDS(paste0("data/",species,"_circumpolar_GBM_100_ensemble model evaluation_Dec-Mar.rds"))


es_tss_cutoff <- mean(es_eval[2,2,1,,1]) 
ec_tss_cutoff <- mean(ec_eval[2,2,1,,1]) 

es_bin <- es_mean
es_bin[es_bin< es_tss_cutoff] <- 0
es_bin[es_bin>0] <- 1

ec_bin <- ec_mean
ec_bin[ec_bin< ec_tss_cutoff] <- 0
ec_bin[ec_bin>0] <- 2


file_directory <- "D:/OC-CCI/8-day_Chla"
oc_cci <- list.files(file_directory, pattern = ".nc")
date  <- as.Date(strptime(str_split_fixed(oc_cci,"-",8)[,7], "%Y%m%d"))

# year_array <- array(NA , dim = c(874,874,length(oc_cci)))
year_array <- array(NA , dim = c(438,438,length(oc_cci)))
pp=4
for(pp in 1:4){

  if(pp==1) ext_p <- extent(-10,6400,  -10,6400)
  if(pp==2) ext_p <- extent(-10,6400,  -6400,10)
  if(pp==3) ext_p <- extent(-6400,10, -6400,10)
  if(pp==4) ext_p <- extent(-6400,10, -10,6400)

  
  for(yy in 1:length(date)){
    cat(yy,' - ',length(date),'\r')
    dat1      <- raster(paste(file_directory ,oc_cci[yy], sep="/"), varname = "chlor_a")
    dat1      <- crop(dat1, extent(-180,180,-90,-50))
    dat1_proj <- projectRaster(dat1, background_raster)
    # dat1_proj <- projectRaster(dat1, background_raster)
    dat1_proj <- projectRaster(dat1, crop(background_raster, ext_p))
    dat1_proj_int <- focal(dat1_proj, w=matrix(rep(1,9),ncol=3), fun=mean,na.rm=T,NAonly=T)
    # dat1_proj_int < crop(dat1_proj_int, extent(0,6400,0,6400))
    year_array[,,yy] <- raster::as.matrix(dat1_proj_int)
  }
  
  start <- Sys.time()
  na_mat <- apply(year_array, c(1,2), FUN = function(x){
    if(any(!is.na(x))){
      y <- length(x[!is.na(x)])
    } else {
      y <- 0
    }
    return(y)
  })
  cat("\r",Sys.time()-start, "    ")
  # saveRDS(na_mat,"data/env covariates/OCCCI number of cells with NA.rds")
  saveRDS(na_mat,paste0("data/env covariates/OCCCI number of cells with NA part ",pp,".rds"))
}


na_mat <- vector(mode = "list",length = 4)
for(pp in 1:4){
  na_mat[[pp]] <- readRDS(paste0("data/env covariates/OCCCI number of cells with NA part ",pp,".rds"))
}


bm_ras <- raster("data/env covariates/Chla bm mean climatology 98_18.tif")
bm_ras <- projectRaster(bm_ras, background_raster)
na_ras <- bm_ras
x1 <- y1 <- 2:438
x2 <- y2 <- 1:437
na2 <-  (cbind(rbind(na_mat[[4]][x2, y2],   # top left
                     na_mat[[3]][x1, y2]),  # bottom left
               rbind(na_mat[[1]][x2, y1],   # top right
                     na_mat[[2]][x1, y1]))) # bottom right

values(na_ras) <- na2
na_ras[is.na(bm_ras)] <- 0
na_ras <- mask(na_ras, as_Spatial(circumpolar))
na_ras <- na_ras / length(date) * 100

es_na <- ec_na <- na_ras
es_na[es_bin==0]     <- NA
es_na[is.na(es_bin)] <- NA
ec_na[ec_bin==0]     <- NA
ec_na[is.na(ec_bin)] <- NA

species_na <- data.frame(na = c(values(es_na), values(ec_na)),
                         spe= c(rep("Antarctic krill",length(es_na)), rep("Ice krill",length(ec_na))))
species_na <- species_na[!is.na(species_na$na),]

ggplot(species_na, aes(spe, na)) + geom_boxplot(aes(colour = spe))
tapply(species_na$na, species_na$spe, summary)
summary(values(na_ras))

png("figures/OC_CCI available data.png", units="cm",res=500,width=18, height=18)
par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(st_geometry(circumpolar))
plot(na_ras, col=brewer.pal(11,"Spectral"),
     main="",axes=F,box=F, zlim = c(0,100), legend = F,add=T)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(na_ras, legend.only=T,
     col=brewer.pal(11,"Spectral"),zlim =c(0, 100),add=T,
     axis.args=list(cex.axis=0.8),
     legend.args=list(text='% of days with data', side=4, font=2, line=2.5, cex=0.8),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
dev.off()