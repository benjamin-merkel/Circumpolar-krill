library(raster)
library(sf)
library(RColorBrewer)
library(viridis)
library(spatialEco)
library(ecospat)

# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 

sf_use_s2(F)

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

model.domain       <- circumpolar

# ccamlr_domains      <- st_read("data/CCAMLR/bm_mpa_planningDomains.shp")
# ccamlr_domains      <- st_transform(ccamlr_domains, south_pole_equal_area.proj)
# ccamlr_domains_crop <- st_intersection(ccamlr_domains, circumpolar)
# rownames(ccamlr_domains) <- ccamlr_domains$OBJECTID <- ccamlr_domains$Index

domains <- readRDS("data/map data/Marine Protected Area planning domains.rds")
domains <- st_transform(domains, south_pole_equal_area.proj)


env <-readRDS("data/Environmental covariate stack.rds")
env[["NSIDC_ice_retreat"]][env[["NSIDC_ice_retreat"]]==46] <- NA # 15 Feb
env[["NSIDC_ice_advance"]][env[["NSIDC_ice_advance"]]==365+45] <- NA # 14 Feb
icer <- env[["NSIDC_ice_retreat"]]
iced <- env[["NSIDC_ice_duration"]]
icer[iced<7]<-NA
icer <- mask(icer, ice_shelf, inverse = T)

plot(icer)

shelf <- env[['bath']]
shelf[shelf < -1000] <- 1
shelf[shelf < 1]     <- 0

## load
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spatial cross validation  ----

species                 <- "E.crystallorophias"
species                 <- "E.superba"

e_mean <- raster(paste0("data/",species,"_circumpolar_GBM_MEDIAN_3000_ensemble_TSS_weighted_Dec-Mar.tif"))
e_sd   <- raster(paste0("data/",species,"_circumpolar_GBM_SD_3000_ensemble_Dec-Mar.tif"))
e_dat  <- readRDS(paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
e_eval <- readRDS(paste0("data/",species,"_circumpolar_GBM_100_ensemble model evaluation_Dec-Mar.rds"))


if(species == "E.crystallorophias"){
  icer <- crop(icer, e_mean)
  e_mean[is.na(icer)]<-NA
  e_sd[is.na(icer)]<-NA
}


presence.test        <- st_as_sf(e_dat)
presence.test$domain <- as.character(st_intersects(presence.test, domains))
presence.test$domain <- as.numeric(substr(presence.test$domain,1,1))
presence.test$predict<- extract(e_mean, as_Spatial(presence.test))
pt <- presence.test  <- presence.test[!is.na(presence.test$presence_absence),]
presence.test        <- presence.test[!is.na(presence.test$predict),]
presence.test        <- presence.test[!is.na(presence.test$domain),]

out <- data.frame(domain          = 1:9, 
                  presence.unique = as.numeric(table(factor(pt$domain[pt$presence_absence==1],levels=c(1:9)))), 
                  absence.unique  = as.numeric(table(factor(pt$domain[pt$presence_absence==0],levels=c(1:9)))), 
                  cbi             = NA, 
                  auc             = NA)
out <- out[out$presence.unique > 0,]

#### CBI
#It is continuous and varies between -1 and +1. 
#Positive values indicate a model which present predictions are 
#consistent with the distribution of presences in the evaluation dataset, 
#values close to zero mean that the model is not different from a random model, 
#negative values indicate counter predictions, i.e., predicting poor quality 
#areas where presences are more frequent (Hirzel et al. 2006).
for(i in out$domain){
  cat("\r",i)
  domain.result          <- mask(e_mean, domains[domains$domain == i,])
  presence.domain        <- presence.test[presence.test$domain==i & presence.test$presence_absence==1 & presence.test$domain==i,]
  out$cbi[out$domain==i] <- ecospat.boyce(values(domain.result)[!is.na(values(domain.result))], presence.domain$predict)$cor
  
  auc_data               <- data.frame(Observed = (presence.test$presence_absence[presence.test$domain==i]))
  auc_data               <- cbind(plotID=1,auc_data, Predicted1=presence.test$predict[presence.test$domain==i])
  out$auc[out$domain==i] <- as.numeric(PresenceAbsence::auc(auc_data,which.model=1)[1])
}

eb       <- ecospat.boyce(values(e_mean)[!is.na(values(e_mean))], presence.test$predict[presence.test$presence_absence==1])$cor
auc_data <- data.frame(Observed = presence.test$presence_absence)
auc_data <- cbind(plotID=1,auc_data, Predicted1=presence.test$predict)
au       <- PresenceAbsence::auc(auc_data,which.model=1)

new.line                 <- out[1,]
new.line$domain          <- "all"
new.line$presence.unique <- nrow(pt[presence.test$presence_absence==1,])
new.line$absence.unique  <- nrow(pt[presence.test$presence_absence==0,])
new.line$cbi             <- eb
new.line$auc             <- au$AUC

out <- rbind(out, new.line)
out
write.csv(out, file = paste0("data/",species,"_CCAMLR domain cross validation.csv"))



