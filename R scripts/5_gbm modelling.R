#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages

library(sf)
library(sp)
library(raster)
library(terra)
library(concaveman)
library(biomod2)

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

model.domain       <- circumpolar


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load env covariates ------

env          <- readRDS("data/Environmental covariate stack.rds")
ice_duration <- env[["NSIDC_ice_duration"]]
env.selected <- readRDS("data/Environmental covariates selected.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define species ------

species     <- "E.crystallorophias"
species     <- "E.superba"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load response data -------

setwd(base_dir)
data <- as_Spatial(readRDS(paste0("data/",species," response data.rds")))
data <- data[!is.na(data$year),]
data <- data[!is.na(data$month),]
data <- data[!is.na(data$presence_absence),]
data <- data[!is.na(data$doy),]

# remove data outside the model domain and time period considered
data <- data[data$month %in% c(12,1:3),] # 1 December - 30 March
# data <- data[data$month %in% 1:3,] 
data <- data[data$year > 1970,]
data <- as_Spatial(st_intersection(st_as_sf(data), model.domain))

# extract environmental covariates for each observation and remove locations with NA values
data <- cbind(data, extract(env.selected, data))
datax <- data[!complete.cases(data.frame(data[,names(env.selected)])),]
if(species == "E.superba") env.rows <- 8:length(names(data))
if(species == "E.crystallorophias") env.rows <- 7:length(names(data))

# extract environmental covariates for each observation with NA from neighbour cell
for(i in env.rows){
  if(length(data@data[,i][is.na(data@data[,i])])>0)   data@data[,i][is.na(data@data[,i])] <- unlist(lapply(extract(env.selected[[names(data[,i])]], data[is.na(data@data[,i]),], method="bilinear", buffer = 10),mean,na.rm=T))
}

# remove locations with NA values
data <- data[complete.cases(data.frame(data[,names(env.selected)])),]



# thin out presence and absence points by removing spatial duplicates at the resolution of the covariates (i.e. background_raster)
presence               <- data[data$presence_absence %in% 1,]
presence$background.bin<- extract(background_raster, presence)
presence2              <- presence[!duplicated(presence$background.bin),]
absence                <- data[data$presence_absence %in% 0,]
absence$background.bin <- extract(background_raster, absence)
absence                <- absence[!absence$background.bin %in% unique(presence$background.bin),]
absence2               <- absence[!duplicated(absence$background.bin),]

data <- rbind(presence2, absence2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODELLING -------

# define model settings
myBiomodOption <- BIOMOD_ModelingOptions(GBM=list(n.trees = 3000, 
                                                  cv.folds= 10,
                                                  n.cores = 3,
                                                  interaction.depth = 5,
                                                  shrinkage = 0.01))
biomod_gbm <- gbm_Eval <- dat.used <- NULL

nruns=10
dat.used        <- data[,"presence_absence"]
names(dat.used) <- substr(species, 1, 3)


biomod_data <- BIOMOD_FormatingData(resp.var = dat.used,
                                    resp.name= substr(species, 1, 3),
                                    expl.var = data[,names(env.selected)])

setwd(gbm_output_dir)
biomod_gbm <- BIOMOD_Modeling(bm.format         = biomod_data,
                              models            = "GBM",
                              modeling.id       = 'GBM Dec-Mar',
                              nb.rep            = nruns,
                              data.split.perc   = 70,
                              prevalence        = NULL,
                              weights           = NULL,
                              var.import        = 10,
                              metric.eval       = c('TSS','ROC'),
                              save.output       = T,
                              scale.models      = F,
                              do.full.models    = F,
                              bm.options        = myBiomodOption,
                              data.split.table  = NULL,
                              nb.cpu            = 4,
                              seed.val          = NULL,
                              do.progress       = TRUE) 


gbm_Eval <- get_evaluations(biomod_gbm)

biomod_ensemble <- BIOMOD_EnsembleModeling(bm.mod = biomod_gbm,
                                           models.chosen = 'all',
                                           em.by = "all",
                                           metric.eval = 'ROC',
                                           prob.median = TRUE,
                                           prob.cv = T,
                                           prob.ci = F,
                                           nb.cpu = 4,
                                           prob.mean.weight	= T)

ensemble_Eval <- get_evaluations(biomod_ensemble)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load models

biomod_gbm <- NULL
setwd(gbm_output_dir)
temp.path <- list.files(pattern="Dec-Mar.models.out", recursive = TRUE)
temp.path <- temp.path[grepl(paste0(substr(species,1,3),"."), temp.path)]
if(exists("myBiomodOuttmp")) { rm(myBiomodOuttmp) }
biomod_gbm <- get(load(temp.path))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model evaluation metric ---------

gbm_Eval<- get_evaluations(biomod_gbm)
summary(gbm_Eval[2,1,,,])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check variable importance across all models -----

var_imp <- get_variables_importance(biomod_gbm)
vi3 <- apply(var_imp,1,c)
vi4 <- apply(vi3, 2, mean)

setwd(base_dir)
png(paste0("figures/",species ," GBM circumpolar variable importance_Dec-Mar.png"), res = 800, width=17, height = 17, units="cm")
opar <- par(mfrow=c(1,1),mar=c(4,10,1,1))
boxplot(vi3[,order(vi4)], horizontal=T, las=1,
        col = grey(1), border = "skyblue3", lty = 1, lwd = 2,
        ylab = "", xlab = "variable importance")
points(vi4[order(vi4)],1:length(vi4),cex=1.5,pch=19,col="darkred")
par(opar)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict model output over model domain for each model and ensemble ----

setwd(gbm_output_dir)
proj1 <- BIOMOD_Projection(bm.mod = biomod_gbm,
                           new.env = env.selected,
                           proj.name ='current',
                           selected.models = "all",
                           binary.meth =NULL,
                           compress =F,
                           clamping.mask = F,
                           output.format ='.grd')

biomod_gbm_Proj <- proj1@proj.out@val

proj_ensemble <- BIOMOD_EnsembleForecasting(bm.em  = biomod_ensemble,
                                            bm.proj = proj1,
                                            selected.models = 'all',
                                            compress = T)
mod_proj_ensemble <- get_predictions(proj_ensemble)
ensemble_median   <- mod_proj_ensemble[[7]] # median model ensemble
ensemble_mean     <- mod_proj_ensemble[[8]] # mean model ensemble
ensemble_cv       <- mod_proj_ensemble[[6]] # mean model ensemble

plot(ensemble_mean, zlim=c(0,1000))
plot(ensemble_mean_cv)

# # standard deviation of model prediction based on 100 models using 70% of available data
gbm_ras_sd   <- calc(raster::subset(biomod_gbm_Proj, grep('RUN', names(biomod_gbm_Proj), value = T)), sd)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE ----

setwd(base_dir)
writeRaster(ensemble_mean,   filename=paste0("data/",species,"_circumpolar_GBM_MEAN_ensemble_ROC_weighted_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)
writeRaster(ensemble_median, filename=paste0("data/",species,"_circumpolar_GBM_MEDIAN_ensemble_ROC_weighted_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)
writeRaster(ensemble_cv,     filename=paste0("data/",species,"_circumpolar_GBM_CV_ensemble_ROC_weighted_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)
writeRaster(gbm_ras_sd,      filename=paste0("data/",species,"_circumpolar_GBM_SD_ensemble_unweighted_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)
writeRaster(biomod_gbm_Proj, filename=paste0("data/",species,"_circumpolar_GBM_10model_prediction_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)

saveRDS(biomod_gbm_Proj,     paste0("data/",species,"_circumpolar_GBM_raster_predictions_Dec-Mar.rds"))

saveRDS(data,                paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
saveRDS(vi3,                 paste0("data/",species,"_circumpolar_GBM_estimated variable importance_Dec-Mar.rds"))
saveRDS(gbm_Eval,            paste0("data/",species,"_circumpolar_GBM_all model evaluation_Dec-Mar.rds"))
saveRDS(ensemble_Eval,       paste0("data/",species,"_circumpolar_GBM_ensemble model evaluation_Dec-Mar.rds"))
