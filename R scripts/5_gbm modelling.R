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
# species     <- "E.superba"


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

# biomod_data <- BIOMOD_FormatingData(resp.var = data.frame(dat.used)[,1],
#                                     resp.xy  = data.frame(dat.used)[,2:3],
#                                     resp.name= substr(species, 1, 3),
#                                     expl.var = stack(env.selected))


biomod_data <- BIOMOD_FormatingData(resp.var = dat.used,
                                    resp.name= substr(species, 1, 3),
                                    expl.var = data[,names(env.selected)])

setwd(gbm_output_dir)
biomod_gbm <- BIOMOD_Modeling(data              = biomod_data,
                              models            = "GBM",
                              models.options    = myBiomodOption,
                              NbRunEval         = nruns,
                              DataSplit         = 70,
                              Prevalence        = NULL,
                              Yweights          = NULL,
                              VarImport         = 10,
                              models.eval.meth  = c('TSS','ROC'),
                              SaveObj           = T,
                              rescal.all.models = F,
                              do.full.models    = F,
                              modeling.id       = 'GBM Dec-Mar') 

biomod_ensemble <- BIOMOD_EnsembleModeling(modeling.output = biomod_gbm,
                                           chosen.models = 'all',
                                           em.by = 'all',
                                           eval.metric = c('TSS'),
                                           eval.metric.quality.threshold = NULL,
                                           models.eval.meth = c('TSS','ROC'),
                                           prob.median = TRUE )

gbm_Eval <- get_evaluations(biomod_gbm)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load models
# 
# biomod_gbm <- NULL
# setwd(gbm_output_dir)
# temp.path <- list.files(pattern=".models.out", recursive = TRUE)
# temp.path <- temp.path[grepl(paste0(substr(species,1,3),"."), temp.path)]
# temp.path <- temp.path[grepl(paste("E.s.1.GBM Dec-Mar circumpolar.models.out"), temp.path)][1]
# if(exists("myBiomodOuttmp")) { rm(myBiomodOuttmp) }
# biomod_gbm <- get(load(temp.path))
# 
# gbm_Eval<- get_evaluations(biomod_gbm)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check marginal response plots for each covariate across models ------

setwd(gbm_output_dir)
Biomodresponse <- response.plot2(
  models = BIOMOD_LoadModels(biomod_gbm),
  Data = get_formal_data(biomod_gbm, 'expl.var'),
  show.variables = get_formal_data(biomod_gbm,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'mean',
  plot=F,
  legend = F,
  data_species = get_formal_data(biomod_gbm, 'resp.var')
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check variable importance across all models -----

var_imp <- get_variables_importance(biomod_gbm)
vi3 <- t(var_imp[,,1:nruns,"AllData"])
vi4 <- apply(vi3,2,mean)

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
# predict model output over model domain -----

setwd(gbm_output_dir)
proj1 <- BIOMOD_Projection(modeling.output = biomod_gbm,
                           new.env = env.selected,
                           proj.name ='current',
                           selected.models = "all",
                           binary.meth =NULL,
                           compress =F,
                           clamping.mask = F,
                           output.format ='.grd')
biomod_gbm_Proj <- proj1@proj@val

proj_ensemble <- BIOMOD_EnsembleForecasting(biomod_ensemble,
                                            projection.output = proj1,
                                            selected.models = 'all',
                                            compress = 'gzip'
)
mod_proj_ensemble <- get_predictions(proj_ensemble)
EnsembleResult <- mod_proj_ensemble[[2]] #This is the median model ensemble
plot(EnsembleResult, zlim=c(0,1000))

gbm_Ens <- get_evaluations(biomod_ensemble)


# mean model prediction based on 10 models using all available data
gbm_ras_mean <- mean(raster::subset(biomod_gbm_Proj, grep('RUN', names(biomod_gbm_Proj), value = T)))
# standard deviation of model prediction based on 100 models using 70% of available data
gbm_ras_sd   <- calc(raster::subset(biomod_gbm_Proj, grep('RUN', names(biomod_gbm_Proj), value = T)), sd)

# plot model predictions
opar <- par(mfrow=c(2,1))
layer1 <- EnsembleResult
# layer1[env[["NSIDC_ice_duration"]]<10]<-NA
plot(mask(layer1,model.domain), main="Ensemble median")
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)

layer1 <- gbm_ras_mean
# layer1[env[["NSIDC_ice_duration"]]<10]<-NA
plot(mask(layer1,model.domain), main="unweighted mean")
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
par(opar)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE ----

setwd(base_dir)
writeRaster(EnsembleResult, filename=paste0("data/",species,"_circumpolar_GBM_MEDIAN_3000_ensemble_TSS_weighted_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)
writeRaster(gbm_ras_mean  , filename=paste0("data/",species,"_circumpolar_GBM_MEAN_3000_ensemble_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)
writeRaster(gbm_ras_sd    , filename=paste0("data/",species,"_circumpolar_GBM_SD_3000_ensemble_Dec-Mar.tif"), format="GTiff", overwrite=TRUE)

saveRDS(data             , paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
# saveRDS(Biomodresponse   , paste0("data/",species,"_circumpolar_GBM_100_response_curve_data_Dec-Mar.rds"))
saveRDS(vi3              , paste0("data/",species,"_circumpolar_GBM_100_estimated variable importance_Dec-Mar.rds"))
saveRDS(gbm_Eval         , paste0("data/",species,"_circumpolar_GBM_100_ensemble model evaluation_Dec-Mar.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting ----

layer1 <- EnsembleResult 
if(species == "E.crystallorophias") layer1[ice_duration < 10] <- NA
legend.axis.labels = c("0","1")
legend.labels = "probability"
setwd(base_dir)
png(paste0("figures/",species ," GBM circumpolar model prediction mean map_Dec-Mar.png"), 
    res = 800,   width=20, height = 20, units="cm")
par(mar=rep(0,4))
plot(st_geometry(model.domain),lty=2)  
plot(mask(layer1, model.domain),  zlim=c(0,1000), add=T, legend =F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(layer1,legend.only = T, zlim=c(0,1000),  
     axis.args=list(at=c(0, 1000), #brks,
                    labels=legend.axis.labels,#brks/10 ,
                    cex.axis=1,tick=F,line=-0.2),
     legend.args=list(text=legend.labels, side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
dev.off()



# plot marginal response curves for each covariate
ncols <- ceiling((nlayers(env.selected)+1)/4)
setwd(base_dir)
png(paste0("figures/",species ," GBM circumpolar marginal response curves_Dec-Mar.png"),
    units="cm",res=500,width=4*ncols, height=15)
par(mfrow = c(4,ncols),mar=c(2,2,2,0.1))
for(i in unique(Biomodresponse$expl.name)){
  xx <- Biomodresponse[Biomodresponse$expl.name==i,]
  xx <- xx[complete.cases(xx),]
  xx <- xx[!grepl("Full", xx$pred.name),]
  # xx$expl.val <- round(xx$expl.val,0)
  response <- data.frame(mean   = tapply(xx$pred.val, xx$expl.val, mean),
                         median = tapply(xx$pred.val, xx$expl.val, median),
                         sd     = tapply(xx$pred.val, xx$expl.val, sd),
                         se     = tapply(xx$pred.val, xx$expl.val, sd)/sqrt(length(unique(xx$pred.name))),
                         q0.025 = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.025)),
                         q0.25  = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.25)),
                         q0.75  = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.75)),
                         q0.975 = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.975)))
  response$x <- as.numeric(as.character(rownames(response)))
  response$lower.ci <- response$mean - 2 * response$se
  response$upper.ci <- response$mean + 2 * response$se
  
  plot(response$x, response$mean,type="l",lwd=2,ylim=c(0.4,1),main=i,las=1,ylab="",xlab="",cex.axis=0.8,col="white")
  for(j in unique(Biomodresponse$pred.name)){
    xx2 <- xx[xx$pred.name==j,]
    lines(xx2$expl.val, xx2$pred.val, lwd=1, col=rgb(1,0,0,0.2))
  }
  # lines(response$x, frollmean(response$mean, 5), lwd=3, col=1)
  lines(response$x, response$mean, lwd=3, col=1)
  
  rug(biomod_data@data.env.var[,i])
}
plot(1,1,col="white",ann=F,axes=F)
legend("bottomright", legend = c("ensemble mean","single model"),
       lwd=c(4,1),lty=1,col=c(1,2),cex=0.9)

par(opar)
dev.off()

