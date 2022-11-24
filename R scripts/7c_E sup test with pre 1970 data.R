library(sp)
library(sf)
library(raster)
library(terra)
library(ecospat)

domains <- readRDS("data/map data/Marine Protected Area planning domains.rds")
domains <- st_transform(domains, south_pole_equal_area.proj)


rerun = F
if(rerun){
  test.data <- readRDS("data/E.superba response data.rds")
  test.data <- test.data[!is.na(test.data$year),]
  test.data <- test.data[!is.na(test.data$month),]
  test.data <- test.data[!is.na(test.data$presence_absence),]
  test.data <- test.data[!is.na(test.data$doy),]
  
  # remove test.data outside the model domain and time period considered
  test.data <- test.data[test.data$month %in% c(12,1:3),] # 1 December - 30 March
  test.data <- test.data[test.data$year < 1970,]
  test.data <- as_Spatial(st_intersection(st_as_sf(test.data), model.domain))
  
  # thin out presence and absence points by removing spatial duplicates at the resolution of the covariates (i.e. background_raster)
  presence               <- test.data[test.data$presence_absence %in% 1,]
  presence$background.bin<- extract(background_raster, presence)
  presence2              <- presence[!duplicated(presence$background.bin),]
  absence                <- test.data[test.data$presence_absence %in% 0,]
  absence$background.bin <- extract(background_raster, absence)
  absence                <- absence[!absence$background.bin %in% unique(presence$background.bin),]
  absence2               <- absence[!duplicated(absence$background.bin),]
  
  test.data <- rbind(presence2, absence2)
  
  saveRDS(test.data, "data/E.superba_circumpolar_pre1970_data_Dec-Mar.rds")
} else {test.data <- readRDS("data/E.superba_circumpolar_pre1970_data_Dec-Mar.rds")}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spatial cross validation  ----

es_mean <- raster(paste0("data/E.superba_circumpolar_GBM_MEAN_ensemble_ROC_weighted_Dec-Mar.tif"))


presence.test        <- st_as_sf(test.data)
presence.test$domain <- as.character(st_intersects(presence.test,domains))
presence.test$domain <- as.numeric(substr(presence.test$domain,1,1))
presence.test$predict<- extract(es_mean, as_Spatial(presence.test))
presence.test        <- presence.test[!is.na(presence.test$predict),]
presence.test        <- presence.test[!is.na(presence.test$presence_absence),]
presence.test        <- presence.test[!is.na(presence.test$domain),]


out <- data.frame(domain = 1:9, 
                      presence.unique=as.numeric(table(factor(presence.test$domain[presence.test$presence_absence==1],levels=c(1:9)))), 
                      absence.unique =as.numeric(table(factor(presence.test$domain[presence.test$presence_absence==0],levels=c(1:9)))), 
                      # presence.unique=as.numeric(table(factor(presence.test$domain[presence.test$'E.s'==1 & !duplicated(presence.test$background.bin)],levels=c(1:9)))), 
                      # absence.unique =as.numeric(table(factor(presence.test$domain[presence.test$'E.s'==0 & !duplicated(presence.test$background.bin)],levels=c(1:9)))), 
                      #pseudo  =as.numeric(table(presence.test$domain[presence.test$E.superba==2])), 
                      cbi = NA, 
                      auc = NA)
out <- out[out$presence.unique>0,]

#### CBI
#It is continuous and varies between -1 and +1. 
#Positive values indicate a model which present predictions are 
#consistent with the distribution of presences in the evaluation dataset, 
#values close to zero mean that the model is not different from a random model, 
#negative values indicate counter predictions, i.e., predicting poor quality 
#areas where presences are more frequent (Hirzel et al. 2006).
for(i in out$domain){
  cat("\r",i," of ",9,"  ")
  domain.result <- mask(es_mean, domains[domains$domain==i,])
  
  presence.domain <- presence.test[presence.test$domain==i & presence.test$presence_absence==1,]
  
  out$cbi[out$domain==i] <- ecospat.boyce(values(domain.result)[!is.na(values(domain.result))], presence.domain$predict)$cor
  
  auc_data<-data.frame(Observed = (presence.test$presence_absence[presence.test$domain==i]))
  auc_data<-cbind(plotID=1,auc_data,Predicted1=presence.test$predict[presence.test$domain==i])
  
  out$auc[out$domain==i] <- as.numeric(PresenceAbsence::auc(auc_data,which.model=1)[1])
}
out



eb       <- ecospat.boyce(values(es_mean)[!is.na(values(es_mean))], presence.test$predict[presence.test$presence_absence==1])$cor
auc_data <- data.frame(Observed = presence.test$presence_absence)
auc_data <- cbind(plotID=1,auc_data, Predicted1=presence.test$predict)
au       <- PresenceAbsence::auc(auc_data,which.model=1)

new.line                 <- out[1,]
new.line$domain          <- "all"
new.line$presence.unique <- nrow(presence.test[presence.test$presence_absence==1,])
new.line$absence.unique  <- nrow(presence.test[presence.test$presence_absence==0,])
new.line$cbi             <- eb
new.line$auc             <- au$AUC

out <- rbind(out, new.line)
out
write.csv(out, file = paste0("data/E.superba_CCAMLR domain validation with pre 1970 data.csv"))

