library(raster)

PATH <- "data/files to publish/"

species    <- "E.crystallorophias"
species    <- "E.superba" 

for(species in c("E.crystallorophias", "E.superba")){
  ## load
  ec_mean <- raster(paste0("data/",species,"_circumpolar_GBM_MEAN_ensemble_ROC_weighted_Dec-Mar.tif"))
  ec_sd   <- raster(paste0("data/",species,"_circumpolar_GBM_SD_ensemble_unweighted_Dec-Mar.tif"))
  ec_cv   <- raster(paste0("data/",species,"_circumpolar_GBM_CV_ensemble_ROC_weighted_Dec-Mar.tif"))
  ec_en_eval <- readRDS(paste0("data/",species,"_circumpolar_GBM_ensemble model evaluation_Dec-Mar.rds"))
  
  ## binary habitat
  ec_tss_cutoff <- ec_en_eval[,2,8]
  ec_bin <- ec_mean
  ec_bin[ec_bin< ec_tss_cutoff] <- 0
  ec_bin[ec_bin>0] <- 2
  
  writeRaster(ec_mean, paste0(PATH, species, "_mean_prediction.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(ec_sd,   paste0(PATH, species, "_standard_deviation.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(ec_cv,   paste0(PATH, species, "_coefficient_of_variation.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(ec_bin,  paste0(PATH, species, "_binary_habitat.tif"), format="GTiff", overwrite=TRUE)
}
