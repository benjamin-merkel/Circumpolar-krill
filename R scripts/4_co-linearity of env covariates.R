#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages

library(sf)
library(concaveman)
library(raster)
library(corrplot)
library(usdm)
library(factoextra)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load mapping data and create polygon for pseudo abscences  -------

sf_use_s2(F)

# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- ("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 

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


model.domain <- circumpolar


# load env covariates
env <-readRDS("data/Environmental covariate stack.rds")

# reduce dimension of ox_0 temp_0 & temp_200 with PCA
# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab10_PCA.html#Today%E2%80%99s_session
env.subset.list <- c("WOA_temp_0", "WOA_temp_200", "WOA_ox_0")
env.pca <- subset(env, env.subset.list)
env.pca <- mask(env.pca, circumpolar)
env.pca <- scale(env.pca)
names(env.pca) <- c("SST","T200","SOX")
# correlation btw temp_0, temp_200 and ox_0
cor <- cor(as.data.frame(env.pca), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))
summary(abs(cor[upper.tri(cor)]))

# PCA btw temp_0, temp_200 and ox_0
clim_samp <- sampleRandom(env.pca, size = 50000)
clim_samp <- as.data.frame(clim_samp)
pca <- prcomp(clim_samp, scale. = TRUE)
clim_pca <- raster::predict(env.pca, pca, index = 1)

# scale env covariates and add PCA axis1
env     <- mask(env, circumpolar)
env     <- scale(env)
env$pca <- clim_pca



png(paste0("figures/PCA of SST-T200-SOX.png"), units="cm",res=500,width=15, height=15)
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()


# set initial selection of relevant covariates
env.init.selected <- subset(env, c(
  "NSIDC_ice_duration",'NSIDC_ice_retreat',"NSIDC_ice_conc","NSIDC_ice_pers","NSIDC_spring_edge",
  "WOA_sal_0","WOA_sal_200",
  "OCCCI_bd","OCCCI_bm",
  "WOA_si_0",
  "bath","dis_1000","dis_ice.pol",
  "Pellichero_ml_depth", 'pca',
  env.subset.list))

names(env.init.selected) <- c("ICEDUR", "ICEMELT", "ICECONC", "ICEPERS", "MIZ", "SSS", "S200", 
                              "BD", "BM", "SSI", "BATH", "DIS","ICEPOL", "MLD", "PCA",
                              "SST","T200","SOX")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test colinearity of covariates across model domain --------

# env.cor.test <- crop(env,model.domain)
env.cor.test <- crop(env.init.selected,model.domain)


# test colinearity of covariates in model domain using VIF
vif3.select <- vifstep(env.cor.test, th = 10)
vif3.select

M  <- cor(as.data.frame(env.cor.test), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))
hc <- hclust(as.dist(1-abs(M)), method = "complete")

png(paste0("figures/Clustering of colinearity of env covariates across model domain.png"), res = 800, width=20, height = 20, units="cm")
opar <- par(mfrow=c(1,1),mar=c(5,4,2,0))
plot(hc,hang=-1,las=1,main='',xlab="",ylab="pearson correlation")
rect.hclust(hc,h=0.2)
# axis(1, at = c(1:16)[hc$labels[hc$order] %in% vif3.select@results$Variables], labels = round(vif3.select@results$VIF, 2), tick = F, las = 2)
par(opar)
dev.off()

png(paste0("figures/Colinearity of env covariates across model domain.png"), res = 800, width=20, height = 18, units="cm")
corrplot.mixed(M, order = 'hclust', upper = "square", tl.pos = "lt", 
               cl.cex = 0.9, tl.cex = 0.9, number.cex = 0.7,
               tl.col = c(1,grey(0.8),1,grey(0.8),rep(1,6),grey(0.8),1,grey(0.8),grey(0.8),rep(1,3),grey(0.8)))
dev.off()

# test colinearity of covariates in model domain using VIF
vif3.select <- vifstep(subset(env.cor.test, names(env.init.selected)[!names(env.init.selected) %in% c("SST","T200","SOX")]), th = 10)
vif3.select


# remove colinear parameters at VIF cutoff th = 10
env.selected <- env.init.selected[[vif3.select@results$Variables]]


png(paste0("figures/Colinearity of selected env covariates across model domain.png"), res = 800, width=20, height = 20, units="cm")
MS  <- cor(as.data.frame(env.selected), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))
corrplot.mixed(MS, order = 'hclust', upper = "square", tl.pos = "lt", cl.cex = 0.9, tl.cex = 0.9, number.cex = 0.7)
dev.off()

# save
saveRDS(env.selected, file = "data/Environmental covariates selected.rds")


