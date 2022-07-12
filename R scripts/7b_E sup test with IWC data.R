library(terra)
library(sf)
library(sp)
library(raster)


sf_use_s2(F)

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

domains <- readRDS("data/map data/Marine Protected Area planning domains.rds")
domains <- st_transform(domains, south_pole_equal_area.proj)

ccamlr_domains      <- st_read("data/CCAMLR/bm_mpa_planningDomains.shp")
ccamlr_domains      <- st_transform(ccamlr_domains, south_pole_equal_area.proj)
ccamlr_domains_crop <- st_intersection(ccamlr_domains, circumpolar)
rownames(ccamlr_domains) <- ccamlr_domains$OBJECTID <- ccamlr_domains$Index


env       <- readRDS("data/Environmental covariate stack.rds")


path_used <- "data/IWC/v7.1"
file.list <- list.files(path_used, pattern = ".csv", recursive = T)
file.list <- file.list[!grepl("SHL.csv", file.list)]
# These are coded as follows:
#   01 Pilot      	06 Sperm 	        11 Right 	          16 Pygmy Right 
#   02 Bottlenose 	07 Humpback 	    12 Gray 	          17 Cuvier's Beaked 
# 	03 Killer     	08 Sei 	          13 Baird's Beaked 	18 Bowhead 
#   04 Blue 	      09 Common Minke 	14 Baleen          	19 Beaked (unspecified)
#   05 Fin 	        10 Bryde's 	      15 Pygmy Blue 	    20 Antarctic Minke

columns_used <- c("Day" ,"Mon", "Year", "Sp", "Len" , "L.u" , "Sx", "NoF",                                                   
                  "F1.L", "F1.S", "F.u", "Lat", "Mn", "X", "Ac", "Lon",                                                   
                  "Mn.1", "X.1", "Ac.1", "Nt", "Fem", "Mat", "Txt")

for(i in 1:length(file.list)){
  dat1 <- read.csv(paste(path_used, file.list[i],sep="/"))
  dat1 <- dat1[,columns_used]
  if(i==1) data <- dat1 else data <- rbind(dat1,data)
}
data                            <- data[complete.cases(data),]
data$date                       <- ISOdate(data$Year,data$Mon,data$Day)
data$doy                        <- as.numeric(strftime(data$date, "%j")) 
data$latitude                   <- data$Lat + data$Mn/60
data$latitude[data$X == "S"]    <- - data$latitude[data$X == "S"]
data$longitude                  <- data$Lon + data$Mn.1/60
data$longitude[data$X.1 == "W"] <- - data$longitude[data$X.1 == "W"]
data                            <- data[data$latitude <= -30,]
data2                           <- st_as_sf(data, coords = c("longitude","latitude"), crs = 4326)
data2$count                     <- 1


# number of digits after dot for locations
decimalplaces <- function(x) {
  ifelse(abs(x - round(x)) > .Machine$double.eps^0.5,
         nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(x)))),
         0)
}

data$digits.lon <- decimalplaces(data$longitude)
data$digits.lat <- decimalplaces(data$latitude)

round(tapply(data$digits.lon, data$Year, mean),1)
round(tapply(data$digits.lat, data$Year, mean),1)





background_raster_lonlat         <- raster(xmn = -180, xmx=180, ymn=-90, ymx=-40, res=1, crs = 4326)
values(background_raster_lonlat) <- 1:ncell(background_raster_lonlat)


# # mean timimng of whale catch 
data3 <- data2[data2$Sp %in% c(4,5,20) & data2$Lat >= 40 & data2$X == "S",]
data.doy <- data3[!is.na(data3$doy),] #  & data3$Mon %in% c(12,1:4)
data.doy$doy2 <- data.doy$doy
data.doy$doy2[data.doy$doy2>200] <- data.doy$doy2[data.doy$doy2>200] - 365


data.doy2 <- st_as_sf(data.doy[data.doy$Sp %in% c(4,5,20)  & data.doy$Mon %in% c(1:3,12),], coords = c("longitude","latitude"), crs=lonlat.proj)
ras_doy <- rasterize(data.doy2, background_raster_lonlat, field = "doy2",  fun = mean)
pol_doy <- rasterToPolygons(ras_doy)
pol_doy <- st_transform(st_as_sf(pol_doy), south_pole_equal_area.proj)


png("figures/IWC historical catch timing Dec-Mar.png",units="cm",res=800,width=20, height=20)
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(st_geometry(circumpolar))
plot(st_intersection(pol_doy, circumpolar), border="transparent",add=T)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(ras_doy,legend.only=T,col=bpy.colors(n = 11, cutoff.tails = 0.1, alpha = 1.0),add=T,
     axis.args=list(cex.axis=0.8),
     legend.args=list(text='days since 1 Jan', side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
par(opar)
dev.off()


# whale catches Dec-Mar
data3                   <- data2[data2$Sp %in% c(4,5,20) & data2$Mon %in% c(1:3,12) & data2$Lat >= 40 & data2$X == "S",]
dat_ras                 <- rasterize(data3, background_raster_lonlat, field = "count",  fun = "sum")
dat_ras[is.na(dat_ras)] <- 0
pol_ras                 <- rasterToPolygons(dat_ras)
pol_ras                 <- st_transform(st_as_sf(pol_ras), south_pole_equal_area.proj)
pol_ras                 <- st_crop(pol_ras, circumpolar)
pol_ras                 <- st_intersection(pol_ras, circumpolar)
pol_ras$log <- log(pol_ras$layer)
pol_ras$log[pol_ras$log == -Inf] <- NA
pol_ras$catch <- pol_ras$layer
pol_ras$catch[pol_ras$catch>1000] <- 1000
dat_ras2 <- dat_ras
dat_ras2[dat_ras2>1000] <- 1000


png("figures/IWC historical catch Dec-Mar.png",units="cm", res=800, width=20, height=20)
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(st_geometry(circumpolar))
cols_used <- colorRampPalette(c("white","red","purple"))(41)
plot(st_geometry(pol_ras), border="transparent", col=cols_used[(1+round(pol_ras$catch/25,0))],add=T)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(dat_ras2,legend.only=T, col=cols_used,
     axis.args=list(cex.axis=0.8, at=seq(0,1000,200), labels=c(seq(0,800,200),"\u22651000")),
     legend.args=list(text='whales caught', side=4, font=2, line=3, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
par(opar)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spatial correlation in lon lat ------------

dat_ras[dat_ras==0 & is.na(es_mean_lonlat)] <- NA
dat_ras             <- crop(dat_ras, extent(-180, 180, -90, -50))
es_mean             <- raster(paste0("data/E.superba_circumpolar_GBM_MEDIAN_3000_ensemble_TSS_weighted_Dec-Mar.tif"))
es_mean_lonlat      <- projectRaster(es_mean, dat_ras)
spat_cor_lonlat     <- raster::corLocal(es_mean_lonlat, dat_ras, ngb = c(11,21), type = "spearman")
spat_cor2           <- spat_cor_lonlat
spat_cor2[spat_cor2>1 | spat_cor2 < -1] <- NA
spat_cor3           <- rasterToPolygons((spat_cor2))
spat_cor3           <- st_transform(st_as_sf(spat_cor3), south_pole_equal_area.proj)
spat_cor3           <- st_crop(spat_cor3, circumpolar)
spat_cor3           <- st_intersection(spat_cor3, circumpolar)


png("figures/E superba vs IWC spatial correlation Dec-Mar.png",units="cm",res=800,width=20, height=20)
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(st_geometry(circumpolar))
plot(spat_cor3, border="transparent", breaks =seq(-1,1,0.4), pal=brewer.pal(5,"Spectral"),add=T)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot((spat_cor2),legend.only=T,
     col=brewer.pal(5,"Spectral"),breaks=seq(-1,1,0.4),add=T,
     axis.args=list(cex.axis=0.8),
     legend.args=list(text='correlation', side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
par(opar)
dev.off()


png("figures/E.superba GBM mean prediction with ice duration.png", res = 800,   width=20, height = 20, units="cm")
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(st_geometry(model.domain),lty=2)  
plot(mask(es_mean, model.domain),zlim=c(0,1000), add=T, legend =F)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(es_mean,legend.only = T, zlim=c(0,1000), 
     axis.args=list(at=c(0,1000), labels= c(0, 1),
                    cex.axis=1,tick=F,line=-0.2),
     legend.args=list(text="probability", side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
par(opar)
dev.off()

# # test catch timing vs spatial correlation
# dat <- data.frame(doy=values(ras_doy),
#                   cor=values(spat_cor2))
# names(dat) <- c("doy","cor")
# dat <- dat[complete.cases(dat),]
# cor(dat)
# 
# plot(dat, ylim=c(-1,1))
# 
# library(mgcv)
# g1 <-gam(cor~s(doy), data=dat, family="gaussian")
# l1 <-lm(cor~doy, data=scale(dat))
# summarY(l1)
# plot(g1)
# 
# plot(dat$doy, dat$cor)
# lines(-100:150, predict(g1, newdata= data.frame(doy = -100:150), type="response"),col=2,lwd=2)
# 
# 
# ggplot(dat, aes(x = doy, y = cor) ) +
#   geom_point() +
#   geom_smooth(method = "gam",alpha = .3, aes())
# 


## correlation between E sup GBM and IWC within different domains

pts        <- as.points(rast(es_mean_lonlat))
names(pts) <- "pred"
pts        <- st_as_sf((pts), coords = c('x','y'),crs=4326)
pts$iwc    <- extract(dat_ras, as_Spatial(pts))
pts        <- pts[!is.na(pts$pred) & !is.na(pts$iwc),]
pts        <- st_transform(pts, st_crs(domains))

pts$domain <- as.character(st_intersects(pts, domains))
pts$domain <- as.numeric(substr(pts$domain,1,1))
pts$pred   <- data.frame(pts[,1])[,1]
pts1        <- pts[!is.na(pts$iwc),]

iwc_cor <- data.frame(domain=c(1:9,"all"),
                      spearman=NA,
                      pearson=NA)
for(i in 1:9){
  cat("\r",i," of ",9,"  ")
  
  pts2 <- pts1[pts1$domain==i,]
  pts2 <- pts2[!is.na(pts2$domain),]
  
  iwc_cor$spearman[i] <- cor(pts2$pred, pts2$iwc,     method = "spearman")
  iwc_cor$pearson[i]  <- cor(pts2$pred, pts2$iwc,     method = "pearson")
  
}  
iwc_cor$spearman[10] = cor(pts1$pred, pts1$iwc, method = "spearman")
iwc_cor$pearson[10] = cor(pts1$pred, pts1$iwc, method = "pearson")

iwc_cor
write.csv(iwc_cor, file = paste0("data/E.superba vs IWC by CCAMLR domain.csv"))





# 4 panel figure

png("figures/E.superba with IWC 4 panel figure.png", res = 800, width=40, height = 40, units="cm")
opar <- par(mfrow=c(2,2), mar=c(0,0,0,0))

plot(st_geometry(model.domain),lty=2)  
plot(mask(es_mean, model.domain),zlim=c(0,1000), add=T, legend =F)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(es_mean,legend.only = T, zlim=c(0,1000), 
     axis.args=list(at=c(0,1000), labels= c(0, 1),
                    cex.axis=1,tick=F,line=-0.2),
     legend.args=list(text="probability", side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
mtext('a)', side=3, line= -2.5, at=-4350, cex=2.5)


plot(st_geometry(circumpolar))
cols_used <- colorRampPalette(c("white","red","purple"))(41)
plot(st_geometry(pol_ras), border="transparent", col=cols_used[(1+round(pol_ras$catch/25,0))],add=T)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(dat_ras2,legend.only=T, col=cols_used,
     axis.args=list(cex.axis=0.8, at=seq(0,1000,200), labels=c(seq(0,800,200),"\u22651000")),
     legend.args=list(text='whales caught', side=4, font=2, line=3, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
mtext('b)', side=3, line= -2.5, at=-4350, cex=2.5)


plot(st_geometry(circumpolar))
plot(spat_cor3, border="transparent", breaks =seq(-1,1,0.4), pal=brewer.pal(5,"Spectral"),add=T)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot((spat_cor2),legend.only=T,
     col=brewer.pal(5,"Spectral"),breaks=seq(-1,1,0.4),add=T,
     axis.args=list(cex.axis=0.8),
     legend.args=list(text='correlation', side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
mtext('c)', side=3, line= -2.5, at=-4350, cex=2.5)


plot(st_geometry(circumpolar))
plot(st_intersection(pol_doy, circumpolar), border="transparent",add=T)
contour(rast(env$NSIDC_ice_duration), levels = c(10/365,9/12,11/12)*365, lty = c(3,2,1), add=T, lwd = 2, drawlabels = F)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(ras_doy,legend.only=T,col=bpy.colors(n = 11, cutoff.tails = 0.1, alpha = 1.0),add=T,
     axis.args=list(cex.axis=0.8),
     legend.args=list(text='days since 1 Jan', side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
mtext('d)', side=3, line= -2.5, at=-4350, cex=2.5)

par(opar)
dev.off()