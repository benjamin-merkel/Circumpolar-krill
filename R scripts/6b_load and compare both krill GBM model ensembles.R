library(raster)
library(sf)
library(RColorBrewer)
library(viridis)
library(spatialEco)
library(concaveman)
library(cartography)
library(dplyr)
library(scales)

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


PF_fronts    <- data.frame(lon = lonPF, lat = latPF, name = "Polar Front (PF)")
SACCF_fronts <- data.frame(lon = lonSACCF, lat = latSACCF, name = "Southern Antarctic Circumpolar Current Front (SACCf)")
SAF_fronts   <- data.frame(lon = lonSAF, lat = latSAF, name = "Subantarctic Front (SAF)")
SB_fronts    <- data.frame(lon = lonSB, lat = latSB, name = "Southern Boundary (SB)")

PF_fronts <- st_as_sf(PF_fronts[!is.na(PF_fronts$lat),], coords = c('lon','lat'), crs =4326, agr = "name")
PF_fronts <- st_transform(st_combine(PF_fronts), south_pole_equal_area.proj)
PF_fronts <- st_cast(PF_fronts, "MULTILINESTRING")
SACCF_fronts <- st_as_sf(SACCF_fronts[!is.na(SACCF_fronts$lat),], coords = c('lon','lat'), crs =4326, agr = "name")
SACCF_fronts <- st_transform(st_combine(SACCF_fronts), south_pole_equal_area.proj)
SACCF_fronts <- st_cast(SACCF_fronts, "MULTILINESTRING")
SAF_fronts <- st_as_sf(SAF_fronts[!is.na(SAF_fronts$lat),], coords = c('lon','lat'), crs =4326, agr = "name")
SAF_fronts <- st_transform(st_combine(SAF_fronts), south_pole_equal_area.proj)
SAF_fronts <- st_cast(SAF_fronts, "MULTILINESTRING")
SB_fronts <- st_as_sf(SB_fronts[!is.na(SB_fronts$lat),], coords = c('lon','lat'), crs =4326, agr = "name")
SB_fronts <- st_transform(st_combine(SB_fronts), south_pole_equal_area.proj)
SB_fronts <- st_cast(SB_fronts, "MULTILINESTRING")


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
ec_rast <- readRDS(paste0("data/",species,"_circumpolar_GBM_10raster_predictions_Dec-Mar.rds"))
ec_dat  <- readRDS(paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
ec_eval <- readRDS(paste0("data/",species,"_circumpolar_GBM_100_ensemble model evaluation_Dec-Mar.rds"))


# icer <- crop(icer, ec_mean)
# ec_mean[is.na(icer)]<-NA
# ec_sd[is.na(icer)]<-NA

species                 <- "E.superba" 
es_mean <- raster(paste0("data/",species,"_circumpolar_GBM_MEDIAN_3000_ensemble_TSS_weighted_Dec-Mar.tif"))
es_sd   <- raster(paste0("data/",species,"_circumpolar_GBM_SD_3000_ensemble_Dec-Mar.tif"))
es_rast <- readRDS(paste0("data/",species,"_circumpolar_GBM_10raster_predictions_Dec-Mar.rds"))
es_dat  <- readRDS(paste0("data/",species,"_circumpolar_response_data_Dec-Mar.rds"))
es_eval <- readRDS(paste0("data/",species,"_circumpolar_GBM_100_ensemble model evaluation_Dec-Mar.rds"))





## binary habitat
es_tss_cutoff <- mean(es_eval[2,2,1,,1]) 
ec_tss_cutoff <- mean(ec_eval[2,2,1,,1]) 

es_bin <- es_mean
es_bin[es_bin< es_tss_cutoff] <- 0
es_bin[es_bin>0] <- 1

ec_bin <- ec_mean
ec_bin[ec_bin< ec_tss_cutoff] <- 0
ec_bin[ec_bin>0] <- 2

e_overlap <- sum(stack(es_bin,crop(ec_bin,es_bin)),na.rm=T)
e_overlap[e_overlap == 0 & is.na(es_mean)] <- NA

e_data <- e_overlap
e_data[!is.na(e_data)] <- 1

es_pol <- rasterToPolygons(es_bin, fun = function(x) {x==1}, n = 16, dissolve = T)
ec_pol <- rasterToPolygons(ec_bin, fun = function(x) {x==2}, n = 16, dissolve = T)
# names(es_pol) <- names(ec_pol) <- "habitat"
# st_write(st_as_sf(es_pol), "data/Antarctic krill habitat.shp")
# st_write(st_as_sf(ec_pol), "data/Ice krill habitat.shp")

es_rast_bin <- ec_rast_bin <- vector(mode="list")
for(i in 1:10){
  es_rast_bin[[i]] <- es_rast[[1]]
  es_rast_bin[[i]][es_rast_bin[[i]] < es_eval[2,2,1,,1][i]] <- 0
  es_rast_bin[[i]][es_rast_bin[[i]] > 0] <- 1
  
  ec_rast_bin[[i]] <- ec_rast[[1]]
  ec_rast_bin[[i]][ec_rast_bin[[i]] < ec_eval[2,2,1,,1][i]] <- 0
  ec_rast_bin[[i]][ec_rast_bin[[i]] > 0] <- 1
}
es_rast_bin <- stack(es_rast_bin)
ec_rast_bin <- stack(ec_rast_bin)
ec_rast_bin[is.na(icer)]<-NA



patterns <- c("diamond","grid","hexagon","horizontal", "vertical",
              "zigzag","left2right","right2left","circle")

es_pol_hatch <- st_as_sf(es_pol) %>%
  hatchedLayer(mode = "sfc", pattern = patterns[7], density = 15)
es_pol_hatch <- st_sf(geometry = es_pol_hatch)

ec_pol_hatch <- st_as_sf(ec_pol) %>%
  hatchedLayer(mode = "sfc", pattern = patterns[8], density = 15)
ec_pol_hatch <- st_sf(geometry = ec_pol_hatch)


png("figures/krill comparison map hatched 2.png",units="cm",res=800,width=20, height=20)
cols_used <- c(grey(0.9),brewer.pal(5,"Set3")[c(2,5,4)])
cols_used <- c(grey(0.9),brewer.pal(5,"Set3")[c(4,5,2)])
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(circumpolar)
plot(e_data,  col = grey(0.95), axes=F,legend=F,add=T)
plot(es_pol_hatch$geometry, col = cols_used[2], add=T, lwd =0.5)
plot(ec_pol_hatch$geometry, col = cols_used[3], add=T, lwd =0.5)
plot(es_pol, border = cols_used[2], add=T, lwd =0.5)
plot(ec_pol, border = cols_used[3], add=T, lwd =0.5)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(circumpolar, add=T)
# plot(grat.40s,add=T,lwd=1.5)
legend("bottomright",legend=c("No data","Antarctic krill","Ice krill"),pch=22,
       pt.bg = c("white",cols_used[2:3]), col = 1, box.col = "transparent", cex = 1.3)
par(opar)
dev.off()



# png("figures/circumpolar/krill comparison map with CCAMLR fish.png",units="cm",res=800,width=20, height=20)
# cols_used <- c(grey(0.9),brewer.pal(5,"Set3")[c(2,5,4)])
# opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
# plot(grat.50S,lwd=1.5)
# plot(e_overlap,  col = cols_used, axes=F,legend=F,add=T)
# plot(ccamlr_fish$geometry ,add=T,lty=2)
# # plot(ccamlr_mpa$geometry ,add=T, border="darkred",lty=3)
# plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
# plot(st_geometry(land50),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
# # plot(st_centroid(ccamlr_fish)$geometry ,add=T, border="darkred",lty=3,pch=as.character(ccamlr_fish$ShortLabel),cex=2)
# plot(grat.50S,lwd=1.5, add=T)
# par(opar)
# dev.off()



png("figures/krill comparison map with CCAMLR mpa.png",units="cm",res=800,width=20, height=20)
cols_used <- c(grey(0.9),brewer.pal(5,"Set3")[c(4,5,1)])
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(circumpolar)
# plot(e_overlap,  col = cols_used, axes=F,legend=F,add=T)
plot(e_data,  col = grey(0.95), axes=F,legend=F,add=T)
plot(es_pol_hatch$geometry, col = cols_used[2], add=T, lwd =0.5)
plot(ec_pol_hatch$geometry, col = cols_used[3], add=T, lwd =0.5)
plot(es_pol, border = cols_used[2], add=T, lwd =0.5)
plot(ec_pol, border = cols_used[3], add=T, lwd =0.5)

plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(circumpolar, add=T)
plot(st_geometry(st_point_on_surface(ccamlr_domains_crop)), add=T, col=grey(0.2), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1.5)

legend("bottomright",legend=c("No data","Antarctic krill","Ice krill"),pch=22,
       pt.bg = c("white",cols_used[2:3]), col = 1, box.col = "transparent", cex = 1.3)
par(opar)
dev.off()




png("figures/krill available data with CCAMLR mpa.png",units="cm",res=800,width=20, height=20)
cols_used <- c(grey(0.9),brewer.pal(5,"Set3")[c(4,5,1)])
opar <- par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(circumpolar)
plot(es_dat,  col = cols_used[2], pch=3, axes=F,legend=F,add=T, cex=0.5)
plot(ec_dat,  col = cols_used[3], pch=4, axes=F,legend=F,add=T, cex=0.5)
plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(circumpolar, add=T)
plot(st_geometry(st_point_on_surface(ccamlr_domains_crop)), add=T, col=grey(0.2), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1.5)
# legend("bottomright",legend=c("Ant krill","Ice krill"),pch=21,
#        pt.bg = c(cols_used[2:3]), col = "white", box.col = "transparent")
par(opar)
dev.off()





# percent of total area south of 50S
tot_area  <- sum(as.numeric(table(values(e_overlap))[1:4]))
over_area <- as.numeric(table(values(e_overlap))[4])
esup_area <- as.numeric(table(values(es_bin))[2])
ecry_area <- as.numeric(table(values(ec_bin))[2])

stats <- data.frame(species    = c("Esup","Ecry","total","overlap"),
                    area.miokm = c(esup_area,ecry_area,tot_area,over_area)*100/1000000,
                    per.total  = c(esup_area,ecry_area,tot_area,over_area)/tot_area,
                    per.over   = over_area/c(esup_area,ecry_area,tot_area,over_area))
stats[,2:4] <- round(stats[,2:4], 3)
stats


# precent of total krill habitat area
tot_area  <- sum(as.numeric(table(values(e_overlap))[2:4]))
over_area <- as.numeric(table(values(e_overlap))[4])
esup_area <- as.numeric(table(values(es_bin))[2])
ecry_area <- as.numeric(table(values(ec_bin))[2])

stats <- data.frame(species    = c("Esup","Ecry","total","overlap"),
                    area.miokm = c(esup_area,ecry_area,tot_area,over_area)*100/1000000,
                    per.total  = c(esup_area,ecry_area,tot_area,over_area)/tot_area,
                    per.over   = over_area/c(esup_area,ecry_area,tot_area,over_area))
stats[,2:4] <- round(stats[,2:4], 3)
stats



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## split by on and off shelf -------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


e_overlap_shelf <- e_overlap_deep <- e_overlap
e_overlap_shelf[shelf==1] <- NA
e_overlap_deep [shelf==0] <- NA


total<- as.numeric(table(values(e_overlap)))
deep <- as.numeric(table(values(e_overlap_deep)))
shelf<- as.numeric(table(values(e_overlap_shelf)))
stats<- data.frame(what = c("no", "Esup","Ecry","both"),
                   total, shelf, deep)
stats[2,2:4] <- stats[2,2:4] + stats[4,2:4]
stats[3,2:4] <- stats[3,2:4] + stats[4,2:4]
# stats <- stats[1:3,]
stats$rel_shelf <- stats$shelf/stats$total
stats$rel_deep  <- stats$deep /stats$total

stats


png("figures/krill proportion of habitat on and off shelf.png",units="cm", res=800, width=8, height=18)

opar <- par(mfrow=c(1,1), mar=c(3,5,2,1))
cols_bar <- c(brewer.pal(6,"Accent"))[4:5]
barplot(t(as.matrix(stats[2:4,5:6])),las=1,xaxt="n",ylab="Proportion of habitat", col=cols_bar, border=1)
axis(4, at =c(0.25, 0.8), labels = c("SHELF", "OPEN OCEAN") ,tick=F, line = -1)
axis(1, at =c(0.7, 1.9, 3.1), labels = c("Antarctic","Ice","area of"),tick=F, line = -0.5)
axis(1, at =c(0.7, 1.9, 3.1), labels = c("krill","krill","overlap"),tick=F, line = 0.5)
par(opar)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## split by MPA planning domain -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dom.out <- data.frame(domain = 1:9, 
                      e.sup.area = 0,
                      e.cry.area = 0,
                      over.area  = 0,
                      total.area = 0)
for(i in 1:9){
  cat("\r",i," of ",9,"  ")
  domain.result <- mask(e_overlap,domains[domains$domain==i,])
  dom.out$total.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[2:4]))*100/1000000
  dom.out$over.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[4]))*100/1000000
  dom.out$e.sup.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[2]))*100/1000000
  dom.out$e.cry.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[3]))*100/1000000
}
dom.out[is.na(dom.out)]<-0
dom.out
apply(dom.out,2,sum)

round(dom.out$e.sup.area/sum(dom.out$e.sup.area),2)


dom.out3 <- t(as.matrix(dom.out)[,c(2,4,3)])
colnames(dom.out3)<-1:9

png("figures/krill habitat area by planning area v2.png",units="cm",res=800,width=17, height=10)
par(mfrow=c(1,2),mar=c(2,5,0.5,0.5))
cols_used <- c(brewer.pal(5,"Set3")[c(2,4,5)])
cols_used <- c(brewer.pal(5,"Set3")[c(4,2,5)])
barplot(dom.out3,col=cols_used, border=grey(0.2), 
        ylab=expression(paste(10^6," ", km^2,sep="")),
        main="",las=1)
par(mar=c(0.5,0.5,4,0.5))
plot(st_geometry(circumpolar))
plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3, lwd=1.5)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(circumpolar, add=T)
plot(st_geometry(st_point_on_surface(ccamlr_domains_crop)), add=T, col=grey(0.2), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1.5)
legend("topleft" ,pch=22,pt.cex=3,col=grey(0.2),pt.bg=rev(cols_used), inset=c(0,-0.23), xpd=TRUE, 
       legend=c("Ice krill","Overlap","Antarctic krill"), box.col ="transparent")

par(opar)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## split by fishing boxes -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# dom.out <- data.frame(sub.area = ccamlr_fish$ShortLabel, 
#                       e.sup.area = 0,
#                       e.cry.area = 0,
#                       over.area  = 0,
#                       total.area = 0)
# for(i in 1:nrow(dom.out)){
#   cat("\r",i," of ",19,"  ")
#   domain.result <- mask(e_overlap,ccamlr_fish[ccamlr_fish$ShortLabel==ccamlr_fish$ShortLabel[i],])
#   dom.out$total.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[2:4]))*100/1000000
#   dom.out$over.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[4]))*100/1000000
#   dom.out$e.sup.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[2]))*100/1000000
#   dom.out$e.cry.area[i] <- sum(as.numeric(table(factor(values(domain.result),levels=c(0:3)))[3]))*100/1000000
# }
# dom.out[is.na(dom.out)]<-0
# dom.out$area <- substr(dom.out$sub.area,1,2)
# dom.out <- dom.out[order(dom.out$sub.area),]
# 
# #Ant krill stats
# test <- dom.out[dom.out$area==48,2]+dom.out[dom.out$area==48,4]
# round(test/sum(test),2)
# test <- dom.out[dom.out$area==58,2]+dom.out[dom.out$area==58,4]
# round(test/sum(test),2)
# test <- dom.out[dom.out$area==88,2]+dom.out[dom.out$area==88,4]
# round(test/sum(test),2)
# 
# 
# #Ice krill stats
# test <- dom.out[dom.out$area==48,3]+dom.out[dom.out$area==48,4]
# round(test/sum(test),2)
# test <- dom.out[dom.out$area==58,3]+dom.out[dom.out$area==58,4]
# round(test/sum(test),2)
# test <- dom.out[dom.out$area==88,3]+dom.out[dom.out$area==88,4]
# round(test/sum(test),2)
# 
# 
# dom.out3 <- t(as.matrix(dom.out)[,c(2,4,3)])
# colnames(dom.out3)<- dom.out$sub.area
# 
# png("figures/circumpolar/krill habitat area by fishing area.png",units="cm",res=800,width=20, height=15)
# par(mfrow=c(1,1),mar=c(4,5,1,1))
# cols_used <- c(brewer.pal(5,"Set3")[c(2,4,5)])
# barplot(dom.out3,col=cols_used, border="white",
#         ylab="10^6 km^2",main="",las=2)
# # legend("topright",pch=20,pt.cex=3,col=cols_used, legend=c("E.cry","overlap","E.sup"))
# dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## split by longitude and latitude sections ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lons5  <- data.frame(x1 = c(seq(-70,175,5),seq(-180,-65,5)),
                     x2 = c(seq(-75,180,5),seq(-175,-70,5)))
lats1 <- -90
lats2 <- -50

lons <- seq(-180, 180, 0.01)
# lons2 <- 180

### split into 5 degree  longitude chunks
poly5 <- lapply(X = 1:nrow(lons5),
                FUN = function(x) st_polygon(list(rbind(c(lons5[x, 1], lats1),
                                                        c(lons5[x, 2], lats1),
                                                        c(lons5[x, 2], lats2),
                                                        c(lons5[x, 1], lats2),
                                                        c(lons5[x, 1], lats1)))))
poly5 <- st_sfc(poly5,crs = 4326)
poly5 <- st_transform(poly5, south_pole_equal_area.proj)


poly.out <- data.frame(poly = 1:72,
                       e.sup.area = 0,
                       e.sup.area.min = 0,
                       e.sup.area.max = 0,
                       e.sup.area.sd = 0,
                       e.cry.area = 0,
                       e.cry.area.min = 0,
                       e.cry.area.max = 0,
                       e.cry.area.sd = 0,
                       over.area  = 0,
                       total.area = 0)

for(i in 1:72){
  cat("\r",i," of ",72,"  ")
  poly.result <- mask(e_overlap, as_Spatial(st_zm(poly5[i])))
  poly.out$total.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[2:4]))*100/1000000
  poly.out$over.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[4]))*100/1000000
  poly.out$e.sup.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[2]))*100/1000000
  poly.out$e.cry.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[3]))*100/1000000
  
  for(j in 1:10) {
    rr <- sum(as.numeric(table(factor(values(mask(es_rast_bin[[j]], as_Spatial(st_zm(poly5[i])))),levels=c(0:1)))[2]))*100/1000000
    if(j == 1) es_rr <- rr else es_rr <- c(es_rr, rr)
  }
  for(j in 1:10) {
    rr <- sum(as.numeric(table(factor(values(mask(ec_rast_bin[[j]], as_Spatial(st_zm(poly5[i])))),levels=c(0:1)))[2]))*100/1000000
    if(j == 1) ec_rr <- rr else ec_rr <- c(ec_rr, rr)
  }
  
  poly.out$e.sup.area.min[i] <- min(es_rr)
  poly.out$e.sup.area.max[i] <- max(es_rr)
  poly.out$e.sup.area.sd[i]  <- sd(es_rr)
  
  poly.out$e.cry.area.min[i] <- min(ec_rr)
  poly.out$e.cry.area.max[i] <- max(ec_rr)
  poly.out$e.cry.area.sd[i]  <- sd(ec_rr)
  
}
poly.out[is.na(poly.out)]<-0
poly.out$e.sup.area.se <- poly.out$e.sup.area.sd/sqrt(10)
poly.out$e.cry.area.se <- poly.out$e.cry.area.sd/sqrt(10)

# relative use of each chunk
poly.out2     <- poly.out
poly.out2$e.sup.area <- poly.out2$e.sup.area + poly.out2$over.area
poly.out2$e.cry.area <- poly.out2$e.cry.area + poly.out2$over.area
poly.out2$e.sup.rel  <- poly.out2$e.sup.area/sum(poly.out2$e.sup.area)
poly.out2$e.cry.rel  <- poly.out2$e.cry.area/sum(poly.out2$e.cry.area)
poly.out2$over.rel   <- poly.out2$over.area/sum(poly.out2$over.area)

poly.out2$e.sup.rel.uci  <- (poly.out2$e.sup.area + (1.96 * poly.out2$e.sup.area.se))/sum(poly.out2$e.sup.area)
poly.out2$e.sup.rel.lci  <- (poly.out2$e.sup.area - (1.96 * poly.out2$e.sup.area.se))/sum(poly.out2$e.sup.area)
poly.out2$e.cry.rel.uci  <- (poly.out2$e.cry.area + (1.96 * poly.out2$e.cry.area.se))/sum(poly.out2$e.cry.area)
poly.out2$e.cry.rel.lci  <- (poly.out2$e.cry.area - (1.96 * poly.out2$e.cry.area.se))/sum(poly.out2$e.cry.area)
poly.out2$e.cry.rel.lci[poly.out2$e.cry.rel.lci < 0] <- 0

poly.out2$e.sup.min.rel  <- poly.out2$e.sup.area.min/sum(poly.out2$e.sup.area)
poly.out2$e.sup.max.rel  <- poly.out2$e.sup.area.max/sum(poly.out2$e.sup.area)
poly.out2$e.cry.min.rel  <- poly.out2$e.cry.area.min/sum(poly.out2$e.cry.area)
poly.out2$e.cry.max.rel  <- poly.out2$e.cry.area.max/sum(poly.out2$e.cry.area)

# poly.out2 <- poly.out2[,c(2,3,4)]
# names(poly.out2) <- c("E superba","E cryst","overlap")


lons5.lab <- lons5
lons5.lab$x1[lons5.lab$x1 > 0] <- paste0(lons5.lab$x1[lons5.lab$x1 > 0], "E") 
lons5.lab$x1[lons5.lab$x1 < 0] <- paste0(lons5.lab$x1[lons5.lab$x1 < 0], "W") 
lons5.lab$x1 <- gsub("-","",lons5.lab$x1)
lons5.lab$x2[lons5.lab$x2 > 0] <- paste0(lons5.lab$x2[lons5.lab$x2 > 0], "E") 
lons5.lab$x2[lons5.lab$x2 < 0] <- paste0(lons5.lab$x2[lons5.lab$x2 < 0], "W") 
lons5.lab$x2 <- gsub("-","",lons5.lab$x2)


png("figures/krill proportion of habitat longitudinal around the continent smoothed 2.png",units="cm",
    res=800, width=20, height=12)

opar <- par(mfrow=c(1,1), mar=c(3,5,2,1))
col_used <- c(brewer.pal(5,"Set3")[c(4,5,2)])
plot(movingFun(poly.out2$e.sup.rel,6,circular = T),type="l",
     xaxt="n",las=1,xlab="",ylab="Proportion of habitat",ylim=c(0,0.06),col=col_used[1],lwd=3)
# polygon(x = c(1:72, rev(1:72)), y = c(movingFun(poly.out2$e.sup.min.rel,6,circular = T), rev(movingFun(poly.out2$e.sup.max.rel,6,circular = T))),
#         col = alpha(col_used[1], 0.2), border= "transparent")
polygon(x = c(1:72, rev(1:72)), y = c(movingFun(poly.out2$e.sup.rel.lci,6,circular = T), rev(movingFun(poly.out2$e.sup.rel.uci,6,circular = T))),
        col = alpha(col_used[1], 0.2), border= "transparent")
lines(movingFun(poly.out2$e.sup.rel,6,circular = T),type="l",col=col_used[1],lwd=3)

# polygon(x = c(1:72, rev(1:72)), y = c(movingFun(poly.out2$e.cry.min.rel,6,circular = T), rev(movingFun(poly.out2$e.cry.max.rel,6,circular = T))),
#         col = alpha(col_used[2], 0.2), border= "transparent")
polygon(x = c(1:72, rev(1:72)), y = c(movingFun(poly.out2$e.cry.rel.lci,6,circular = T), rev(movingFun(poly.out2$e.cry.rel.uci,6,circular = T))),
        col = alpha(col_used[2], 0.2), border= "transparent")
lines(movingFun(poly.out2$e.cry.rel,6,circular = T),type="l",col=col_used[2],lwd=3)
# lines(movingFun(poly.out2[,3],6,circular = T),type="l",col=col_used[3],lwd=3)
abline(v=c(21, 45),lty=3)
axis(1,at=seq(1,74,4),labels = lons5.lab[,1][seq(1,74,4)],las=2)
# invisible(lapply(seq(1,74,4), function(i) {axis(1, at=i, labels = parse(text=lons5.lab[,1][i]),las=2)}))
axis(3,at=c(11,33,60),labels = c("ATLANTIC SECTOR","INDIAN SECTOR","PACIFIC SECTOR"),las=1, tick=F,line = -1)
legend("topright",legend=c("Antarctic krill","Ice krill"),pt.bg=col_used[1:2],pch=22,box.col="transparent",bg="transparent")
par(opar)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# divide by ocean sector -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


lons  <- data.frame(x1 = c(-70,30,150),
                    x2 = c(30,150,-70))
lats1 <- -90
lats2 <- 88

### split into 5 degree  longitude chunks
poly_section <- lapply(X = 1:nrow(lons),
                       FUN = function(x) st_polygon(list(rbind(c(lons[x, 1], (lats1)),
                                                               c(lons[x, 2], lats1),
                                                               c(lons[x, 2], lats2),
                                                               c(lons[x, 1], lats2),
                                                               c(lons[x, 1], lats1)))))
poly_section <- st_sfc(poly_section, crs = 4326)
poly_section <- st_transform(poly_section, south_pole_equal_area.proj)
# poly_section <- st_intersection(poly_section, circumpolar)
poly.out     <- data.frame(section = c("ATL","IND","PAC"),
                           e.sup.area = 0,
                           e.cry.area = 0,
                           over.area  = 0,
                           total.area = 0)

for(i in 1:3){
  poly.result <- mask(e_overlap,as_Spatial(st_zm(poly_section[i])))
  poly.out$total.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[2:4]))*100/1000000
  poly.out$over.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[4]))*100/1000000
  poly.out$e.sup.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[2]))*100/1000000
  poly.out$e.cry.area[i] <- sum(as.numeric(table(factor(values(poly.result),levels=c(0:3)))[3]))*100/1000000
}
poly.out[,2] <- poly.out[,2]+poly.out[,4]
poly.out[,3] <- poly.out[,3]+poly.out[,4]
poly.out[,2] <- poly.out[,2]/sum(poly.out[,2])
poly.out[,3] <- poly.out[,3]/sum(poly.out[,3])
poly.out[,4] <- poly.out[,4]/sum(poly.out[,4])


png("figures/krill proportion of habitat by ocean sectors.png",units="cm",
    res=800, width=8, height=18)

opar <- par(mfrow=c(1,1), mar=c(3,5,2,1))
cols_bar <- c(brewer.pal(6,"Accent"))[3:5]
barplot((as.matrix(poly.out[,2:4])),las=1,xaxt="n",
        ylab="Proportion of habitat", col=cols_bar, border=1)
axis(4, at =c(0.15, 0.45, 0.8), labels = c("ATLANTIC SECTOR", "INDIAN SECTOR", "PACIFIC SECTOR"),tick=F, line = -1)
axis(1, at =c(0.7, 1.9, 3.1), labels = c("Antarctic","Ice","area of"),tick=F, line = -0.5)
axis(1, at =c(0.7, 1.9, 3.1), labels = c("krill","krill","overlap"),tick=F, line = 0.5)
par(opar)

dev.off()







# 4 panel figure

png("figures/E.superba and E.crystallorophias 4 panel figure v2.png", res = 800, width=40, height = 40, units="cm")
opar <- par(mfrow=c(2,2), mar=c(0,0,0,0))

plot(st_geometry(model.domain),lty=2)  
plot(mask(es_mean, model.domain),zlim=c(0,1000), add=T, legend =F)
# plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)

plot(st_geometry(st_intersection(PF_fronts,circumpolar)), add=T, border=1,lty=2, lwd=0.8)
text(-1600, -3450, "PF",cex=0.7, col=1)
plot(st_geometry(SACCF_fronts), add=T, border=1,lty=1, lwd=0.8)
text(-300, -3050, "SACCF",cex=0.7, col=1)
contour(rast(env$NSIDC_ice_duration), levels = c(10), lty = c(3), add=T, lwd = 0.8, drawlabels = F)
text(2000, -2800, "ICE",cex=0.7, col=1)

plot(st_geometry(circumpolar),add=T)
plot(es_mean,legend.only = T, zlim=c(0,1000), 
     axis.args=list(at=c(0,1000), labels= c(0, 1),
                    cex.axis=1,tick=F,line=-0.2),
     legend.args=list(text="probability", side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
mtext('a)', side=3, line= -2.5, at=-4350, cex=2.5)


plot(st_geometry(circumpolar))
plot(mask(ec_mean, model.domain),zlim=c(0,1000), add=T, legend =F)
# plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)

plot(st_geometry(st_intersection(PF_fronts,circumpolar)), add=T, border=1,lty=2, lwd=0.8)
text(-1600, -3450, "PF",cex=0.7, col=1)
plot(st_geometry(SACCF_fronts), add=T, border=1,lty=1, lwd=0.8)
text(-300, -3050, "SACCF",cex=0.7, col=1)
contour(rast(env$NSIDC_ice_duration), levels = c(10), lty = c(3), add=T, lwd = 0.8, drawlabels = F)
text(2000, -2800, "ICE",cex=0.7, col=1)

plot(st_geometry(circumpolar),add=T)
mtext('b)', side=3, line= -2.5, at=-4350, cex=2.5)


cols_used <- c(grey(0.9),brewer.pal(5,"Set3")[c(4,5,1)])
plot(st_geometry(circumpolar))
plot(e_data,  col = grey(0.95), axes=F,legend=F,add=T)
plot(es_pol_hatch$geometry, col = cols_used[2], add=T, lwd =0.5)
plot(ec_pol_hatch$geometry, col = cols_used[3], add=T, lwd =0.5)
plot(es_pol, border = cols_used[2], add=T, lwd =0.5)
plot(ec_pol, border = cols_used[3], add=T, lwd =0.5)
plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(circumpolar, add=T)
plot(st_geometry(st_point_on_surface(ccamlr_domains_crop)), add=T, col=grey(0.2), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1.5)
legend("bottomright",legend=c("No data","Antarctic krill","Ice krill"),pch=22,
       pt.bg = c("white",cols_used[2:3]), col = 1, box.col = "transparent", cex = 1.3)
mtext('c)', side=3, line= -2.5, at=-4350, cex=2.5)


plot(st_geometry(circumpolar))
plot(es_dat,  col = cols_used[2], pch=3, axes=F,legend=F,add=T, cex=0.5)
plot(ec_dat,  col = cols_used[3], pch=4, axes=F,legend=F,add=T, cex=0.5)
plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(circumpolar, add=T)
plot(st_geometry(st_point_on_surface(ccamlr_domains_crop)), add=T, col=grey(0.2), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1.5)
mtext('d)', side=3, line= -2.5, at=-4350, cex=2.5)

par(opar)
dev.off()





# 2 panel figure

png("figures/E.superba and E.crystallorophias sd.png", res = 800, width=40, height = 20, units="cm")
opar <- par(mfrow=c(1,2), mar=c(0,0,0,0))

col_sd <- brewer.pal(9, "Purples")
plot(st_geometry(model.domain),lty=2)  
plot(mask(es_sd, model.domain),zlim=c(0,300), add=T, legend =F, col=col_sd)
# plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)

plot(st_geometry(st_intersection(PF_fronts,circumpolar)), add=T, border=1,lty=2, lwd=0.8)
text(-1600, -3450, "PF",cex=0.7, col=1)
plot(st_geometry(SACCF_fronts), add=T, border=1,lty=1, lwd=0.8)
text(-300, -3050, "SACCF",cex=0.7, col=1)
contour(rast(env$NSIDC_ice_duration), levels = c(10), lty = c(3), add=T, lwd = 0.8, drawlabels = F)
text(2000, -2800, "ICE",cex=0.7, col=1)

plot(st_geometry(circumpolar),add=T)
plot(es_sd/1000, legend.only = T, zlim=c(0,0.3), col=col_sd, 
     axis.args=list(cex.axis=1,tick=F,line=-0.2),
     legend.args=list(text="standard deviation", side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.85,0.87, 0.03,0.17))
mtext('a)', side=3, line= -2.5, at=-4350, cex=2.5)


plot(st_geometry(circumpolar))
plot(mask(ec_sd, model.domain),zlim=c(0,300), add=T, legend =F, col=col_sd)
# plot(st_geometry(ccamlr_domains_crop), add=T, border=grey(0.4),lty=3)
plot(st_geometry(ice_shelf),add=T,border=grey(0.7),col=grey(0.7),lwd=0.1)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)

plot(st_geometry(st_intersection(PF_fronts,circumpolar)), add=T, border=1,lty=2, lwd=0.8)
text(-1600, -3450, "PF",cex=0.7, col=1)
plot(st_geometry(SACCF_fronts), add=T, border=1,lty=1, lwd=0.8)
text(-300, -3050, "SACCF",cex=0.7, col=1)
contour(rast(env$NSIDC_ice_duration), levels = c(10), lty = c(3), add=T, lwd = 0.8, drawlabels = F)
text(2000, -2800, "ICE",cex=0.7, col=1)

plot(st_geometry(circumpolar),add=T)
# plot(st_geometry(st_point_on_surface(ccamlr_domains_crop)), add=T, col=grey(0.2), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1.5)
mtext('b)', side=3, line= -2.5, at=-4350, cex=2.5)

par(opar)
dev.off()
