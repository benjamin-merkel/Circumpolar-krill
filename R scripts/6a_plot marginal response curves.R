library(raster)
library(biomod2)
library(concaveman)
library(sf)

sf_use_s2(F)
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 

grat.50S           <- st_read("data/map data/ne_50m_graticules_5.shp")
grat.50S           <- st_transform(grat.50S[grat.50S$direction=="S" & grat.50S$degrees==50,][,1], south_pole_equal_area.proj)
grat.50S           <- rbind(grat.50S, grat.50S[1,])
grat.50S           <- concaveman(grat.50S) # works
grat.50S           <- st_buffer(grat.50S,dist = 0)
circumpolar        <- st_zm(grat.50S, drop = T,what = "ZM")



dat_Ec <- readRDS(paste0("data/E.crystallorophias_circumpolar_response_data_Dec-Mar.rds"))
dat_Es <- readRDS(paste0("data/E.superba_circumpolar_response_data_Dec-Mar.rds"))

base_dir       <- getwd()
gbm_output_dir <- paste0(base_dir, "/data/gbm")

setwd(gbm_output_dir)
temp.path     <- list.files(pattern="E.s.GBM Dec-Mar.models.out", recursive = TRUE)
biomod_gbm_Es <- get(load(temp.path))
temp.path     <- list.files(pattern="E.c.GBM Dec-Mar.models.out", recursive = TRUE)
biomod_gbm_Ec <- get(load(temp.path))

Biomodresponse_Es <- bm_PlotResponseCurves(bm.out = biomod_gbm_Es, fixed.var = "median", do.plot = F)
Biomodresponse_Ec <- bm_PlotResponseCurves(bm.out = biomod_gbm_Ec, fixed.var = "median", do.plot = F)


setwd(base_dir)




env.selected <- readRDS("data/Environmental covariates selected.rds")
chosen.pred <- names(env.selected)

env          <- readRDS("data/Environmental covariate stack.rds")
env.subset   <- subset(env, c("NSIDC_ice_retreat", "NSIDC_ice_conc", "NSIDC_spring_edge",
                              "WOA_sal_0", "WOA_sal_200", "OCCCI_bd","OCCCI_bm", "bath", 
                              "dis_1000","dis_ice.pol","Pellichero_ml_depth"))
names(env.subset) <- chosen.pred[-length(chosen.pred)]
env.subset   <- mask(env.subset, circumpolar)

env.sd <- env.mean <- vector(mode="list")
for(i in 1:length(chosen.pred[chosen.pred!="PCA"])) {
  env.mean[[i]] <- mean(values(raster::subset(env.subset, chosen.pred[chosen.pred!="PCA"][i])), na.rm=T)
  env.sd[[i]]   <- sd(values(raster::subset(env.subset, chosen.pred[chosen.pred!="PCA"][i]) - env.mean[[i]]), na.rm=T)
}
names(env.mean) <- names(env.sd) <- chosen.pred[chosen.pred!="PCA"]



setwd(base_dir)
png(paste0("figures/both krill GBM circumpolar marginal response curves_Dec-Mar.png"),
    units="cm",res=500,width=30, height=22)
par(mfrow = c(3, 4),mar=c(4, 0.2, 0.1, 0.2), oma = c(0, 4, 0, 0))
for(i in 1:length(chosen.pred)){
  xx_Es <- Biomodresponse_Es$tab[Biomodresponse_Es$tab$expl.name==chosen.pred[i],c("expl.val","pred.val")]
  xx_Es <- xx_Es[complete.cases(xx_Es),]
  response_Es <- data.frame(mean   = tapply(xx_Es$pred.val, xx_Es$expl.val, mean))
  response_Es$x <- as.numeric(as.character(rownames(response_Es)))
  
  xx_Ec <- Biomodresponse_Ec$tab[Biomodresponse_Ec$tab$expl.name==chosen.pred[i],c("expl.val","pred.val")]
  xx_Ec <- xx_Ec[complete.cases(xx_Ec),]
  response_Ec <- data.frame(mean   = tapply(xx_Ec$pred.val, xx_Ec$expl.val, mean))
  response_Ec$x <- as.numeric(as.character(rownames(response_Ec)))
  
  if(chosen.pred[i] != "PCA") x.Es <- (response_Es$x * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else x.Es <- response_Es$x
  if(chosen.pred[i] != "PCA") x.Ec <- (response_Ec$x * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else x.Ec <- response_Ec$x
  
  plot(x.Es, response_Es$mean,
       xlim=c(range(c(x.Es,x.Ec))),ylim=c(0,1),yaxt = "n", xaxt="n",
       main='',ylab="",xlab="",col="white")
  
  if(i %in% c(1,5,9)) axis(2, las = 1, cex.axis=0.8)
  if(!chosen.pred[i] %in% c("PCA","DIS","BM")) axis(1, las = 1, cex.axis=0.8)
  if(chosen.pred[i] == "DIS") axis(1, las=1, cex.axis=0.8, at = seq(-400,1200,400))
  if(chosen.pred[i] == "BM")  axis(1, las=1, cex.axis=0.8, at = log(c(0.2,0.5,1,2,5)),lab=c(0.2,0.5,1,2,5))
  if(chosen.pred[i]  == "PCA") {
    mtext("Dissolved oxygen", side = 1, line = 1, cex=0.8)
    mtext("Temperature", side = 1, line = 2, cex=0.8)
    axis(1, at = c(-3, 3), labels = c("+", "-"))
    axis(1, at = c(-3, 3), labels = c("-", "+"), line = 1, tick=F)
  }
  polygon(c(x.Ec, rev(x.Ec)), 
          c(tapply(xx_Ec$pred.val, xx_Ec$expl.val,min), rev(tapply(xx_Ec$pred.val, xx_Ec$expl.val,max))), 
          border="transparent", col=rgb(0,0,1,0.1))
  polygon(c(x.Es, rev(x.Es)), 
          c(tapply(xx_Es$pred.val, xx_Es$expl.val,min), rev(tapply(xx_Es$pred.val, xx_Es$expl.val,max))), 
          border="transparent", col=rgb(1,0,0,0.1))
  
  lines(x.Es, response_Es$mean, lwd=3, col=2)
  lines(x.Ec, response_Ec$mean, lwd=3, col=4)
  
  if(chosen.pred[i] != "PCA") r.Es <- (data.frame(dat_Es)[,chosen.pred[i]] * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else r.Es <- data.frame(dat_Es)[,chosen.pred[i]]
  if(chosen.pred[i] != "PCA") r.Ec <- (data.frame(dat_Ec)[,chosen.pred[i]] * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else r.Ec <- data.frame(dat_Ec)[,chosen.pred[i]]
  
  rug(r.Ec, side = 1, col = 4)
  rug(r.Es, side = 3, col = 2)
  
  
  if(chosen.pred[i] == "ICEMELT") mtext("Timing of sea ice retreat [days since 15 Feb]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "ICECONC") mtext("Sea ice concentration [%]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "MIZ")     mtext("Frequency of marginal ice zone [days]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "SSS")     mtext("Surface salinity [ppt]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "S200")    mtext("Salinity at 200m [ppt]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "BD")      mtext("Bloom duration [days]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "BM")      mtext(expression(Mean~Chl~a~during~the~bloom~'['~mg/m^3~']'), side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "BATH")    mtext("Bathymetry [m]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "DIS")     mtext("Distance to the shelf edge [km]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "ICEPOL")  mtext("Distance to coastal polynyas [km]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "MLD")     mtext("Mixed layer depth [m]", side = 1, line = 2.2, cex=0.8)
  
}
legend(0.8, legend = c("Antarctic krill","Ice krill"), lwd=c(4,4),lty=1,col=c(2,4),cex=0.9)
mtext("Probability", side = 2, outer =T, cex = 0.8, line = 2.5)
par(opar)
dev.off()





chosen.pred <- c("MLD", "PCA", "DIS", "SSS", "BM", "MIZ")

setwd(base_dir)
png(paste0("figures/both krill GBM circumpolar 6 most important marginal response curves_Dec-Mar.png"),
    units="cm",res=500,width=17, height=21)
par(mfrow = c(3, 2),mar=c(4, 0.2, 1.1, 0.2), oma = c(0, 4, 0, 0))
for(i in 1:length(chosen.pred)){
  xx_Es <- Biomodresponse_Es$tab[Biomodresponse_Es$tab$expl.name==chosen.pred[i],c("expl.val","pred.val")]
  xx_Es <- xx_Es[complete.cases(xx_Es),]
  response_Es <- data.frame(mean   = tapply(xx_Es$pred.val, xx_Es$expl.val, mean))
  response_Es$x <- as.numeric(as.character(rownames(response_Es)))
  
  xx_Ec <- Biomodresponse_Ec$tab[Biomodresponse_Ec$tab$expl.name==chosen.pred[i],c("expl.val","pred.val")]
  xx_Ec <- xx_Ec[complete.cases(xx_Ec),]
  response_Ec <- data.frame(mean   = tapply(xx_Ec$pred.val, xx_Ec$expl.val, mean))
  response_Ec$x <- as.numeric(as.character(rownames(response_Ec)))
  
  if(chosen.pred[i] != "PCA") x.Es <- (response_Es$x * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else x.Es <- response_Es$x
  if(chosen.pred[i] != "PCA") x.Ec <- (response_Ec$x * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else x.Ec <- response_Ec$x
  
  plot(x.Es, response_Es$mean,
       xlim=c(range(c(x.Es,x.Ec))),ylim=c(0,1),yaxt = "n", xaxt="n",
       main='',ylab="",xlab="",col="white")
  
  mtext(paste0(c("a","b","c","d","e","f")[i],") ", chosen.pred[i]), side = 3, at = min(c(x.Es,x.Ec)), adj=0, cex = 0.8)
  
  if(i %in% c(1,3,5)) axis(2, las = 1, cex.axis=0.8)
  if(!chosen.pred[i] %in% c("PCA","DIS","BM")) axis(1, las = 1, cex.axis=0.8)
  if(chosen.pred[i] == "DIS") axis(1, las=1, cex.axis=0.8, at = seq(-400,1200,400))
  if(chosen.pred[i] == "BM")  axis(1, las=1, cex.axis=0.8, at = log(c(0.2,0.5,1,2,5)),lab=c(0.2,0.5,1,2,5))
  if(chosen.pred[i] == "BM")  axis(1, las=1, cex.axis=0.8, at = log(c(seq(0.1,1,0.1),seq(2,10,1))),lab=NA)
  if(chosen.pred[i]  == "PCA") {
    mtext("Dissolved oxygen", side = 1, line = 1, cex=0.8)
    mtext("Temperature", side = 1, line = 2, cex=0.8)
    axis(1, at = c(-3, 3), labels = c("+", "-"))
    axis(1, at = c(-3, 3), labels = c("-", "+"), line = 1, tick=F)
  }
  polygon(c(x.Ec, rev(x.Ec)), 
          c(tapply(xx_Ec$pred.val, xx_Ec$expl.val,min), rev(tapply(xx_Ec$pred.val, xx_Ec$expl.val,max))), 
          border="transparent", col=rgb(0,0,1,0.1))
  polygon(c(x.Es, rev(x.Es)), 
          c(tapply(xx_Es$pred.val, xx_Es$expl.val,min), rev(tapply(xx_Es$pred.val, xx_Es$expl.val,max))), 
          border="transparent", col=rgb(1,0,0,0.1))
  
  lines(x.Es, response_Es$mean, lwd=3, col=2)
  lines(x.Ec, response_Ec$mean, lwd=3, col=4)
  
  if(chosen.pred[i] != "PCA") r.Es <- (data.frame(dat_Es)[,chosen.pred[i]] * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else r.Es <- data.frame(dat_Es)[,chosen.pred[i]]
  if(chosen.pred[i] != "PCA") r.Ec <- (data.frame(dat_Ec)[,chosen.pred[i]] * env.sd[[chosen.pred[i]]]) + env.mean[[chosen.pred[i]]] else r.Ec <- data.frame(dat_Ec)[,chosen.pred[i]]
  
  rug(r.Ec, side = 1, col = 4)
  rug(r.Es, side = 3, col = 2)
  
  if(chosen.pred[i] == "PCA") legend(1,0.8, legend = c("Antarctic krill","Ice krill"), lwd=c(4,4),lty=1,col=c(2,4),cex=0.9)

  
  if(chosen.pred[i] == "ICEMELT") mtext("Timing of sea ice retreat [days since 15 Feb]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "ICECONC") mtext("Sea ice concentration [%]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "MIZ")     mtext("Frequency of marginal ice zone [days]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "SSS")     mtext("Surface salinity [ppt]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "S200")    mtext("Salinity at 200m [ppt]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "BD")      mtext("Bloom duration [days]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "BM")      mtext(expression(Mean~Chl~a~during~the~bloom~'['~mg/m^3~']'), side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "BATH")    mtext("Bathymetry [m]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "DIS")     mtext("Distance to the shelf edge [km]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "ICEPOL")  mtext("Distance to coastal polynyas [km]", side = 1, line = 2.2, cex=0.8)
  if(chosen.pred[i] == "MLD")     mtext("Mixed layer depth [m]", side = 1, line = 2.2, cex=0.8)
  
}
mtext("Probability", side = 2, outer =T, cex = 0.8, line = 2.5)

dev.off()


