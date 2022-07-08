library(raster)
library(biomod2)

env.selected <- readRDS("data/Environmental covariates selected.rds")

dat_Ec <- readRDS(paste0("data/E.crystallorophias_circumpolar_response_data_Dec-Mar.rds"))
dat_Es <- readRDS(paste0("data/E.superba_circumpolar_response_data_Dec-Mar.rds"))

base_dir       <- getwd()
gbm_output_dir <- paste0(base_dir, "/data/gbm")

setwd(gbm_output_dir)
temp.path <- list.files(pattern="E.s.GBM Dec-Mar.models.out", recursive = TRUE)
biomod_gbm_Es <- get(load(temp.path))
temp.path <- list.files(pattern="E.c.GBM Dec-Mar.models.out", recursive = TRUE)
biomod_gbm_Ec <- get(load(temp.path))

# gbm_Eval<- get_evaluations(biomod_gbm_Es)


Biomodresponse_Es <- response.plot2(
  models = BIOMOD_LoadModels(biomod_gbm_Es),
  Data = get_formal_data(biomod_gbm_Es, 'expl.var'),
  show.variables = get_formal_data(biomod_gbm_Es,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'mean', plot=F, legend = F,
  data_species = get_formal_data(biomod_gbm_Es, 'resp.var')
)


Biomodresponse_Ec <- response.plot2(
  models = BIOMOD_LoadModels(biomod_gbm_Ec),
  Data = get_formal_data(biomod_gbm_Ec, 'expl.var'),
  show.variables = get_formal_data(biomod_gbm_Ec,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'mean', plot=F, legend = F,
  data_species = get_formal_data(biomod_gbm_Ec, 'resp.var')
)



ncols <- ceiling((nlayers(env.selected)+1)/4)
setwd(base_dir)
png(paste0("figures/both krill GBM circumpolar marginal response curves_Dec-Mar.png"),
    units="cm",res=500,width=4*ncols, height=15)
par(mfrow = c(4,ncols),mar=c(2,2,2,0.1))
for(i in unique(Biomodresponse_Es$expl.name)){
  xx_Es <- Biomodresponse_Es[Biomodresponse_Es$expl.name==i,]
  xx_Es <- xx_Es[complete.cases(xx_Es),]
  xx_Es <- xx_Es[!grepl("Full", xx_Es$pred.name),]
  response_Es <- data.frame(mean   = tapply(xx_Es$pred.val, xx_Es$expl.val, mean))
  response_Es$x <- as.numeric(as.character(rownames(response_Es)))
  
  xx_Ec <- Biomodresponse_Ec[Biomodresponse_Ec$expl.name==i,]
  xx_Ec <- xx_Ec[complete.cases(xx_Ec),]
  xx_Ec <- xx_Ec[!grepl("Full", xx_Ec$pred.name),]
  response_Ec <- data.frame(mean   = tapply(xx_Ec$pred.val, xx_Ec$expl.val, mean))
  response_Ec$x <- as.numeric(as.character(rownames(response_Ec)))
  
  
  plot(response_Es$x, response_Es$mean,type="l",lwd=2,ylim=c(0,1),main=i,las=1,ylab="",xlab="",cex.axis=0.8,col="white")
  for(j in unique(Biomodresponse_Es$pred.name)){
    xx2 <- xx_Es[xx_Es$pred.name==j,]
    lines(xx2$expl.val, xx2$pred.val, lwd=1, col=rgb(1,0,0,0.2))
  }
  for(j in unique(Biomodresponse_Ec$pred.name)){
    xx2 <- xx_Ec[xx_Ec$pred.name==j,]
    lines(xx2$expl.val, xx2$pred.val, lwd=1, col=rgb(0,0,1,0.2))
  }
  lines(response_Es$x, response_Es$mean, lwd=3, col=2)
  lines(response_Ec$x, response_Ec$mean, lwd=3, col=4)
  
  rug(data.frame(dat_Ec)[,i], side = 1, col = 4)
  rug(data.frame(dat_Es)[,i], side = 3, col = 2)
}
plot(1,1,col="white",ann=F,axes=F)
legend("bottomright", legend = c("Antarctic krill","Ice krill","ensemble mean","single model"),
       lwd=c(4,4,4,1),lty=1,col=c(2,4,1,1),cex=0.9)

par(opar)
dev.off()








chosen.pred <- c("MLD","PCA","DIS",'SSI')

env          <- readRDS("data/Environmental covariate stack.rds")
env.subset   <- subset(env, c("Pellichero_ml_depth","dis_1000",'WOA_si_0'))
names(env.subset) <- c("MLD","DIS","SSI")

env.sd <- env.mean <- vector(mode="list")
for(i in 1:length(chosen.pred[chosen.pred!="PCA"])) {
  env.mean[[i]] <- mean(values(raster::subset(env.subset, chosen.pred[chosen.pred!="PCA"][i])), na.rm=T)
  env.sd[[i]]   <- sd(values(raster::subset(env.subset, chosen.pred[chosen.pred!="PCA"][i])), na.rm=T)
}
names(env.mean) <- names(env.sd) <- chosen.pred[chosen.pred!="PCA"]



setwd(base_dir)
png(paste0("figures/both krill GBM circumpolar 4 most important marginal response curves_Dec-Mar.png"),
    units="cm",res=500,width=17, height=16)
par(mfrow = c(2, 2),mar=c(4, 0.2, 0.1, 0.2), oma = c(0, 2, 0, 0))
for(i in chosen.pred){
  xx_Es <- Biomodresponse_Es[Biomodresponse_Es$expl.name==i,]
  xx_Es <- xx_Es[complete.cases(xx_Es),]
  xx_Es <- xx_Es[!grepl("Full", xx_Es$pred.name),]
  response_Es <- data.frame(mean= tapply(xx_Es$pred.val, xx_Es$expl.val, mean),
                            min = tapply(xx_Es$pred.val, xx_Es$expl.val, min),
                            max = tapply(xx_Es$pred.val, xx_Es$expl.val, max))
  response_Es$x <- as.numeric(as.character(rownames(response_Es)))
  
  xx_Ec <- Biomodresponse_Ec[Biomodresponse_Ec$expl.name==i,]
  xx_Ec <- xx_Ec[complete.cases(xx_Ec),]
  xx_Ec <- xx_Ec[!grepl("Full", xx_Ec$pred.name),]
  response_Ec <- data.frame(mean= tapply(xx_Ec$pred.val, xx_Ec$expl.val, mean),
                            min = tapply(xx_Ec$pred.val, xx_Ec$expl.val, min),
                            max = tapply(xx_Ec$pred.val, xx_Ec$expl.val, max))
  response_Ec$x <- as.numeric(as.character(rownames(response_Ec)))
  
  if(i!="PCA") x.Es <- (response_Es$x*env.sd[[i]])+env.mean[[i]] else x.Es <- response_Es$x
  if(i!="PCA") x.Ec <- (response_Ec$x*env.sd[[i]])+env.mean[[i]] else x.Ec <- response_Ec$x
  
  plot(x.Es, response_Es$mean,type="l",lwd=2,ylim=c(0,1),
       xlim=c(min(c(x.Es, x.Ec)), max(c(x.Es, x.Ec))), yaxt="n", xaxt="n",
       main="",las=1, ylab="",xlab="",cex.axis=0.8, col="white")
  if(i %in% c("MLD", "DIS")) axis(2, las=1, cex.axis=0.8)
  if(!i %in% c("MLD", "PCA")) axis(1, las=1, cex.axis=0.8)
  polygon(c(x.Ec, rev(x.Ec)), c(response_Ec$min, rev(response_Ec$max))-0.1, 
          border="transparent", col=rgb(0,0,1,0.1))
  polygon(c(x.Es, rev(x.Es)), c(response_Es$min, rev(response_Es$max))+0.1, 
          border="transparent", col=rgb(1,0,0,0.1))
  lines(x.Ec, response_Ec$mean-0.1, lwd=2, col=4)
  lines(x.Es, response_Es$mean+0.1, lwd=2, col=2)
  
  # if(i == 'pca'){
  #   arrows(-4, 0.1, 0, 0.1, length = 0.1)
  #   text(-2, 0.115, "dissolved oxygen")
  #   arrows(0, 0.05, -4, 0.05, length = 0.1)
  #   text(-2, 0.065, "temperature")
  # }
  if(i == "PCA") {
    mtext("Dissolved oxygen", side = 1, line = 1, cex=0.8)
    mtext("Temperature", side = 1, line = 2, cex=0.8)
    axis(1, at = c(-3, 2.6), labels = c("+", "-"))
    axis(1, at = c(-3, 2.6), labels = c("-", "+"), line = 1, tick=F)
  }
  if(i == "MLD") axis(1, at = seq(40,300,20),las=1, cex.axis=0.8)
  
  if(i == 'SSI') {
    legend("right", legend = c("Antarctic krill","Ice krill"),
           lwd=c(4,4),lty=1,col=c(2,4),cex=0.9)
  }
  
  if(i == "MLD") mtext("Mixed layer depth [m]", side = 1, line = 2.2, cex=0.8)
  # if(i == "pca")         mtext("PCA", side = 1, line = 2.2, cex=0.8)
  if(i == "DIS")    mtext("Distance to the shelf edge [km]", side = 1, line = 2.2, cex=0.8)
  if(i == "SSI")    mtext(expression(paste("Surface silicate [", mu, "mol/kg]")), side = 1, line = 2.1, cex=0.8)
  if(i == "SSS")   mtext("Surface salinity [psu]", side = 1, line = 2.2, cex=0.8)
  if(i == "S200") mtext("Salinity at 200m [psu]", side = 1, line = 2.2, cex=0.8)
  
  if(i!="PCA") x.Es <- (data.frame(dat_Es)[,i]*env.sd[[i]])+env.mean[[i]] else x.Es <- data.frame(dat_Es)[,i]
  if(i!="PCA") x.Ec <- (data.frame(dat_Ec)[,i]*env.sd[[i]])+env.mean[[i]] else x.Ec <- data.frame(dat_Ec)[,i]
  
  rug(x.Ec, side = 1, col = 4)
  rug(x.Es, side = 3, col = 2)
}

par(opar)
dev.off()


