library(raster)
library(RColorBrewer)
library(beepr)
library(stringr)
library(oce)
library(sf)
library(concaveman)
library(beepr)
library(lubridate)
library(viridis)

#######################################################################################################
#EPSG:4326 lat lon unprojected -> same as Yoshi's inputcrs_csv
lonlat.proj <- CRS('+proj=longlat +datum=WGS84 +no_defs') 
# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 



### load yearly files in year bins

file_directory <- "D:/OC-CCI/8-day_Chla"

oc_cci <- list.files(file_directory, pattern = ".nc")
date  <- as.Date(strptime(str_split_fixed(oc_cci,"-",8)[,7], "%Y%m%d"))


cols_used <- colorRampPalette(c(brewer.pal(3,"YlGn")[1:3],
                                brewer.pal(3,"BuGn")[c(2,1)],
                                brewer.pal(3,"BuPu")[2:3],
                                brewer.pal(3,"RdPu")[c(2,1)],
                                brewer.pal(3,"Reds")[2:3]))


date.list <- c("0101", "0210", "0322", "0501", "0610", "0720", "0829", "1008", "1117")
dl = 1

be2_var<- bi2_var<- bd2_var<- bt2_var <- bs2_var <- bm2_var <- bv2_var <- ba2_var <- vector(mode="list",length=length(date.list))
be2_mean <- bi2_mean <- bd2_mean <- bt2_mean <- bs2_mean <- bm2_mean <- bv2_mean <- ba2_mean <- vector(mode="list",length=length(date.list))
years = 1998:2018


for(dl in 1:length(date.list)){
  
  years = 1998:2018
  # if(dl==9) years = 1998:2013
  # stats include : 
  # 1. initiation (bi), 
  # 2. end (be), 
  # 3. duration(bd), 
  # 4. timing of max (bt), 
  # 5. sum of Chla btw bi and be (bs)
  # 6. mean of Chla btw bi and be (bm)
  # 7. var of Chla btw bi and be (bm)
  # 8. amplitude (ba), 
  # 9. yearly minimum (year_min), 
  # 10.yearly maximum (year_max)
  
  be <- bi <- bd <- ba <- bt <- bs <- bm <- bv <- max <- min <- vector(mode="list",length=length(years))
  for(y in years){
    cat('\r',date.list[dl]," - ",y,"   ")
    # filename <- paste0('data/OC-CCI/yearly stats/Chla stat_',y,"_start_",date.list[dl],".grd")
    filename <- paste0('data/env covariates/OC-CCI/yearly stats/Chla stat double bin cutoff_',y,"_start_",date.list[dl],".grd")
    year_stat <- stack(filename)
    names(year_stat) <- c("bi","be","bd","bt","bs","bm","bv","ba","min","max")
    # names(year_stat) <- c("bi","be","bd")
    bt[[y-1997]] <- year_stat[["bt"]]
    bi[[y-1997]] <- year_stat[["bi"]]
    be[[y-1997]] <- year_stat[["be"]]
    bd[[y-1997]] <- year_stat[["bd"]]
    bs[[y-1997]] <- year_stat[["bs"]]
    bm[[y-1997]] <- year_stat[["bm"]]
    bv[[y-1997]] <- year_stat[["bv"]]
    ba[[y-1997]] <- year_stat[["ba"]]
    # max[[y-1997]] <- year_stat[["max"]]
    # min[[y-1997]] <- year_stat[["min"]]
    # 
    # date.select <- as.Date(c(as.Date(strptime(paste0(y,date.list[dl]),"%Y%m%d")):as.Date(strptime(paste0(y+1,date.list[dl]),"%Y%m%d")-1)),origin="1970-01-01")
    # date.select <- date[date %in% date.select]
    # doy.select <- as.numeric(strftime(date.select,"%j"))
    # 
    
    
  }
  # 
  bt <- stack(bt)
  bi <- stack(bi)
  be <- stack(be)
  bd <- stack(bd)
  bs <- stack(bs)
  bm <- stack(bm)
  bv <- stack(bv)
  ba <- stack(ba)
  # max <- stack(max)
  # min <- stack(min)
  #
  # 
  bi_mean <- calc(bi, mean, na.rm=T)
  be_mean <- calc(be, mean, na.rm=T)
  bd_mean <- calc(bd, mean, na.rm=T)
  bt_mean <- calc(bt, mean, na.rm=T)
  bs_mean <- calc(bs, mean, na.rm=T)
  bm_mean <- calc(bm, mean, na.rm=T)
  bv_mean <- calc(bv, mean, na.rm=T)
  ba_mean <- calc(ba, mean, na.rm=T)
  # max_mean <- calc(max, mean, na.rm=T)
  # min_mean <- calc(min, mean, na.rm=T)
  # # # 
  # bd[abs(bd) >= 200] <- sample(seq(-100000,100000,10),length(raster::as.array(bd)[raster::as.array(bd) >= 200]))
  # bi[bi == doy.select[1]] <- sample(seq(-100000,100000,10),length(raster::as.array(bi)[raster::as.array(bi) == doy.select[1]]))
  # be[be == rev(doy.select)[1]] <- sample(seq(-100000,100000,10),length(raster::as.array(be)[raster::as.array(be) == rev(doy.select)[[1]]))
  # 
  bi_var <- calc(bi, var, na.rm=T)
  be_var <- calc(be, var, na.rm=T)
  bd_var <- calc(bd, var, na.rm=T)
  bt_var  <- calc(bt, var, na.rm=T)
  #   # bt_var2 <- bt_var
  #   
  # bs_var <- calc(bs, var, na.rm=T)
  # bm_var <- calc(bm, var, na.rm=T)
  # ba_var <- calc(ba, var, na.rm=T)
  # max_var <- calc(max, var, na.rm=T)
  # min_var <- calc(min, var, na.rm=T)
  # # 
  # 
  #   # bt_var2[bd_mean > 200] <- 1000000
  # 
  # bs2_var[[dl]] <- bs_var
  # bm2_var[[dl]] <- bm_var
  bi2_var[[dl]] <- bi_var
  be2_var[[dl]] <- be_var
  bd2_var[[dl]] <- bd_var
  bt2_var[[dl]] <- bt_var
  # ba2_var[[dl]] <- ba_var
  # # 
  bi2_mean[[dl]] <- bi_mean
  be2_mean[[dl]] <- be_mean
  bd2_mean[[dl]] <- bd_mean
  bt2_mean[[dl]] <- bt_mean
  bs2_mean[[dl]] <- bs_mean
  bm2_mean[[dl]] <- bm_mean
  bv2_mean[[dl]] <- bv_mean
  ba2_mean[[dl]] <- ba_mean
  
}



doys <- c(as.numeric(strftime(strptime(paste0(2018,date.list),"%Y%m%d"),"%j")),400)

bt_mean <- bt2_mean[[1]]
for(k in 1:(length(doys)-1)){
  bt_mean[bt_mean >= doys[k] & bt_mean < doys[k+1]]  <- k-4
}
bt_mean <- round(focal(bt_mean, w=matrix(rep(1/9,9),ncol=3),na.rm=F),0)
bt_mean[bt_mean<1] <- bt_mean[bt_mean<1]+9



# coefficient of variation
mean_ba <- max_mean
sd_ba <- sqrt(max_var)
cov <- sqrt(max_var)/max_mean
# cov <- focal(cov, w=matrix(rep(1,9),ncol=3),median,na.rm=T)
# cov[cov>2] <- 2
# plot(cov, col=viridis(20))

bd_weight <- stack(bd2_mean)
bd_weight[bd_weight>=250] <- 14000
bd_weight[bd_weight>200 & bd_weight<250] <- 7000
bd_weight[bd_weight<400]<-0

# 
# plot(stack(bd2_mean),col=cols_used,zlim=c(0,366))
# plot(stack(bd_weight),col=cols_used)
# plot(min(stack(bd_weight)),col=cols_used)
# # plot(stack(bt2_var2),zlim=c(0,1000000))
# # plot(stack(bd2_var),zlim=c(0,70000))
# 
# plot(log(mean(stack(bm2_var),na.rm=T)))
# 
# 


btie_index<- apply(raster::as.array((stack(bt2_var)+stack(bi2_var)+stack(be2_var))/3 + bd_weight), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# btie_index<- apply(raster::as.array(stack(bt2_var)+stack(bi2_var)+stack(be2_var)), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# bti_index<- apply(raster::as.array(stack(bt2_var)+stack(bi2_var)), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# btied_index<- apply(raster::as.array(stack(bt2_var)+stack(bi2_var)+stack(be2_var)+stack(bd2_var)), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# bie_index <- apply(raster::as.array((stack(bi2_var)+stack(be2_var))/2 + bd_weight), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# bti_index <- apply(raster::as.array((stack(bt2_var)+stack(bi2_var))/2), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# bt_index  <- apply(raster::as.array(stack(bt2_var)), c(1,2), FUN = function(x) which(x == min(x, na.rm=T))[1])
# 
# bras <- focal(raster(bt_index), w=matrix(rep(1,9),ncol=3),median,na.rm=F)


plot(raster(bti_index),col=rainbow(9))
plot(raster(btie_index),col=rainbow(9))
plot(raster(btied_index),col=rainbow(9))

btie_ras <- background_raster
values(btie_ras) <- c(t(btie_index))

btie_ras_a4 <- mask(btie_ras,grat.50s)
btie_ras_a4 <- btie_ras_a4
btie_ras_a4[btie_ras_a4!=4] <- NA
btie_pol4 <- rasterToPolygons(btie_ras_a4,fun = function(x) x==4, dissolve = T)
btie_pol4 <- st_as_sf(btie_pol4)
btie_pol4 <- st_cast(btie_pol4, "POLYGON")
area    <- st_area(btie_pol4)
btie_pol4 <- btie_pol4[area == max(area),]
btie_pol4 <- nngeo::st_remove_holes(btie_pol4)

btie_ras_a <- mask(btie_ras,btie_pol4)
btie_ras_a[btie_ras_a!=4] <- 4
btie_ras_b <- mask(btie_ras,btie_pol4,inverse=T)

btie_ras_v2 <- sum(stack(btie_ras_a,btie_ras_b),na.rm = T)
btie_ras_v2[btie_ras_v2==0] <- NA

png("figures/bloom phenology yearly bin used.png",width=25,height=28,res=600,units="cm")
opar <- par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
plot((btie_ras),col=rainbow(9),legend=F,axes=F)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(dml.box),add=T,border=1,lwd=1.5)
plot(btie_ras,legend.only=T, col=rainbow(9),
     axis.args=list(at=c(1:9),
                    labels=strftime(strptime(date.list,"%m%d"),"%d-%b"), 
                    cex.axis=1),
     #legend.args=list(text=legend.label, side=4, font=2, line=2.5, cex=1),
     legend.width=0.2, legend.shrink=0.75,
     smallplot=c(0.81,0.83, 0.02,0.17)); par(mar = par("mar"))
par(new=T,fig=c(0,.23,0,.20),mar=c(0.5,2.3,0,0))
btie_table <- table(values(btie_ras))
names(btie_table) <- ""
btie_table <- btie_table/sum(btie_table)
barplot(btie_table,las=2,col=rainbow(9),line=-0.2,cex.axis=0.8)
par(opar)
dev.off()

# plot(raster(bti_index),col=rainbow(9))

weight_mat <- matrix(rep(1,25),ncol=5)
weight_mat[c(1,5),] <- 0.5
weight_mat[,c(1,5)] <- 0.5
weight_mat <-weight_mat/sum(weight_mat)

index_mat <- bie_index
index_mat <- btie_index
index_mat <- raster::as.matrix(bt_mean)
index_mat <- raster::as.matrix(btie_ras_v2)
index_mat <- btie_index


arr <- raster::as.array(abs(stack(bd2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
bd_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
bd_ras <- background_raster
values(bd_ras) <- c(t(bd_comp))

#bd_comp <- focal((bd_comp), w=weight_mat,na.rm=F)
plot(bd_ras,zlim=c(0,366),col=cols_used(52))


arr <- raster::as.array(abs(stack(bi2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
bi_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
bi_ras <- background_raster
values(bi_ras) <- c(t(bi_comp))
#bi_comp <- focal(bi_comp, w=weight_mat, na.rm=F)
bi_ras[bi_ras > 365 & !is.na(bi_ras)] <- bi_ras[bi_ras > 365 & !is.na(bi_ras)] -365
plot(bi_ras,col=brewer.pal(11,"Spectral"))


arr <- raster::as.array(abs(stack(be2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
be_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
be_ras <- background_raster
values(be_ras) <- c(t(be_comp))
#bi_comp <- focal(bi_comp, w=weight_mat, na.rm=F)
be_ras[be_ras > 365 & !is.na(be_ras)] <- be_ras[be_ras > 365 & !is.na(be_ras)] -365
plot(be_ras,col=rainbow(52))


arr <- raster::as.array(abs(stack(bt2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
bt_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
bt_ras <- background_raster
values(bt_ras) <- c(t(bt_comp))
#bt_comp <- focal(bt_comp, w=weight_mat, na.rm=F)
bt_ras[bt_ras > 365 & !is.na(bt_ras)] <- bt_ras[bt_ras > 365 & !is.na(bt_ras)] -365
plot(bt_ras,col=rainbow(52))



arr <- raster::as.array(abs(stack(bm2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
bm_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
bm_ras <- background_raster
values(bm_ras) <- c(t(bm_comp))
#bm_comp <- focal(log(raster(bm_comp)), w=weight_mat,na.rm=F)
plot(log(bm_ras),col=cols_used(200))


arr <- raster::as.array(abs(stack(bs2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
bs_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
bs_ras <- background_raster
values(bs_ras) <- c(t(bs_comp))
#bs_comp <- focal(log(raster(bs_comp)), w=weight_mat,na.rm=F)
plot(log(bs_ras),col=cols_used(200))


arr <- raster::as.array(abs(stack(bv2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
bv_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
bv_ras <- background_raster
values(bv_ras) <- c(t(bv_comp))
#bv_comp <- focal(log(raster(bv_comp)), w=weight_mat,na.rm=F)
plot(log(bv_ras),col=cols_used(200))


arr <- raster::as.array(abs(stack(ba2_mean)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
ba_comp <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
ba_ras <- background_raster
values(ba_ras) <- c(t(ba_comp))
#ba_comp <- focal(log(raster(ba_comp)), w=weight_mat,na.rm=F)
plot(log(ba_ras),col=cols_used(200))

arr <- raster::as.array(abs(stack(ba2_var)))#[,,c(1,4,7)]
i1 <- dim(arr)[1]
j1 <- dim(arr)[2]
ba_comp_var <- matrix(arr[cbind(rep(seq_len(i1),  j1),rep(seq_len(j1), each = i1), c(index_mat))], ncol=i1)
ba_ras_var <- background_raster
values(ba_ras_var) <- c(t(ba_comp_var))
plot(log(ba_ras_var),col=cols_used(200))


ba_cov <- sqrt((ba_ras_var))/(ba_ras)
cov2 <- sqrt((bv_ras))/(bm_ras)



writeRaster(bd_ras, "data/env covariates/Chla bd mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(bi_ras, "data/env covariates/Chla bi mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(be_ras, "data/env covariates/Chla be mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(bt_ras, "data/env covariates/Chla bt mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(bm_ras, "data/env covariates/Chla bm mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(bv_ras, "data/env covariates/Chla bv mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(bs_ras, "data/env covariates/Chla bs mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
writeRaster(ba_ras, "data/env covariates/Chla ba mean climatology 98_18.tif", format="GTiff", overwrite=TRUE)
