# NSIDC-0051: Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS Passive Microwave Data, Version 1
# https://nsidc.org/data/NSIDC-0051/versions/1


library(RCurl)
library(stringr)
library(lubridate)
library(raster)
library(RColorBrewer)


h <- curl::new_handle()
curl::handle_setopt(
  handle = h,
  httpauth = 1,
  userpwd = "BM_earthdata:3arthData"
)

r.mask <- raster(ncol = 316,nrow = 332, xmn = -3950, xmx = 3950, ymn = -3950, ymx = 4350)
inputcrs_ice <- CRS('+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=km +no_defs')
projection(r.mask) <- inputcrs_ice
directory <- "E:\\environmental data\\sea ice\\sea ice 25km 79-20 south"


files <- read.table("E:/environmental data/sea ice/4363554859-download.txt")[,1]
files <- files[nchar(files) == 90]

files.in.folder <- list.files(path = directory, pattern = ".tif")
files.in.folder <- gsub(".tif","",files.in.folder)

files <- files[gsub(".bin","",str_split_fixed(files, '/',8)[,8]) %in% files.in.folder]


dates <- str_split_fixed(files,"/",10)[,7]
dates <- as.Date(strptime(dates,"%Y.%m.%d"))

files <- files[order(dates)]

dates <- str_split_fixed(files,"/",10)[,7]
dates <- as.Date(strptime(dates,"%Y.%m.%d"))

months<- as.numeric(strftime(dates, "%m"))
years <- as.numeric(strftime(dates, "%Y"))
doy   <- as.numeric(strftime(dates,"%j"))
years2 <- years
years2[doy < 46] <- years2[doy < 46]-1 # 15 Feb
for(i in unique(years)) {
  if(leap_year(i-1)) add.doy = 366 else add.doy = 365
  doy[years==i & doy<46] <- doy[years==i & doy<46] + add.doy
}



years.used <- 1979:2019

rrr <- vector(mode = "list", length = length(years.used))
names(rrr) <- years.used
advance <- retreat <- duration <- rrr

for (j in 1:length(years.used)) {
  
  iyear <- years.used[j]
  
  files.selected <- files[years2 == iyear]
  dates.selected <- dates[years2 == iyear]
  doy.selected   <- doy[years2 == iyear]
  
  rra <- rrb <- vector(mode = "list", length = length(files.selected))
  names(rra) <- names(rrb) <- dates.selected
  
  
  for(i in 1:length(files.selected)){
    filename <-  tail(strsplit(files.selected[i], '/')[[1]], n = 1) # Keep original filename
    filename <- gsub(".bin","",filename)
    
    cat("\r",iyear,' - ',doy.selected[i],'  ')
    
    if(!filename %in% files.in.folder){
      resp <- curl::curl_fetch_memory(files.selected[i], handle = h)
      
      xx <- as.numeric(resp$content[301:length(resp$content)])
      
      r <- setValues(r.mask, as.numeric(resp$content[301:length(resp$content)]))
      r[r==251]<-250 # circular data hole around north pole, assumed to be 100%
      r[r==253]<-NA # coastline
      r[r==254]<-NA # land mask
      r[r==255]<-NA # missing data
      r <- r/250
      
      rr[[i]] <- r
      
      #writeRaster(r, paste0(directory,'\\',filename,".tif"),format="GTiff", overwrite=TRUE,options=c("COMPRESS=NONE", "TFW=NO"))
    } else {
      
      r <- raster(paste0(directory,'\\',filename,".tif"))
      r[r<0.15] <- 0
      # r[r<0.15] <- NA
      
      r2 <- r
      r2[!r2 %in% 0 & !is.na(r2)] <- 1
      rra[[i]] <- r2
      
      r2 <- r
      r2[r2>0.5]  <- 0
      r2[!r2 %in% 0 & !is.na(r2)] <- 1
      rrb[[i]] <- r2
      
    }
  }
  
  if(leap_year(iyear)) add.doy = 366 else add.doy = 365
  
  rr2 <- stack(rra)
  names(rr2) <- doy.selected
  rr3 <- sum(rr2,na.rm=T)
  
  duration[[j]] <- rr3/length(files.selected)*add.doy
  duration[[j]][duration[[j]]<3] <- NA
  
  rr3 <- calc(rr2, function(x) ifelse(x==1,as.numeric(gsub("X","",names(x))),x))
  retreat[[j]]  <- max(rr3,na.rm=T)
  advance[[j]]  <- min(rr3,na.rm=T)
  retreat[[j]] [is.na(duration[[j]])] <- NA
  advance[[j]] [is.na(duration[[j]])] <- NA
  
  rr2 <- stack(rrb)
  rr3 <- sum(rr2,na.rm=T)
  
  rrr[[j]] <- rr3/length(files.selected)*add.doy
  
  
}

#rrb <- rr[order(doy.selected)]
saveRDS(rrr,     file="data/env covariates/sea_ice_edge_duration.rds")
saveRDS(duration,file="data/env covariates/sea_ice_duration.rds")
saveRDS(advance, file="data/env covariates/sea_ice_advance.rds")
saveRDS(retreat, file="data/env covariates/sea_ice_retreat.rds")

plot(mean(stack(retreat),na.rm=T),col=brewer.pal(11,"Spectral"))
plot(mean(stack(advance),na.rm=T),col=brewer.pal(11,"Spectral"))
plot(mean(stack(duration),na.rm=T),col=brewer.pal(11,"Spectral"))
plot(mean(stack(rrr),na.rm=T),col=brewer.pal(11,"Spectral"))

ox <- mean(stack(retreat),na.rm=T)
writeRaster(ox,"data/env covariates/sea ice retreat climatology.tif", format="GTiff", overwrite=TRUE)
ox <- mean(stack(advance),na.rm=T)
writeRaster(ox,"data/env covariates/sea ice advance climatology.tif", format="GTiff", overwrite=TRUE)
ox <- mean(stack(duration),na.rm=T)
writeRaster(ox,"data/env covariates/sea ice duration climatology.tif", format="GTiff", overwrite=TRUE)

rr5 <- mean(stack(rrr))
writeRaster(rr5,"data/env covariates/sea ice edge duration climatology.tif", format="GTiff", overwrite=TRUE)



r <- setValues(r.mask, as.numeric(resp$content[301:length(resp$content)]))
r[r<251]<-NA # circular data hole around north pole, assumed to be 100%
r[!is.na(r)]<-1
plot(r)

writeRaster(r,"data/env covariates/sea ice NA.tif", format="GTiff", overwrite=TRUE)
