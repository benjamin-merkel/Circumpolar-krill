# NSIDC-0051: Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS Passive Microwave Data, Version 1
# https://nsidc.org/data/NSIDC-0051/versions/1


library(RCurl)
library(stringr)
library(lubridate)
library(raster)
library(RColorBrewer)
library(sf)

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

files.in.folder <- list.files(path = directory, pattern = ".tif")
files.in.folder <- gsub(".tif","",files.in.folder)
years.used <- 1979:2019
doys.used  <- 304:396

rrr <- vector(mode = "list", length = length(years.used))
names(rrr) <- years.used
# advance <- retreat <- duration <- rrr

for (j in 1:length(years.used)) {
  
  iyear <- years.used[j]
  
  files.selected <- files[years2 == iyear & doy %in% doys.used]
  dates.selected <- dates[years2 == iyear & doy %in% doys.used]
  doy.selected   <- doy  [years2 == iyear & doy %in% doys.used]
  
  
  rr <- vector(mode = "list", length = length(files.selected))
  names(rr) <- dates.selected
  
  
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
    } else {
      
      r <- raster(paste0(directory,'\\',filename,".tif"))
      rr[[i]] <- r
      
    }
  }
  
  rr2 <- stack(rr)
  rr2[rr2<0.15] <- NA
  rr2[rr2>0.5]  <- NA
  rr2[!is.na(rr2)] <- 1
  rr3 <- sum(rr2,na.rm=T)
  
  rrr[[j]] <- rr3
  
  
}

# saveRDS(rrr,file="C:/Users/bme/OneDrive - Akvaplan-niva AS/Maud MPA/data/sea_ice_edge_duration_Nov-Jan.rds")
# rrr <- readRDS("C:/Users/bme/OneDrive - Akvaplan-niva AS/Maud MPA/data/sea_ice_edge_duration_Nov-Jan.rds")

rrr1 <- rrr

for (j in 1:length(years.used)) {
  
  iyear <- years.used[j]
  
  files.selected <- files[years2 == iyear & doy %in% doys.used]
  dates.selected <- dates[years2 == iyear & doy %in% doys.used]
  doy.selected   <- doy  [years2 == iyear & doy %in% doys.used]
  
  
  rrr1[[j]] <- rrr1[[j]]/length(files.selected)*92
  
}

rr5 <- mean(stack(rrr1))
rr5[is.na(r)] <- NA
writeRaster(rr5,"data/env covariates/sea ice edge duration climatology Nov-Jan.tif", format="GTiff", overwrite=TRUE)

