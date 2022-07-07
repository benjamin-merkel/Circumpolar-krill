library(ncdf4)
library(raster)

# WOA 2018 Statistical mean 
# Winter (Dec-Mar) -> 13
# https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18documentation.pdf 

# Objectively analyzed climatology - an
# Statistical mean - mn
# Number of observations - dd
# Seasonal or climatology minus annual climatology - ma
# Standard deviation from statistical mean - sd 
# Standard error of the statistical mean - se 
# Statistical mean minus objectively analyzed climatology - oa
# Number of mean values within radius of influence - gp

#EPSG:4326 lat lon unprojected
lonlat.proj <- CRS('+proj=longlat +datum=WGS84 +no_defs') 


### dissolved oxygen
ox <- raster('data/env covariates/WOA 18/woa18_all_o12_01.nc')
projection(ox) <- lonlat.proj
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_o12_01.nc')
depth <- ncvar_get(nc4,"depth")
tgt12 <- ncvar_get(nc4,varid = 'o_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_o01_01.nc')
tgt01 <- ncvar_get(nc4,varid = 'o_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_o02_01.nc')
tgt02 <- ncvar_get(nc4,varid = 'o_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_o03_01.nc')
tgt03 <- ncvar_get(nc4,varid = 'o_an')
nc_close(nc4)

X <- list(tgt12[,,depth==0], tgt01[,,depth==0], tgt02[,,depth==0], tgt03[,,depth==0])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_ox0 <- ox

X <- list(tgt12[,,depth==200], tgt01[,,depth==200], tgt02[,,depth==200], tgt03[,,depth==200])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_ox200 <- ox

### silicate
ox <- raster('data/env covariates/WOA 18/woa18_all_i12_01.nc')
projection(ox) <- lonlat.proj
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_i12_01.nc')
depth <- ncvar_get(nc4,"depth")
tgt12 <- ncvar_get(nc4,varid = 'i_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_i01_01.nc')
tgt01 <- ncvar_get(nc4,varid = 'i_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_i02_01.nc')
tgt02 <- ncvar_get(nc4,varid = 'i_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_i03_01.nc')
tgt03 <- ncvar_get(nc4,varid = 'i_an')
nc_close(nc4)

X <- list(tgt12[,,depth==0], tgt01[,,depth==0], tgt02[,,depth==0], tgt03[,,depth==0])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_si0 <- ox

X <- list(tgt12[,,depth==200], tgt01[,,depth==200], tgt02[,,depth==200], tgt03[,,depth==200])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_si200 <- ox

### nitrate
ox <- raster('data/env covariates/WOA 18/woa18_all_n12_01.nc')
projection(ox) <- lonlat.proj
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_n12_01.nc')
depth <- ncvar_get(nc4,"depth")
tgt12 <- ncvar_get(nc4,varid = 'n_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_n01_01.nc')
tgt01 <- ncvar_get(nc4,varid = 'n_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_n02_01.nc')
tgt02 <- ncvar_get(nc4,varid = 'n_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_all_n03_01.nc')
tgt03 <- ncvar_get(nc4,varid = 'n_an')
nc_close(nc4)

X <- list(tgt12[,,depth==0], tgt01[,,depth==0], tgt02[,,depth==0], tgt03[,,depth==0])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_ni0 <- ox

X <- list(tgt12[,,depth==200], tgt01[,,depth==200], tgt02[,,depth==200], tgt03[,,depth==200])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_ni200 <- ox


### temperature
ox <- raster('data/env covariates/WOA 18/woa18_decav81B0_t01_04.nc')
projection(ox) <- lonlat.proj
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_t12_04.nc')
depth <- ncvar_get(nc4,"depth")
tgt12 <- ncvar_get(nc4,varid = 't_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_t01_04.nc')
tgt01 <- ncvar_get(nc4,varid = 't_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_t02_04.nc')
tgt02 <- ncvar_get(nc4,varid = 't_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_t03_04.nc')
tgt03 <- ncvar_get(nc4,varid = 't_an')
nc_close(nc4)

X <- list(tgt12[,,depth==0], tgt01[,,depth==0], tgt02[,,depth==0], tgt03[,,depth==0])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_temp0 <- ox

X <- list(tgt12[,,depth==200], tgt01[,,depth==200], tgt02[,,depth==200], tgt03[,,depth==200])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_temp200 <- ox

### salinity
ox <- raster('data/env covariates/WOA 18/woa18_decav81B0_s01_04.nc')
projection(ox) <- lonlat.proj
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_s12_04.nc')
depth <- ncvar_get(nc4,"depth")
tgt12 <- ncvar_get(nc4,varid = 's_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_s01_04.nc')
tgt01 <- ncvar_get(nc4,varid = 's_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_s02_04.nc')
tgt02 <- ncvar_get(nc4,varid = 's_an')
nc_close(nc4)
nc4 <- nc_open('data/env covariates/WOA 18/woa18_decav81B0_s03_04.nc')
tgt03 <- ncvar_get(nc4,varid = 's_an')
nc_close(nc4)

X <- list(tgt12[,,depth==0], tgt01[,,depth==0], tgt02[,,depth==0], tgt03[,,depth==0])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_sal0 <- ox

X <- list(tgt12[,,depth==200], tgt01[,,depth==200], tgt02[,,depth==200], tgt03[,,depth==200])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
mean_tgt <- apply(Y, c(1, 2), mean, na.rm = TRUE)
values(ox) <-apply(t(mean_tgt),2,rev) 
WOA_sal200 <- ox


writeRaster(WOA_ox0,    "data/env covariates/WOA 18 dissloved oxygen 0m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_ox200,  "data/env covariates/WOA 18 dissloved oxygen 200m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_si0,    "data/env covariates/WOA 18 silicate 0m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_si200,  "data/env covariates/WOA 18 silicate 200m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_ni0,    "data/env covariates/WOA 18 nitrate 0m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_ni200,  "data/env covariates/WOA 18 nitrate 200m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_temp0,  "data/env covariates/WOA 18 temperature 0m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_temp200,"data/env covariates/WOA 18 temperature 200m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_sal0,   "data/env covariates/WOA 18 salinity 0m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
writeRaster(WOA_sal200, "data/env covariates/WOA 18 salinity 200m Dec-Mar climatology.tif", format="GTiff", overwrite=TRUE)
