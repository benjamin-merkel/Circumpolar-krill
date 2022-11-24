#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages

library(RColorBrewer)
library(raster)
library(sf)
library(concaveman)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

fronts             <- st_read("data/map data/SouthernOceanFronts.shp")
fronts             <- st_transform(fronts, south_pole_equal_area.proj)
fronts             <- st_intersection(fronts, circumpolar)
PF_fronts          <- fronts[fronts$NAME == "Polar Front (PF)",]
SAF_fronts         <- fronts[fronts$NAME == "Subantarctic Front (SAF)",]
sACCf_fronts       <- fronts[fronts$NAME == "Southern Antarctic Circumpolar Current Front (sACCf)",]

env.selected <- readRDS("data/Environmental covariates selected.rds")
env          <- readRDS("data/Environmental covariate stack.rds")
names(env) 
names(env)   <- c("SST","T200","SSS","S200","MLD","MLT","MLS",
                  "BD","BI","BE","BT","BM","BS","BV","BA",
                  "ICEMELT","ICEADVANCE","ICEDUR","ICEEDGE","MIZ",
                  "ICECONC","ICEPERS","BATH","DIS","DISCOAST","SOX",
                  "OX200","SSI","SI200","SSN","N200","ICEPOL")
env2 <- env[[names(env.selected)]]


png("figures/Environmental covariates.png", units="cm",res=500,width=20, height=15)
par(mfrow=c(3,4), mar=c(0,0,0,0))
for(i in 1:12){
  plot(st_geometry(circumpolar))
  if(i<12) {x <- env2[[i]]} else {x <- env.selected[[i]]}
  xx <- mask(x, as_Spatial(circumpolar))
  xx <- mask(xx, as_Spatial(ice_shelf),inverse = T)
  plot(xx, col=brewer.pal(11,"Spectral"),add=T, legend=F)
  plot(st_geometry(ice_shelf),add=T,border=grey(0.3),col=grey(1),lwd=0.9)
  plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
  plot(st_geometry(circumpolar),add=T)
  title(names(env.selected[[i]]), line=-8.2, cex=0.8)
  plot(xx, legend.only=T,
       col=brewer.pal(11,"Spectral"),add=T,
       axis.args=list(cex.axis=0.5),
       legend.args=list(text='', side=4, font=2, line=1.5, cex=0.7),
       legend.width=0.2, legend.shrink=0.75,
       smallplot=c(0.85,0.87, 0.03,0.17))
}
dev.off()