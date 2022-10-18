#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages

pckg <- c('raster','fasterize','ncdf4','sf','sp','rgdal','rgeos','terra',
          'nngeo','stringr','readxl','classInt','smoothr','usdm',
          'ecospat',"e1071",'spatialEco', "ggplot2",
          'viridis','RColorBrewer',"scales","remote","angstroms","data.table","fmsb",
          'adehabitatHR','concaveman','corrplot','RStoolbox',"PresenceAbsence")

for(i in 1:length(pckg)) do.call("library", list(pckg[i]))

library(sf)
library(raster)
library(plotrix)

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

background_polygon <- st_difference((grat.50S),st_geometry(ice_shelf))
background_polygon <- st_difference(background_polygon,st_geometry(land))
background_raster  <- raster(as_Spatial(background_polygon),res=10)
values(background_raster) <- 1:ncell(background_raster)


model.domain <- circumpolar

ccamlr_mpa <- st_read("data/CCAMLR/bm_mpa_planningDomains.shp")
ccamlr_mpa <- st_transform(ccamlr_mpa, south_pole_equal_area.proj)
ccamlr_mpa_crop <- st_intersection(ccamlr_mpa, circumpolar)



ibsco <- raster('data/env covariates/IBCSO_v2_ice-surface.tif')
# gebco <- raster('data/gebco 30S and below.tif')
# gebco <- resample(gebco, raster(crs=lonlat.proj,resolution=0.01,xmn=-180,xmx=180,ymn=-90,ymx=-40))
bath  <- terra::project(rast(ibsco), rast(background_raster))
bath[bath>0] <- 0
bath <- raster(bath)


# set up graticules
lons <- c(-70,30,150)#seq(-180, 180, 60)
lats <- -50#c(-80, -60)
grat <- st_graticule(st_as_sf(circumpolar), lat = lats, lon = lons)
grat <- st_transform(grat, projection(bath))
g <- grat <- st_intersection(grat, circumpolar)




bb <- brewer.pal(9, 'Blues')
bb <- colorRampPalette(bb)(11)
blue <- colorRampPalette(c(bb[3:9], 1))(100)
blue <- colorRampPalette(c(bb[2:8]))(100)
island.col <- brewer.pal(9,"YlOrRd")[8] #brewer.pal(2,"Set3")[2]
sea.col <- brewer.pal(9,"Purples")[9] 




# setwd(maud_directory)
png("figures/study area v3.png", res = 800, width=20, height = 20, units="cm")
par(mar=rep(0,4))
plot(st_geometry(grat),lty=2)  
plot(mask(bath, circumpolar),  zlim=c(-1000,0), add=T, legend =F, col=rev(bb[2]))
plot(mask(bath, circumpolar),  zlim=c(-2000,-1000), add=T, legend =F, col=rev(bb[3]))
plot(mask(bath, circumpolar),  zlim=c(-3000,-2000), add=T, legend =F, col=rev(bb[4]))
plot(mask(bath, circumpolar),  zlim=c(-4000,-3000), add=T, legend =F, col=rev(bb[5]))
plot(mask(bath, circumpolar),  zlim=c(-5000,-4000), add=T, legend =F, col=rev(bb[6]))
plot(mask(bath, circumpolar),  zlim=c(-10000,-5000), add=T, legend =F, col=rev(bb[7]))
plot(st_geometry(ccamlr_mpa_crop), add=T, border=grey(0.4),lty=3, lwd=1.3)

plot(st_geometry(PF_fronts), add=T, border=1,lty=2, lwd=1)
text(-1600, -3350, "PF",cex=0.9, col=1)
plot(st_geometry(sACCf_fronts), add=T, border=1,lty=1, lwd=1)
text(-600, -2950, "sACCf",cex=0.9, col=1)


plot(st_geometry(ice_shelf),add=T,border=grey(0.3),col=grey(1),lwd=0.9)
plot(st_geometry(st_intersection(land, circumpolar)),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(circumpolar),add=T)
plot(st_geometry(grat),add=T,lty=1, lwd=1)
arctext(x = "INDIAN SECTOR", center = c(0, 0), radius = 4500, middle = 0*pi/2)
arctext(x = "ATLANTIC SECTOR", center = c(0, 0), radius = 4500, middle = 1.2*pi/2)
arctext(x = "PACIFIC SECTOR", center = c(0, 0), radius = 4500, middle = -1.5*pi/2)

plot(st_geometry(st_point_on_surface(ccamlr_mpa_crop)), add=T, 
     col=grey(0), lty=3,pch=as.character(c(9,1,7,3,2,4,5,6,8)), cex=1)



for(i in 1:3) {
  direct = ""
  if(g$degree[i]<0) direct = "W"
  if(g$degree[i]>0) direct = "E"
  arctext(x = paste0(abs(g$degree[i]),'\u00B0',direct), center = c(0, 0), radius = 4500, middle = rad(g$angle_start[i]), cex=0.7)
}

# Islands
# text(-2200, 3250, "South Georgia",cex=0.9, col=island.col,font=4)
island.cex = 0.7
text(-2450, 3400, "South",cex=island.cex, col=island.col,font=4)
text(-2450, 3250, "Georgia",cex=island.cex, col=island.col,font=4)

text(-2200, 2550, "South Orkney",cex=island.cex, col=island.col,font=4)
text(-2200, 2400, "islands",cex=island.cex, col=island.col,font=4)

text(-2800, 1200, "WAP", cex=island.cex, col=island.col,font=4)

text(400, 3800, "Bouvet island", cex=island.cex, col=island.col,font=4)

text(-2400, -100, "Peter I",cex=island.cex, col=island.col,font=4)
text(-2400, -250, "Island",cex=island.cex, col=island.col,font=4)

text(-1250, 3300, "South",cex=island.cex, col=island.col,font=4)
text(-1100, 3150, "Sandwich",cex=island.cex, col=island.col,font=4)
text(-1200, 3000, "Islands",cex=island.cex, col=island.col,font=4)

text(700, -2650, "Balleny",cex=island.cex, col=island.col,font=4)
text(700, -2800, "Islands",cex=island.cex, col=island.col,font=4)

text(-250, -2450, "Scott",cex=island.cex, col=island.col,font=4)
text(-250, -2600, "Island",cex=island.cex, col=island.col,font=4)

text(3850, 1050, "Kerguelen", cex=island.cex, col=island.col,font=4)
text(3850, 900, "Islands", cex=island.cex, col=island.col,font=4)


#Seas
sea.cex = 0.9
text(-1500, 1900, "Weddell", cex=sea.cex, col=sea.col)
text(-1500, 1700, "Sea", cex=sea.cex, col=sea.col)

text(100, 3350, "Lazarev", cex=sea.cex, col=sea.col)
text(100, 3150, "Sea", cex=sea.cex, col=sea.col)

text(2000, 2400, "Cosmonauts", cex=sea.cex, col=sea.col)
text(2000, 2200, "Sea", cex=sea.cex, col=sea.col)

text(2600,  800, "Prydz", cex=sea.cex, col=sea.col)
text(2600,  600, "Bay", cex=sea.cex, col=sea.col)

text(2900,  -1000, "Mawson", cex=sea.cex, col=sea.col)
text(2900,  -1200, "Sea", cex=sea.cex, col=sea.col)

text(2300,  -2000, "D'Urville", cex=sea.cex, col=sea.col)
text(2300,  -2200, "Sea", cex=sea.cex, col=sea.col)

text(0, -1800, "Ross", cex=sea.cex, col=sea.col)
text(0, -2000, "Sea", cex=sea.cex, col=sea.col)

text(-1800, -1500, "Amundsen", cex=sea.cex, col=sea.col)
text(-1800, -1700, "Sea", cex=sea.cex, col=sea.col)

text(-3000,  150, "Bellinghausen", cex=sea.cex, col=sea.col)
text(-3000, -50, "Sea", cex=sea.cex, col=sea.col)

# legend("bottomright", legend = rev(c("5000-10000", "4000-5000", "3000-4000","2000-3000","1000-2000","0-1000")),
#        fill = (bb[2:8]),bty="n", title  ="depth [m]", border = "transparent")
legend(3000, -3000, legend = rev(c("5000 - 10000", "4000 - 5000", "3000 - 4000","2000 - 3000","1000 - 2000","0 - 1000")),
       pt.bg = (bb[2:7]),bty="n", title  ="depth [m]", col = "transparent",pt.cex=2.2,pch=22, text.font =1, cex=0.8)
# plot(-bath,legend.only = T, zlim=c(6000,0), col = blue,
#      axis.args=list(#at=seq(0,10000,3000), #brks,
#                     #label=seq(0,10000,3000),
#                     cex.axis=0.9,tick=T),
#      legend.args=list(text="depth [m]", side=4, font=2, line=3.5, cex=1),
#      legend.width=0.2, legend.shrink=0.75,
#      smallplot=c(0.85,0.87, 0.03,0.17))

dev.off()





