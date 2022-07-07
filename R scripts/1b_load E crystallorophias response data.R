library(sf)
library(readxl)
library(sp)

# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 


sf  <- read.csv('data/response/GBIF/Antkrill/occurrence.csv',header=T,sep=";")
sf$decimalLongitude <- as.numeric(sf$decimalLongitude)
sf$decimalLatitude <- as.numeric(sf$decimalLatitude)
sf2  <- sf[!is.na(sf$decimalLatitude) & !is.na(sf$decimalLongitude) & sf$decimalLatitude > -85 & sf$decimalLatitude < -40,]
sf2$longitude <- as.numeric(sf2$decimalLongitude)
sf2$latitude  <- as.numeric(sf2$decimalLatitude)
sf2 <- sf2[!sf2$lifeStage %in% c("LARVA","JUVENILE"),]
sf2$data <- "GBIF"
sf2$date <- as.Date(strptime(sf2$eventDate,"%Y-%m-%d"))
sf2$year <- as.numeric(strftime(sf2$date, "%Y"))
sf2$month <- as.numeric(strftime(sf2$date, "%m"))
sf2$doy <- as.numeric(strftime(sf2$date, "%j"))
sf2 <- sf2[sf2$hasGeospatialIssue==F,]
sf2 <- sf2[sf2$year > 1970,]
sf2 <- sf2[!is.na(sf2$decimalLatitude) & !is.na(sf2$decimalLongitude),]

# sf <- st_as_sf(x = sf2[sf2$month %in% c(12,1:4),], coords = c("longitude", "latitude"),  crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
# sf <- st_as_sf(x = sfc[sfc$month %in% c(12,1:4),], coords = c("longitude", "latitude"),  crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
# sf <- st_transform(sf, south_pole_equal_area.proj)

sf3 <- sf2[sf2$datasetName=="Antarctic Euphausiacea occurrence data from ITALICA 2000 Expedition",]#"Ross Sea Biodiversity Survey 2004 (BioRoss)",]#"Antarctic Euphausiacea occurence data from Norwegian Antarctic Research Expedition 1976-77",]
sf3$source <- "ITALICA"
sf3$abundance <- NA
sf3$presence_absence <- 0
sf3$presence_absence[sf3$scientificName == "Euphausia crystallorophias Holt & Tattersall, 1906"] <- 1
sf3 <- sf3[sf3$scientificName != "Euphausia crystallorophias Holt & Tattersall, 1906",] 

sf4 <- sf2[sf2$scientificName == 'Euphausia crystallorophias Holt & Tattersall, 1906',]
sf4 <- sf4[sf4$datasetName !="Antarctic Euphausiacea occurrence data from ITALICA 2000 Expedition",]

sf4$presence_absence <- NA
sf4$presence_absence[sf4$individualCount==0] <- 0
sf4$presence_absence[sf4$occurrenceStatus %in% c('absent', 'Absent')] <- 0
sf4$presence_absence[sf4$individualCount>0] <- 1
sf4$presence_absence[sf4$occurrenceStatus %in% c('present', 'Present')] <- 1
sf4$data = "GBIF"


sf  <- read.csv('data/response/GBIF/Icekrill/occurrence.csv',header=T)
sf$decimalLongitude <- as.numeric(sf$decimalLongitude)
sf$decimalLatitude  <- as.numeric(sf$decimalLatitude)
sf2 <- sf[!is.na(sf$decimalLatitude) & sf$decimalLatitude > -80 & sf$decimalLatitude < -40,]
sf2$longitude <- as.numeric(sf2$decimalLongitude)
sf2$latitude  <- as.numeric(sf2$decimalLatitude)
sf2 <- sf2[!sf2$lifeStage %in% c("LARVA","JUVENILE"),]
sf2$presence_absence <- NA
sf2$presence_absence[sf2$individualCount==0] <- 0
sf2$presence_absence[sf2$occurrenceStatus %in% c('absent', 'Absent')] <- 0
sf2$presence_absence[sf2$individualCount>0] <- 1
sf2$presence_absence[sf2$occurrenceStatus %in% c('present', 'Present')] <- 1
sf2$data <- "GBIF"
sf2$date <- as.Date(strptime(sf2$eventDate,"%Y-%m-%d"))
sf2$year <- as.numeric(strftime(sf2$date, "%Y"))
sf2$month <- as.numeric(strftime(sf2$date, "%m"))
sf2$doy <- as.numeric(strftime(sf2$date, "%j"))
sf2$source <- sf2$bibliographicCitation

# load thunen data

Biota       <- read_excel("data/Thünen krill/euphausidae_siegel_vers3.xlsx", sheet = "Biota")
Biota$merge <- paste(Biota$Cruise, Biota$Station, Biota$Sample)
Sample      <- read_excel("data/Thünen krill/euphausidae_siegel_vers3.xlsx", sheet = "Sample")
Sample$merge<- paste(Sample$Cruise, Sample$Station, Sample$Sample)

siegel <- merge(Biota, Sample, by="merge")
siegel$date  <- as.Date(siegel$Date)
siegel$year  <- as.numeric(strftime(siegel$date, "%Y"))
siegel$month <- as.numeric(strftime(siegel$date, "%m"))
siegel$doy <- as.numeric(strftime(siegel$date, "%j"))
siegel$presence_absence <- 0
siegel$presence_absence[siegel$Abundance>0] <- 1
siegel$source <- siegel$Cruise.x
siegel$abundance <- siegel$Abundance
siegel$longitude <- as.numeric(siegel$LonStart)
siegel$latitude <- as.numeric(siegel$LatStart)
siegel.subset <- siegel[siegel$ScientificName=="Euphausia crystallorophias" & 
                          siegel$LifeStage=="adult",]
siegel.subset$data <- "Thunen"

#siegel.krill <- siegel.krill[siegel.krill$doy >= 363 | siegel.krill$doy <=117,]


# NAre data
nare       <- read.csv("data/Falk-Petersen 1999/Falk-Petersen 1999.csv")
#nare       <- st_as_sf(nare, coords = c("longitude", "latitude"),  crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
nare$date  <- as.Date(strptime(nare$date,"%d.%m.%y"))
nare$year  <- as.numeric(strftime(nare$date, "%Y"))
nare$month <- as.numeric(strftime(nare$date, "%m"))
nare$doy   <- as.numeric(strftime(nare$date, "%j"))
nare$presence_absence <- nare$E.crystallorophias
#nare       <- st_transform(nare, south_pole_equal_area.proj)
nare$source<- "MV Polarsirkel 1993 (NARE)"
nare$data <- "NARE"


# KPH cruise 2019 ----

kph.trawl <- read.csv("data/KPH cruise 2019/krill stations from cruise report.csv",sep=";")
kph.trawl <- kph.trawl[!is.na(kph.trawl$start.lat),]
# kph.trawl <- st_as_sf(x = kph.trawl, coords = c("start.lon", "start.lat"), 
#                       crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
# kph.trawl <- st_transform(kph.trawl, south_pole_equal_area.proj)
kph.trawl$latitude <- kph.trawl$start.lat
kph.trawl$longitude <- kph.trawl$start.lon
kph.trawl$year <- 2019
kph.trawl$date <- strptime(kph.trawl$start.date, "%d.%m.%Y")
kph.trawl$source <- kph.trawl$data <- "KPH2019"
kph.trawl$abundance <- NA
kph.trawl$doy <- as.numeric(strftime(kph.trawl$date,"%j"))
kph.trawl$month <- as.numeric(strftime(kph.trawl$date,"%m"))
kph.trawl$presence_absence <- kph.trawl$E.cry.registered


# save and combine

picked.columns <- c("longitude","latitude","presence_absence","year","month","doy","source","data")
krill.combined <- rbind(kph.trawl[,picked.columns],sf2[,picked.columns],sf3[,picked.columns],sf4[,picked.columns],siegel.subset[,picked.columns], nare[,picked.columns])
krill.combined <- st_as_sf(x = krill.combined, coords = c("longitude", "latitude"),  crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
krill.combined <- st_transform(krill.combined, south_pole_equal_area.proj)
krill.combined <- krill.combined[!is.na(krill.combined$presence_absence),]
krill.combined <- krill.combined[!duplicated(paste(krill.combined$geometry, krill.combined$presence_absence, krill.combined$year, krill.combined$data)),]

saveRDS(krill.combined,"data/E.crystallorophias response data.rds")
