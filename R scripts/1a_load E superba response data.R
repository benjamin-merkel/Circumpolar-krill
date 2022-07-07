library(readxl)

# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 


# load krill base -----
krill                     <- read.csv("data/response/krillbase.csv", header=T)
krill                     <- krill[order(krill$Season),]
krill$presence_absence    <- 0
krill$presence_absence[krill$No..ofkrill.under.1m2>0.1] <- 1
krill$date                <- as.Date(strptime(krill$Date,"%d-%b-%Y"))
krill$month               <- as.numeric(strftime(krill$date,"%m"))
krill$year                <- as.numeric(strftime(krill$date, "%Y"))
krill$doy                 <- as.numeric(strftime(krill$date, "%j"))
krill$source              <- krill$Source
krill$abundance           <- krill$Standardisedkrillunder.1m2
krill$longitude           <- krill$Longitude
krill$latitude            <- krill$Latitude
krill$country             <- str_split_fixed(krill$source," ",3)[,1]
krill$depthmeters.sampled <- krill$Bottomsamplingdepth..m. - krill$Topsamplingdepth..m.
krill$data                <- "KRILLBASE"

krill$keep                <- "Y"
# remove absence data with starting depth deeper than 50m and depth meters sampled less than 100m
krill$keep[krill$Topsamplingdepth..m. > 50  & krill$presence_absence==0] <- "N"
krill$keep[krill$depthmeters.sampled  < 100 & krill$presence_absence==0] <- "N"

# remove absences during day during summer following Siegel 2016
krill$keep[krill$Day.Night == "day" & krill$presence_absence==0 & !krill$month %in% c(12,1:3)] <- "N"

# remove LAKRIS expedition data as these were provided by AWI
krill$keep[krill$source == "German data (LAKRIS cruise), V. Siegel, pers. comm."] <- "N"

krill <- krill[krill$keep=="Y",]


# load LAKRIS data from Bettina Meyer, AWI -----
lakris           <- data.frame(read_excel("data/response/from AWI/Lakris_Esuperba.xlsx", sheet = "Lakris_Esuperba"))
lakris$date      <- as.Date(strptime(lakris$STATDATUM,"%Y%m%d"))
lakris$year      <- as.numeric(strftime(lakris$date, "%Y"))
lakris$month     <- as.numeric(strftime(lakris$date, "%m"))
lakris$doy       <- as.numeric(strftime(lakris$date, "%j"))
lakris$longitude <- as.numeric(lakris$GLSTART)
lakris$latitude  <- as.numeric(lakris$GBSTART)
lakris$GESAMTSTCK[is.na(lakris$GESAMTSTCK)] <- 0

lakris2           <- lakris[!duplicated(paste(lakris$STATION,lakris$REISENR,lakris$FANGGERAET)),]
lakris2$abundance <- round(lakris2$GESAMTSTCK/abs(lakris2$FILVOL),3)
lakris2$subadult  <- lakris2$adult <- lakris2$juv <- 0
for(i in 1:nrow(lakris2)){
  l2 <- lakris[lakris$STATION==lakris2$STATION[i] & lakris$REISENR==lakris2$REISENR[i] & lakris$FANGGERAET==lakris2$FANGGERAET[i],]
  if(nrow(l2[l2$STAGE=="adult",])>0)     lakris2$adult[i]    <- lakris2$abundance[i] * l2$SumOfLANZAHL[l2$STAGE=="adult"]   /sum(l2$SumOfLANZAHL)
  if(nrow(l2[l2$STAGE=="subadult",])>0)  lakris2$subadult[i] <- lakris2$abundance[i] * l2$SumOfLANZAHL[l2$STAGE=="subadult"]/sum(l2$SumOfLANZAHL)
  if(nrow(l2[l2$STAGE=="juv",])>0)       lakris2$juv[i]      <- lakris2$abundance[i] * l2$SumOfLANZAHL[l2$STAGE=="juv"]     /sum(l2$SumOfLANZAHL)
}
lakris2$subadult[is.na(lakris2$subadult)]   <- lakris2$adult[is.na(lakris2$subadult)] <- lakris2$juv[is.na(lakris2$subadult)] <- lakris2$abundance[is.na(lakris2$subadult)] <- 0
lakris2$presence_absence                    <- 0
lakris2$presence_absence[lakris2$adult > 0] <- 1
lakris.subset                               <- lakris2[!(lakris2$LICHT=="B" & lakris2$presence_absence==0),] # D= DARK, B= Beleuchtet und T= Twilight
lakris.subset$source                        <- lakris.subset$REISENR
lakris.subset$data                          <- "LAKRIS"



# load old Thunen data from Katharina Teschke, AWI -----

Biota       <- read_excel("data/response/from Thünen/euphausidae_siegel_vers3.xlsx", sheet = "Biota")
Biota$merge <- paste(Biota$Cruise, Biota$Station, Biota$Sample)
Sample      <- read_excel("data/response/from Thünen/euphausidae_siegel_vers3.xlsx", sheet = "Sample")
Sample$merge<- paste(Sample$Cruise, Sample$Station, Sample$Sample)

teschke             <- merge(Biota, Sample, by="merge")
teschke$abundance   <- round(teschke$Abundance,3)
teschke$date        <- as.Date(teschke$Date)
teschke$year        <- as.numeric(strftime(teschke$date, "%Y"))
teschke$month       <- as.numeric(strftime(teschke$date, "%m"))
teschke$doy         <- as.numeric(strftime(teschke$date, "%j"))
teschke$starthour   <- as.numeric(substr(teschke$StartTime,1,2)) + as.numeric(substr(teschke$StartTime,4,5))/60
teschke$presence_absence <- 0
teschke$presence_absence[teschke$abundance>0] <- 1
teschke$reisenr     <- 0
teschke$reisenr[teschke$Cruise.x=="ANT-XXI/4"]   <- "21-4"
teschke$reisenr[teschke$Cruise.x=="ANT-XXIII/2"] <- "23-2"
teschke$reisenr[teschke$Cruise.x=="ANT-XXIII/6"] <- "23-6"
teschke$reisenr[teschke$Cruise.x=="ANT-XXIV/2"]  <- "24-2"

teschke$source      <- teschke$Cruise.x
teschke$longitude   <- as.numeric(teschke$LonStart)
teschke$latitude    <- as.numeric(teschke$LatStart)
teschke.subset      <- teschke[teschke$ScientificName=="Euphausia superba" & teschke$reisenr=="0"  & teschke$LifeStage=="adult",]
teschke.subset$data <- "THUNEN"

# NAre data ----
nare                  <- read.csv("data/response/Falk-Petersen 1999/Falk-Petersen 1999.csv")
nare$date             <- as.Date(strptime(nare$date,"%d.%m.%y"))
nare$year             <- as.numeric(strftime(nare$date, "%Y"))
nare$month            <- as.numeric(strftime(nare$date, "%m"))
nare$doy              <- as.numeric(strftime(nare$date, "%j"))
nare$presence_absence <- nare$E.superba
nare$source           <- "MV Polarsirkel 1993 (NARE)"
nare$abundance        <- 0
nare$data             <- "NARE"


# CHINARE data -----
Sys.setlocale("LC_ALL","English")
chinare           <- read_excel("data/response/CHINARE/Circumpolar data of Antarctic krill and zooplankton from CHINARE2013-2014.xlsx")
chinare$source    <- "CHINARE 2013/14"
chinare$abundance <- chinare$`Es adult`
chinare$date      <- as.Date(strptime(chinare$`Sampling date`,"%d-%b-%Y"))
chinare$year      <- as.numeric(strftime(chinare$date, "%Y"))
chinare$month     <- as.numeric(strftime(chinare$date, "%m"))
chinare$doy       <- as.numeric(strftime(chinare$date, "%j"))
chinare$presence_absence <- chinare$abundance
chinare$presence_absence[chinare$presence_absence>0] <- 1
chinare$longitude <- chinare$Longitude
chinare$latitude  <- chinare$Latitude
chinare$data      <- "CHINARE"


# GBIF ----
gbif                  <- read.csv('data/response/GBIF/Antkrill/occurrence.csv',header=T,sep=";")
gbif$decimalLongitude <- as.numeric(gbif$decimalLongitude)
gbif$decimalLatitude  <- as.numeric(gbif$decimalLatitude)
gbif                  <- gbif[!is.na(gbif$decimalLatitude) & !is.na(gbif$decimalLongitude) & gbif$decimalLatitude > -85 & gbif$decimalLatitude < -40,]
gbif$longitude        <- as.numeric(gbif$decimalLongitude)
gbif$latitude         <- as.numeric(gbif$decimalLatitude)
gbif                  <- gbif[!gbif$lifeStage %in% c("LARVA","JUVENILE"),]
gbif                  <- gbif[!gbif$datasetName %in% c("Southern Ocean Continuous Zooplankton Recorder (SO-CPR) Survey", 
                                                       "KRILLBASE_antarctickrillssalpsdensity_1926_2016"),]
gbif$presence_absence <- NA
gbif$presence_absence[gbif$individualCount==0]                               <- 0
gbif$presence_absence[tolower(gbif$occurrenceStatus) %in% c('absent')]       <- 0
gbif$presence_absence[tolower(gbif$occurrenceStatus) %in% c('present')]      <- 1
gbif$presence_absence[gbif$scientificName != "Euphausia superba Dana, 1850"] <- 0
gbif$source          <- "gbif"
gbif$date            <- as.Date(strptime(gbif$eventDate,"%Y-%m-%d"))
gbif$year            <- as.numeric(strftime(gbif$date, "%Y"))
gbif$month           <- as.numeric(strftime(gbif$date, "%m"))
gbif$doy             <- as.numeric(strftime(gbif$date, "%j"))
gbif                 <- gbif[gbif$hasGeospatialIssue==F,]
gbif                 <- gbif[gbif$year > 1970,]
gbif                 <- gbif[!is.na(gbif$decimalLatitude) & !is.na(gbif$decimalLongitude),]
gbif                 <- gbif[gbif$datasetName=="Antarctic Euphausiacea occurrence data from ITALICA 2000 Expedition",]#"Ross Sea Biodiversity Survey 2004 (BioRoss)",]#"Antarctic Euphausiacea occurence data from Norwegian Antarctic Research Expedition 1976-77",]
gbif$data            <- "ITALICA"
gbif$abundance       <- NA
gbif$presence_absence[gbif$scientificName == "Euphausia crystallorophias Holt & Tattersall, 1906"] <- 0


# KPH cruise 2019 ----
kph.trawl           <- read.csv("data/response/KPH cruise 2019/krill stations from cruise report.csv",sep=";")
kph.trawl           <- kph.trawl[!is.na(kph.trawl$start.lat),]
kph.trawl$latitude  <- kph.trawl$start.lat
kph.trawl$longitude <- kph.trawl$start.lon
kph.trawl$year      <- 2019
kph.trawl$date      <- strptime(kph.trawl$start.date, "%d/%m/%Y")
kph.trawl$source    <- kph.trawl$data <- "KPH2019"
kph.trawl$abundance <- NA
kph.trawl$doy       <- as.numeric(strftime(kph.trawl$date,"%j"))
kph.trawl$month     <- as.numeric(strftime(kph.trawl$date,"%m"))
kph.trawl$presence_absence <- kph.trawl$E.sup.registered



# save and combine ----
picked.columns <- c("longitude","latitude","presence_absence","abundance","year","month","doy","source","data")
krill.combined <- rbind(krill[,picked.columns], 
                        teschke.subset[,picked.columns], 
                        lakris.subset[,picked.columns], 
                        nare[,picked.columns],
                        chinare[,picked.columns],
                        gbif[,picked.columns],
                        kph.trawl[,picked.columns])
krill.combined$latitude  <- as.numeric(krill.combined$latitude)
krill.combined$longitude <- as.numeric(krill.combined$longitude)
krill.combined <- krill.combined[!duplicated(paste(krill.combined$longitude,
                                                   krill.combined$latitude,
                                                   krill.combined$year, 
                                                   krill.combined$doy, 
                                                   krill.combined$abundance)),]


krill.combined <- st_as_sf(x = krill.combined, coords = c("longitude", "latitude"),  crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
krill.combined <- st_transform(krill.combined, south_pole_equal_area.proj)
krill.combined <- krill.combined[!is.na(krill.combined$presence_absence),]
krill.combined <- krill.combined[!duplicated(paste(krill.combined$geometry, krill.combined$presence_absence, krill.combined$year, krill.combined$data)),]


saveRDS(krill.combined,"data/E.superba response data.rds")
