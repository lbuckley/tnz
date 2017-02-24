## CHECK DATA QUALITY

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#Read physiology data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRexpansibility_Buckleyetal.csv")

#read data quality data
setwd(paste(mydir,"MRelevation\\Data\\DataChecking\\", sep=""))
qual=read.csv("TNZ_DataQuality.csv")

#load geographic data
geo=read.csv("TNZ_GeoData.csv")

#load UCT quality data
uct.qual=read.csv("McKechnieetal2016.csv")

#--------------------------------
#Add data

#add quality columns
phy$lon=NA; phy$lat=NA; phy$active=NA; phy$fasted=NA;
phy$capture_quality=NA; phy$feeding=NA; phy$activity=NA; phy$omit=NA

match1= match(phy$Species, qual$Spec.syn)
matched= which(!is.na(match1))

phy$lon[matched]= as.numeric( as.character( qual$Longitude[match1[matched]] ) )
phy$lat[matched]= as.numeric( as.character( qual$Latitude[match1[matched]] ) )
phy$active[matched]= as.character( qual$active[match1[matched]] ) 
phy$fasted[matched]= as.character( qual$fasted[match1[matched]] )
phy$capture_quality[matched]= as.character( qual$capture_quality[match1[matched]] ) 
phy$feeding[matched]= as.character( qual$feeding[match1[matched]] )
phy$activity[matched]= as.character( qual$activity[match1[matched]] )
phy$omit[matched]= as.character( qual$omit[match1[matched]] )

#add geographic data columns
phy$lat.center=NA

match1= match(phy$Species, geo$Species)
matched= which(!is.na(match1))

phy$lat.center[matched]= geo$Latitudinal.center...N.[match1[matched]]

#add UCT quality
phy$uct.qual=NA

match1= match(phy$Species, uct.qual$Species)
matched= which(!is.na(match1))

phy$uct.qual[matched]= as.character( uct.qual$Category[match1[matched]] )

#--------------------------------
# ANALYSIS
## OMIT SPECIES

bird=subset(phy, phy$Taxa=="Bird")
mamm=subset(phy, phy$Taxa=="Mammal")

summary(as.factor(bird$active))
summary(as.factor(bird$fasted))
summary(as.factor(bird$capture_quality))
summary(as.factor(bird$feeding))
summary(as.factor(bird$activity))
summary(as.factor(bird$uct.qual))

#-------------------------------------

#check whether sampling location is within range extent
phy$check= ifelse(phy$UpperLat>phy$lat & phy$LowerLat<phy$lat, 1, 0)

#Calculate distance from measurement to range edge
cold= ifelse( abs(phy$UpperLat)> abs(phy$LowerLat), phy$UpperLat, phy$LowerLat)

phy$dist.cold= abs(cold-phy$lat)
phy$dist.cold.perrange= phy$dist.cold/ (phy$UpperLat - phy$LowerLat  )

#distance to center
phy$dist.center= abs(phy$lat.center-phy$lat)

#-------------------------------------
#write out

setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy,"MRexpansibility_Buckleyetal_wQual.csv")



