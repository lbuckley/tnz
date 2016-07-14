#CODE FROM http://rfunctions.blogspot.com/2013/03/spatial-analysis-split-one-shapefile.html TO SPLIT SHAPEFILES

library(sp)
library(rgdal)

setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MammalsNoShape.csv"  )

#Load physiological data
mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"
setwd(paste(mydir,"MRelevation\\Data\\", sep=""))

#----------------------

wd="C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\Data\\"
wd1= paste(wd, "Shapefiles\\TERRESTRIAL_MAMMALS\\", sep="")

setwd(paste(wd, "Shapefiles\\TERRESTRIAL_MAMMALS\\", sep=""))

# Now load the Mammals' shapefile: All_MAMMALS_NOV2013.shp.
data <- readOGR(choose.files(), "TERRESTRIAL_MAMMALS")
##PICK TERRESTRIAL_MAMMALS.shp

# Ok, we need to choose what we want to separate. Type:
names(data)

# As we want to get one shapefile for each species, we will choose the 'binomial' variable. In that way, we need determine the names and the number of species we are using.

unique <- unique(data@data$binomial)

#--------------------

phy.mam= phy$Spec.syn[phy$Taxa=="Mammal"]

match1= match(phy.mam, unique)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy.mam[not.matched]

spec.ext= phy.mam[matched]

unique= as.character(spec.ext)

#Can't match: "Miniopterus schreibersii"

#----------------------

# Finally, we use a loop to save shapefiles for each species.

for (i in 1:length(unique)) {
  tmp <- data[data$binomial == unique[i], ] 
  writeOGR(tmp, dsn=getwd(), unique[i], driver="ESRI Shapefile",
           overwrite_layer=TRUE)
}

#==========================================================
#EXTRACT TEMPERATURE DATA

setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_all.csv")

phy.noT= phy[is.na(phy$Tmin) | is.na(phy$Tmax), ]


#========================================================
#MAMMALS

phy=phy.noT[phy.noT$Taxa=="Mammal"  ,]
#--------------------------------------------
#Load shape files

#specify climate data
clim.min= clim.pmin
clim.max= clim.pmax

sdir= paste(wd, "data\\Shapefiles\\TERRESTRIAL_MAMMALS\\", sep="")
setwd(sdir)

#get list of polygon files
speciesfiles<-list.files(sdir,pattern="\\.shp$")#get shape files
speciesnames<- gsub(".shp", "", speciesfiles)

#Match names
match1= match(phy$Spec.syn,speciesnames)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

species.matched=as.character(phy$Spec.syn[!is.na(match1)])

#---------------------------------------------

#Define function to extract Tmin, Tmax
# TEST: species=species.matched[3]

#function to apply to each polygon

#---------------------------

#TminTmaxFun(species.matched[3])

#Run extraction function
output<-ldply(species.matched,TminTmaxFun)

#write out
wd="C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"
write.csv(output, paste(wd, "OUT\\MammalTminTmax_add.csv", sep=""), row.names=F)

#========================================================
#BIRDS

phy=phy.noT[phy.noT$Taxa=="Bird"  ,]
#--------------------------------------------
sdir= paste(wd, "Data\\Shapefiles\\Birds\\", sep="")
setwd(sdir)

#get list of polygon files
speciesfiles<-list.files(sdir,pattern="\\.shp$")#get shape files
speciesnames<- gsub(".shp", "", speciesfiles)
speciesnames.match<- gsub("_", " ", speciesnames)
speciesnames.match<- substr(speciesnames.match, 1, nchar(speciesnames.match)-9)

#------------------
#Match names
match1= match(as.character(phy$Spec.syn), speciesnames.match )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
#phy$Spec.syn[not.matched]

species.matched=as.character(phy$Spec.syn[matched])
#convert to file names
match1= match(species.matched, speciesnames.match )
species.filename= speciesnames[match1]

#---------------------------

#output3=TminTmaxFun(species.filename[2])

#Run extraction function
#output<-ldply(species.filename[1:length(species.filename)],TminTmaxFun)
#output<-ldply(species.filename[2:50],TminTmaxFun)
output<-ldply(species.filename,TminTmaxFun)

output$genspec= species.matched
output$shapename= species.filename

#write out
write.csv(output, paste(wd, "OUT/BirdTminTmax_add.csv", sep=""), row.names=F)

#====================================
#Combine Tmax Tmin

setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_all.csv")

#add species
setwd(paste(mydir,"\\Out\\", sep=""))
phy.addb=read.csv("BirdTminTmax_add.csv")
phy.addm=read.csv("MammalTminTmax_add.csv")

#birds
match1= match(as.character(phy.addb$genspec) , as.character(phy$Spec.syn) )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy$Tmin[match1]<- phy.addb$Tmin
phy$Tmedian.min[match1]<- phy.addb$Tmedian.min
phy$T5q.min[match1]<- phy.addb$T5q.min
phy$T10q.min[match1]<- phy.addb$T10q.min
phy$Tsd.min[match1]<- phy.addb$Tsd.min

phy$Tmax[match1]<- phy.addb$Tmax
phy$Tmedian.max[match1]<- phy.addb$Tmedian.max
phy$T5q.max[match1]<- phy.addb$T5q.max
phy$T10q.max[match1]<- phy.addb$T10q.max
phy$Tsd.max[match1]<- phy.addb$Tsd.max

phy$UpperLat[match1]<- phy.addb$UpperLat
phy$LowerLat[match1]<- phy.addb$LowerLat
#------------------------------------------
#mammals
match1= match(as.character(phy.addm$genspec) , as.character(phy$Spec.syn) )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy$Tmin[match1]<- phy.addm$Tmin
phy$Tmedian.min[match1]<- phy.addm$Tmedian.min
phy$T5q.min[match1]<- phy.addm$T5q.min
phy$T10q.min[match1]<- phy.addm$T10q.min
phy$Tsd.min[match1]<- phy.addm$Tsd.min

phy$Tmax[match1]<- phy.addm$Tmax
phy$Tmedian.max[match1]<- phy.addm$Tmedian.max
phy$T5q.max[match1]<- phy.addm$T5q.max
phy$T10q.max[match1]<- phy.addm$T10q.max
phy$Tsd.max[match1]<- phy.addm$Tsd.max

phy$UpperLat[match1]<- phy.addm$UpperLat
phy$LowerLat[match1]<- phy.addm$LowerLat

#-----------------------------------------
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy, "MRelevation_all.csv")

