#CODE FROM http://rfunctions.blogspot.com/2013/03/spatial-analysis-split-one-shapefile.html TO SPLIT SHAPEFILES

library(sp)
library(rgdal)

#Load physiological data
mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"
setwd(paste(mydir,"MRelevation\\Data\\", sep=""))

##ESTRACT FILES FOR ADDITIONAL SPECIES
#phy=read.csv("MRelevation_wTraits2.csv")
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
#phy=read.csv("MRelevation_all.csv")
phy=read.csv("Phy_wTraits.csv")
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
#Select species from physiological data
phy.spec= as.character(phy$Species)
phy.spec[which(phy.spec=="Galerella sanguinea")]="Herpestes sanguineus"
phy.spec[which(phy.spec=="Cryptomys bocagei")]="Fukomys bocagei"
phy.spec[which(phy.spec=="Cryptomys damarensis")]="Fukomys damarensis"
phy.spec[which(phy.spec=="Cryptomys mechowi")]="Fukomys mechowi"
phy.spec[which(phy.spec=="Marmosa microtarsus")]="Gracilinanus microtarsus"
phy.spec[which(phy.spec=="Ningaui yvonnae")]="Ningaui yvonneae"

phy.spec[which(phy.spec=="Echymipera kalabu")]="Echymipera kalubu"
phy.spec[which(phy.spec=="Eptesicus vulturnus")]="Vespadelus vulturnus"
#phy.spec[which(phy.spec=="Miniopterus schreibersii")]=""
phy.spec[which(phy.spec=="Anoura caudifera")]="Anoura caudifer"
phy.spec[which(phy.spec=="Choeroniscus perspicillata")]="Carollia perspicillata" #?
phy.spec[which(phy.spec=="Artibeus literatus")]="Artibeus lituratus"
#phy.spec[which(phy.spec=="Pteronotus parnelli")]=""
phy.spec[which(phy.spec=="Rhinonicteris aurantius")]="Rhinonicteris aurantia"
phy.spec[which(phy.spec=="Clethrionomys rufocanus")]="Myodes rufocanus"
phy.spec[which(phy.spec=="Clethrionomys rutilus")]="Myodes rutilus"
phy.spec[which(phy.spec=="Cricetulus triton")]="Tscherskia triton"
phy.spec[which(phy.spec=="Apodemus hermonensis")]="Apodemus witherbyi"
#phy.spec[which(phy.spec=="Mormota flaviventris")]=""
phy.spec[which(phy.spec=="Nannospalax leucodon")]="Spalax leucodon"
#phy.spec[which(phy.spec=="Cavia porcellus")]=""
phy.spec[which(phy.spec=="Procavia johnstoni")]="Procavia capensis"
phy.spec[which(phy.spec=="Hemiechinus aethiopicus")]="Paraechinus aethiopicus"
#phy.spec[which(phy.spec=="Equus caballus")]=""
phy.spec[which(phy.spec=="Ovis orientalis aries")]="Ovis orientalis"
phy.spec[which(phy.spec=="Capra hircus")]="Capra aegagrus"
phy.spec[which(phy.spec=="Cervus elaphus c.")]="Cervus elaphus"
phy.spec[which(phy.spec=="Alopex lagopus")]="Vulpes lagopus"
phy.spec[which(phy.spec=="Proteles cristatus")]="Proteles cristata"
phy.spec[which(phy.spec=="Callithrix pygmaea")]="Cebuella pygmaea"
phy.spec[which(phy.spec=="Dipodillus dasyurus")]="Gerbillus dasyurus"
phy.spec[which(phy.spec=="Micaelamys namaquensis")]="Aethomys namaquensis"

phy.mam= phy.spec[phy$Taxa=="Mammal"]

match1= match(phy.mamm, unique)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy.mam[not.matched]

spec.ext= phy.mam[matched]

unique= spec.ext

#Can't match: "Miniopterus schreibersii"

#----------------------

# Finally, we use a loop to save shapefiles for each species.

for (i in 173:length(unique)) {
  tmp <- data[data$binomial == unique[i], ] 
  writeOGR(tmp, dsn=getwd(), unique[i], driver="ESRI Shapefile",
           overwrite_layer=TRUE)
}

