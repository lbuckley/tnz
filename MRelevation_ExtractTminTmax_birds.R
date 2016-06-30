#########################################
## TNZ analysis
# Project ranges in current and future climates
#BIRDS

#Calculate Tmin and Tmax
#BASED ONLY ON EDGES
#########################################
mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"
wd="C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#libraries needed
library(raster)
library(maptools)
library(proj4)
library(Hmisc)
library(plyr)
library(rgeos)
library(rgdal)
library(dismo)
library(foreach)

#First last of vector
firstlast= function(x){  
  fl= c(NA, NA)
  x=na.omit(x)
  if(length(x)>0) fl= c( x[1], x[length(x)])
  return(fl)
}

#========================
# Load climate data
#USE WORLD CLIM PROJECTIONS
#5 Minutes

# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month

# 2100 RCP6

#current data
clim.p= getData('worldclim', var='bio', res=5)
clim.pmax= clim.p$bio5/10
clim.pmin= clim.p$bio6/10

#future data
clim.f=getData('CMIP5', var='bio', res=5, rcp=60, model='BC', year=70)
clim.fmax= clim.f$bc60bi705/10
clim.fmin= clim.f$bc60bi706/10
#'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO".
#'rcp' should be one of 26, 45, 60, or 85.
#'year' should be 50 or 70

plot(clim.fmax-clim.pmax)

#-----------------------------------------
#Load physiological data

setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_all.csv")
#phy=na.omit(phy)

##ESTRACT FILES FOR ADDITIONAL SPECIES
phy=read.csv("Phy_wTraits.csv")

phy=phy[phy$Taxa=="Bird",]

#Use synonyms to match up
phy$Spec.syn= as.character(phy$Species)
phy$Spec.syn[which(phy$Spec.syn=="Anas rhynchotis")]="Spatula rhynchotis"
phy$Spec.syn[which(phy$Spec.syn=="Collocalia vanikorensis")]="Aerodramus vanikorensis"
phy$Spec.syn[which(phy$Spec.syn=="Aceros plicatus")]="Rhyticeros plicatus"
phy$Spec.syn[which(phy$Spec.syn=="Coturnix chinensis")]="Synoicus chinensis"
phy$Spec.syn[which(phy$Spec.syn=="Aramides cajanea")]="Aramides cajaneus"
phy$Spec.syn[which(phy$Spec.syn=="Gallinula mortierii")]="Tribonyx mortierii"
phy$Spec.syn[which(phy$Spec.syn=="Gallinula ventralis")]="Tribonyx ventralis"
phy$Spec.syn[which(phy$Spec.syn=="Gallirallus owstoni")]="Hypotaenidia owstoni"
phy$Spec.syn[which(phy$Spec.syn=="Porzana cinerea")]="Amaurornis cinerea"
phy$Spec.syn[which(phy$Spec.syn=="Cacatua roseicapilla")]="Eolophus roseicapilla"
phy$Spec.syn[which(phy$Spec.syn=="Otus leucotis")]="Ptilopsis leucotis"
phy$Spec.syn[which(phy$Spec.syn=="Poecile atricapillus")]="Parus atricapillus"
phy$Spec.syn[which(phy$Spec.syn=="Junco hyemalis hyemalis")]="Junco hyemalis"
phy$Spec.syn[which(phy$Spec.syn=="Junco hyemalis oreganus")]="Junco hyemalis"
phy$Spec.syn[which(phy$Spec.syn=="Hesperiphona vespertina")]="Coccothraustes vespertinus"

#--------------------------------------------
#Load shape files

#specify climate data
clim.min= clim.pmin
clim.max= clim.pmax

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

#sdir= paste(wd, "Data\\Shapefiles\\", sep="")
#setwd(sdir)
#write.csv(speciesnames, "BirdSpecies.csv")

#---------------------------------------------

#Define function to extract Tmin, Tmax
# TEST: species=species.filename[15]

#function to apply to each polygon
TminTmaxFun<-function(species){

  speciesdata<-readOGR(dsn=".",layer=species)
  
  #define current projection of polygon )if needed
  proj4string(speciesdata) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"#(as per readOGR)
      out<-speciesdata
  
  #for upper and lower lat
  latlon<-spTransform(out,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  LowerLat<-summary(latlon)$bbox[2,1]#lower latitude - full extent of distribution data
  UpperLat<-summary(latlon)$bbox[2,2]#upper latitude - full extent of distribution data
  
  #EXTRACT MIN, MAX
  for(clim.k in 1:2){
  if(clim.k==1) bio1= clim.min
  if(clim.k==2) bio1= clim.max
    
  #convert species polygon to projection of the climatic raster file
  proj.xy<-spTransform(out,CRS(projection(bio1)))
  #biggest polygon
  #areas= proj.xy$shape_Area
  #shape= proj.xy[which.max(areas),]
  clip= raster::mask(bio1,proj.xy)
  
  #clump
  sc= clump(clip, directions=4)
  
  ## calculate frequency of each clump/patch
  sc.freq <- as.data.frame(freq(sc))
  sc.freq<- sc.freq[!is.na(sc.freq$value),]
  
  rmID <- sc.freq$value[which.max(sc.freq$count)]
  clip[!sc %in% rmID] <- NA
  #clip[sc %in% rmID] <- 1 
  
  #trim
  clip= trim(clip)
  
  #find max and min across each column
  #flip distribution and get first and last value each row
  edges= foreach(i=1:ncol(clip)) %do% firstlast(getValues(t(clip$bio6),i)) 
  edges1= t(as.data.frame(edges)) #one column warm edges, 1 column cold edges
  
  ##extract grid cells at edge
  #bound <- as(proj.xy, "SpatialLines", sum=FALSE) 
  #clip3= extract(bio1, bound)
  ##lengths
  #SpatialLinesLengths(bound, longlat=TRUE)
  
  ##choose biggest section of range
  #vals= clip3[[which.max(lengths(clip3))]]
  #ncells= floor(length(vals)/20) #average across 5% of cells
  
  ##get climatic grid data where it overlaps the species distribution
  #clip<-extract(bio1,proj.xy, weights=T,cellnumbers=TRUE)
  #clip<-ldply(clip)
  #clip <-ddply(clip,.(cell),summarise,value=mean(value),weight=sum(weight))
  #clip<-subset(clip,!is.na(clip[,2]))
  #vals<-clip$value
  #ws<-clip$weight
  
  #coldest / warmest cell
  if(clim.k==1) {
    vals= edges1[,which.min(colMeans(edges1))] #cold edge
    Tmin= min(vals,na.rm=T)
    Tmedian.min= median(vals,na.rm=T)
    T5q.min= quantile(vals, 0.05, na.rm=T)
    T10q.min=quantile(vals, 0.10, na.rm=T)
    Tsd.min= sd(vals,na.rm=T)
  }
  if(clim.k==2){
    vals= edges1[,which.max(colMeans(edges1))] #warm edge
    Tmax= max(vals,na.rm=T)
    Tmedian.max= median(vals,na.rm=T)
    T5q.max= quantile(vals, 0.95, na.rm=T)
    T10q.max=quantile(vals, 0.90, na.rm=T)
    Tsd.max= sd(vals,na.rm=T)
  }
  NumberGrids<-length(vals[!is.na(vals)])
  Area<-gArea(proj.xy)
  } #end clim loop
  
  data.frame(species,UpperLat,LowerLat,Tmin, Tmedian.min, T5q.min,T10q.min,Tsd.min, Tmax, Tmedian.max, T5q.max,T10q.max,Tsd.max, NumberGrids,Area,stringsAsFactors=F)
  
} #end TminTmaxFun
#---------------------------

#TminTmaxFun(species.filename[2])

#Run extraction function
#output<-ldply(species.filename[1:100],TminTmaxFun)
output<-ldply(species.filename[2:length(species.filename)],TminTmaxFun)

output$genspec= species.matched

#write out
write.csv(output, paste(wd, "OUT/BirdTminTmax.csv", sep=""), row.names=F)



 