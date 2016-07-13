#########################################
## TNZ analysis
# Project ranges in current and future climates
#MAMMALS

#Calculate Tmin and Tmax
#########################################

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

#---------------------------------------
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

setwd(paste(wd,"MRelevation\\Out\\", sep=""))
#phy=read.csv("MRelevation_all.csv")
phy=read.csv("Phy_all.csv")
#phy=na.omit(phy)
phy= phy[phy$Taxa=="Mammal",]

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

#get list of saved shapefiles 
speciesfiles.rda<-list.files(sdir,pattern="\\.rda$")#get shape files
speciesnames.rda<- gsub("\\.rda", "", speciesfiles.rda)

#Match names
match1= match(phy$Spec.syn,speciesnames)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

species.matched=as.character(phy$Spec.syn[!is.na(match1)])

#---------------------------------------------

#Define function to extract Tmin, Tmax
# TEST: species=species.matched[1]

#function to apply to each polygon
TminTmaxFun<-function(species){
  print(species)
  
  #check if shapefile saved
  if(species %in% speciesnames.rda) speciesdat= load(paste(species,".rda",sep="") )
  #read and save shapefile otherwise
  if(!(species %in% speciesnames.rda)){
  speciesdata<-readOGR(dsn=".",layer=species)
  #save for future loading
  save(speciesdata, file = paste(species,".rda",sep=""))
  }
  
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
 
  #biggest polygon, #now use all polygons
  #areas= proj.xy$shape_Area
  #shape= proj.xy[which.max(areas),]
  #clip= raster::mask(bio1,shape)
  
  clip= raster::mask(bio1,proj.xy)
 # clip= trim(clip)
  
  #--------
  #find max and min across each column
  clipm=  as.matrix(clip)
  edges1= na.omit(t(apply(clipm, MARGIN=2, firstlast)))
  
  #coldest / warmest cell
  if(clim.k==1) {
    vals= edges1[,which.min(colMeans(edges1, na.rm=TRUE))] #cold edge
    Tmin= min(vals,na.rm=T)
    Tmedian.min= median(vals,na.rm=T)
    T5q.min= quantile(vals, 0.05, na.rm=T)
    T10q.min=quantile(vals, 0.10, na.rm=T)
    Tsd.min= sd(vals,na.rm=T)
    Tmad.min= mad(vals,na.rm=T)
    }
   if(clim.k==2){
     vals= edges1[,which.max(colMeans(edges1))] #warm edge
     Tmax= max(vals,na.rm=T)
     Tmedian.max= median(vals,na.rm=T)
     T5q.max= quantile(vals, 0.95, na.rm=T)
     T10q.max=quantile(vals, 0.90, na.rm=T)
     Tsd.max= sd(vals,na.rm=T)
     Tmad.max= mad(vals,na.rm=T)
   }
  NumberGrids<-length(vals[!is.na(vals)])
  #Area<-gArea(proj.xy)
  } #end clim loop
  
  data.frame(species,UpperLat,LowerLat,Tmin, Tmedian.min, T5q.min,T10q.min,Tsd.min,Tmad.min, Tmax, Tmedian.max, T5q.max,T10q.max,Tsd.max,Tmad.max, NumberGrids,stringsAsFactors=F)

} #end TminTmaxFun
#---------------------------

#TminTmaxFun(species.matched[3])

#Run extraction function
#output<-ldply(species.matched[50:340],TminTmaxFun)
output<-ldply(species.matched[1:length(species.matched)],TminTmaxFun)

#write out
wd="C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"
write.csv(output, paste(wd, "OUT\\MammalTminTmax.csv", sep=""), row.names=F)
