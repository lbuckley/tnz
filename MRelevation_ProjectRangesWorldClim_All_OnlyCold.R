#########################################
## TNZ analysis
# Project ranges in current and future climates
# Estimaates ranges based on physiological limits
# Don't consider west / east range boundaries, use observed limits?
#
#########################################
#pick model
mod="HD"
#mod="CC"

#pick taxa
#tax="Bird"
tax="Mammal"

library(dismo)
library(raster)
library(maptools)
data(wrld_simpl)

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#USE WORLD CLIM PROJECTIONS
#5 Minutes
# BIO6 = Min Temperature of Coldest Month
# 2100 RCP6

#current data
clim.p= getData('worldclim', var='bio', res=5)
clim.pmax= clim.p$bio5/10
clim.pmin= clim.p$bio6/10

#future data
if(mod=="HD") {
  clim.f=getData('CMIP5', var='bio', res=5, rcp=60, model='HD', year=70) #HadGEM2-AO
  clim.fmax= clim.f$hd60bi705/10
  clim.fmin= clim.f$hd60bi706/10
}
  
if(mod=="CC") {
  clim.f=getData('CMIP5', var='bio', res=5, rcp=60, model='CC', year=70) #CCSM4
  clim.fmax= clim.f$cc60bi705/10
  clim.fmin= clim.f$cc60bi706/10
}  
  
#'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO".
#'rcp' should be one of 26, 45, 60, or 85.
#'year' should be 50 or 70
# models: http://worldclim.org/CMIP5_30s

plot(clim.fmax-clim.pmax)
#-----------------------------------------
#Load physiological data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRexpansibility_Buckleyetal.csv")

#Calculate ambient prediction
#Calculate MR elevation

#Calculate conductance (CMIN)
phy$Cmin= abs((0-phy$BMR_mlO2_h)/(phy$Tb-phy$Tlc) )

#assign temperature metric
Tmin= phy$Tmedian.min
Tmax= phy$Tmedian.max

#Calculate MR elevation
NBMR= abs(phy$Tlc- Tmin)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev<- NBMR / phy$BMR_mlO2_h

NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev.hot<- NBMR / phy$BMR_mlO2_h

#Add temp prediction
phy$Tamb_lowSS= phy$Tlc-(phy$MetElev-1)* phy$BMR_mlO2_h / phy$Cmin
phy$Tamb_upSS= phy$Tuc+(phy$MetElev.hot-1)* phy$BMR_mlO2_h / phy$Cmin

#Limit to taxa
if(tax=="Bird") phy= phy[which(phy$Taxa=="Bird"),]
if(tax=="Mammal") phy= phy[which(phy$Taxa=="Mammal"),]

#======================================
#Recode North / South range constraints to Cold / Warm boundaries
#C=1: cold constraint; W=1: warm constraint

phy$Cconstrained= NA 

#Northern hemisphere
inds= which(phy$UpperLat>0 & phy$LowerLat>0  )
phy$Cconstrained[inds]= phy$Nconstrained[inds] 

#Southern hemisphere
inds= which(phy$UpperLat<0 & phy$LowerLat<0  )
phy$Cconstrained[inds]= phy$Sconstrained[inds] 

#Crosses hemispheres 
#Change to both cold boundaries, consider constrained if either constraint
inds= which(phy$UpperLat>0 & phy$LowerLat<0  )
con= phy$Sconstrained + phy$Nconstrained
con[which(con==2)]=1

phy$Cconstrained[inds]= con[inds] 

#drop cold constrained ranges
phy=phy[which(phy$Cconstrained==0),]

#----------

#add hemisphere
phy$hemi=NA
inds= which(phy$LowerLat>0)
phy$hemi[inds]="N"
inds= which(phy$UpperLat<0)
phy$hemi[inds]="S"
inds= which(phy$LowerLat<0 & phy$UpperLat>0)
phy$hemi[inds]="B"

#sort by species name
phy= phy[order(phy$Species),]

#-----------------------------------------
#Estimate range limits in current and future environments

rlim= matrix(NA, nrow=nrow(phy), ncol=13) #ymin past, ymax past, ymin future, ymax future, area past, area future, ymin med past, ymax med past, ymin med future, ymax med future

#plot range maps
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
pdf(paste(tax,"RangePred_",mod,".pdf", sep=""), height = 10, width = 10)

par(mfrow=c(4,4), cex=1.2, mar=c(3, 3, 1.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

#change directory back for shapefiles
if(tax=="Bird") setwd(paste(mydir,"Data\\Shapefiles\\Birds\\", sep=""))
if(tax=="Mammal") setwd(paste(mydir,"Data\\Shapefiles\\TERRESTRIAL_MAMMALS\\", sep=""))

#get list of saved shapefiles 
if(tax=="Bird") sdir= paste(mydir,"Data\\Shapefiles\\Birds\\", sep="")
if(tax=="Mammal") sdir= paste(mydir,"Data\\Shapefiles\\TERRESTRIAL_MAMMALS\\", sep="")
  
speciesfiles.rds<-list.files(sdir,pattern="\\.rds$")#get shape files
speciesnames.rds<- gsub("\\.rds", "", speciesfiles.rds)

#---------
#LOOP SPECIES
for(spec in 1:nrow(phy) ){

  #LOAD SHAPEFILE AND EXTRACT EXT
  #shape= shapefile(paste(phy[spec,"ShapeName"],".shp",sep="")); saveRDS(shape, paste(phy[spec,"ShapeName"],".rds",sep=""))
  if(tax=="Bird") shape= readRDS(paste(phy[spec,"ShapeName"],".rds",sep=""))
  if(tax=="Mammal") shape= readRDS(paste(phy[spec,"Spec.syn"],".rds",sep=""))
  
  extent2= extent(shape)
  
for(clim in 1:2){ #present, future  

  if(clim==1) clim.pmin1= clim.pmin; clim.pmax1= clim.pmax
  if(clim==2) clim.pmin1= clim.fmin; clim.pmax1= clim.fmax
    
  ## CURRENT CLIMATE
  #Subset to above Tmin
  clim.pmin1[clim.pmin1< phy$Tamb_lowSS[spec]]<- NA

  #specify prediction
  clim.spec= clim.pmin1
  clim.spec= trim(clim.spec)

#crop to lon extent
ext.spec= extent(clim.spec)
ext= extent(t(matrix(c( sort(c(extent2@xmin, extent2@xmax), decreasing=FALSE), ext.spec@ymin, ext.spec@ymax), nrow=2)))
#account for no longitudinal range
if(ext@xmin==ext@xmax) ext@xmin=ext@xmin-1; ext@xmax=ext@xmax+1
clim.spec.crop=crop(clim.spec, ext)

if(cellStats(clim.spec.crop,sum, ra.rm=TRUE)==0) break()
  
clim.spec.crop= trim(clim.spec.crop)

#-----------------------------
## find clumps overlapping initial distribution
##find clumps
sc= clump(clim.spec.crop, directions=4)

clump_id <- getValues(sc) 
xy <- xyFromCell(sc,1:ncell(sc))

#clump freqs
freqs= as.data.frame( freq(sc, useNA = 'no') )
#remove small clumps
freqs <- freqs[which(freqs$count> max(freqs$count)/20 ),]
#find clumps overlapping distribution  
df <- data.frame(xy, clump_id, is_clump = sc[] %in% freqs$value)
df= df[df$is_clump == T, ]
df1= aggregate(df, list(df$clump_id),max)
df2=  aggregate(df, list(df$clump_id),min)
df1$min= df2$y

df1= df1[which(df1$y>extent2@ymin & df1$min<extent2@ymax),]

clim.spec.crop[!(sc %in% df1$clump_id)] <- NA

if(!is.na(sum(as.matrix(clim.spec.crop))) ) clim.spec.crop= trim(clim.spec.crop)

#-----------------------------

range.s= extent(clim.spec.crop)

#Control for constraints
lat.range= c(range.s@ymin, range.s@ymax) 

if(phy$hemi[spec]=="N")lat.range[1]=NA
if(phy$hemi[spec]=="S")lat.range[2]=NA

#find range
if(clim==1) rlim[spec,1:2]= lat.range #past lat range
if(clim==2) rlim[spec,3:4]= lat.range #future lat range

#find lat range by column
pred1=clim.spec.crop
pred1[!is.na(pred1)] <- 1

#record area
if(clim==1) range.p= pred1
if(clim==2) range.f= pred1

yr= raster(clim.spec.crop)
yr <- init(yr, 'y')
yr= yr*pred1
yr= as.matrix(yr)

ymax=  median(apply(yr, MARGIN=2, max, na.rm=TRUE),na.rm=TRUE)
ymin=  median(apply(yr, MARGIN=2, min, na.rm=TRUE),na.rm=TRUE)

#Control for constraints
if(phy$hemi[spec]=="N") ymin=NA
if(phy$hemi[spec]=="S") ymax=NA

if(clim==1) rlim[spec,10:11]= c(ymin,ymax) # past range medians
if(clim==2) rlim[spec,12:13]= c(ymin,ymax) # future range medians

} #end clim loop
#----------------
#plot

#set extent
ymin= extent(shape)@ymin
  
extent1= union(extent(range.p),extent(range.f))
extent2= extent(shape)
cext= union(extent1, extent2)

#set warm extent to shape file
#N hemi
if(phy$UpperLat[spec]>0 & phy$LowerLat[spec]>0) cext@ymin= extent(shape)@ymin -5  
#S hemi
if(phy$UpperLat[spec]<0 & phy$LowerLat[spec]<0) cext@ymax= extent(shape)@ymax +5
#both hemispheres
if(phy$UpperLat[spec]>0 & phy$LowerLat[spec]<0 & phy$UpperLat[spec]>abs(phy$LowerLat[spec]) ) cext@ymin= extent(shape)@ymin -5
if(phy$UpperLat[spec]>0 & phy$LowerLat[spec]<0 & phy$UpperLat[spec]<abs(phy$LowerLat[spec]) ) cext@ymax= extent(shape)@ymax +5

#make extent square
xl= cext@xmax - cext@xmin
yl= cext@ymax - cext@ymin
xyl= max(xl, yl)
if(xl<yl){
  mid= mean( c(cext@xmax, cext@xmin) )
  cext@xmin= mid-xyl/2
  cext@xmax= mid+xyl/2
  if(cext@xmin< (-180)) cext@xmin= -180
  if(cext@xmax> 180) cext@xmax= 180
}
if(xl>yl){
  mid= mean(cext@ymax, cext@ymin)
  cext@ymin= mid-xyl/2
  cext@ymax= mid+xyl/2
  if(cext@ymin< (-90)) cext@ymin= -90
  if(cext@ymax> 90) cext@ymax= 90
}

image(range.p, col="blue", main=phy[spec,"Species"], xlim=c(cext@xmin,cext@xmax), ylim=c(cext@ymin,cext@ymax))
plot(wrld_simpl, add=TRUE, border="gray")
#plot(range.p, col="blue", legend=FALSE, main=phy[spec,"Species"], xlim=c(cext@xmin,cext@xmax), ylim=c(cext@ymin,cext@ymax))
#plot(range.p, col="blue", legend=FALSE, main=phy[spec,"Species"], ext=cext)
plot(range.f, col=rainbow(1, alpha=0.5),add=TRUE, legend=FALSE)
plot(shape, add=TRUE)

#-------------------

#Extract stats
rlim[spec,5]=cellStats(range.p, stat='sum', na.rm=TRUE) #past area
rlim[spec,6]=cellStats(range.f, stat='sum', na.rm=TRUE) #future area

range.p2=range.p*-2

#cell stats to get area change N/S, use sum -2,1?
range.b=mosaic(range.p2, range.f, fun=mean, tolerance=0.05)
freqs= raster::freq(range.b, useNA='no')

freq1= freqs[freqs[,"value"]==-2, "count"]
if(length(freq1>0)) rlim[spec,7]= freq1 #past only

freq1= freqs[freqs[,"value"]==-1, "count"]
if(length(freq1>0)) rlim[spec,8]= freq1 #both fast and future

freq1= freqs[freqs[,"value"]==1, "count"]
if(length(freq1>0)) rlim[spec,9]= freq1 #future only

#-------
##range extent
#xy <- xyFromCell(clim.spec.crop,1:ncell(clim.spec.crop))
#if(clim==1) rlim[spec,1:2]= range(xy[,"y"])
#if(clim==2) rlim[spec,3:4]= range(xy[,"y"])

  print(spec)
} #end spec loop

dev.off() #end mapping

#======================
# PLOT
rlim= as.data.frame(rlim)
colnames(rlim)= c('pmin','pmax','fmin','fmax', 'parea','farea','ponly','pandf','fonly', 'pmin.median','pmax.median','fmin.median','fmax.median')

#write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
rlim.out= cbind(phy$Species, rlim)
write.csv(rlim.out, paste(tax,"Rlim_",mod,".csv", sep=""))
#=======================
