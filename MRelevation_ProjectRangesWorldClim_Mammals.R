#########################################
## TNZ analysis
# Project ranges in current and future climates
# Estimaates ranges based on physiological limits
# Don't consider west / east range boundaries, use observed limits?
#
#MAMMALS
#########################################

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
clim.f=getData('CMIP5', var='bio', res=5, rcp=60, model='BC', year=70)
clim.fmax= clim.f$bc60bi705/10
clim.fmin= clim.f$bc60bi706/10
#'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO".
#'rcp' should be one of 26, 45, 60, or 85.
#'year' should be 50 or 70

plot(clim.fmax-clim.pmax)
#-----------------------------------------
#Load physiological data
#setwd(paste(mydir,"MRelevation\\Data\\", sep=""))
#phy=read.csv("MRelevation_wTraits.csv")
#phy=na.omit(phy[,1:30])

#Read physiology data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_allMASTER.csv")
#phy=na.omit(phy[,c(5:32,53:54)])

#Calculate ambient prediction
#Calculate MR elevation

#TRY MEDIAN
Tmin= phy$T10q.min
Tmax= phy$T10q.max

#Calculate MR elevation
NBMR= abs(phy$Tlc- Tmin)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev<- NBMR / phy$BMR_mlO2_h

NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev.hot<- NBMR / phy$BMR_mlO2_h

#Add temp prediction
phy$Tamb_lowSS= phy$Tlc-(phy$MetElev-1)* phy$BMR_mlO2_h / phy$Cmin
phy$Tamb_upSS= phy$Tuc+(phy$MetElev.hot-1)* phy$BMR_mlO2_h / phy$Cmin

#Limit to mammals
phy$Taxa=gsub("Mammals","Mammal",phy$Taxa)
phy= phy[which(phy$Taxa=="Mammal"),]

#--------
#recode constrained
phy$Cconstrained= NA 
phy$Wconstrained= NA
phy$hemi=NA
#N
inds= which(phy$UpperLat>0 & phy$LowerLat>0  )
phy$Cconstrained[inds]= phy$Nconstrained[inds] 
phy$Wconstrained[inds]= phy$Sconstrained[inds]
phy$hemi[inds]= "N"

#S
inds= which(phy$UpperLat<0 & phy$LowerLat<0  )
phy$Cconstrained[inds]= phy$Sconstrained[inds] 
phy$Wconstrained[inds]= phy$Nconstrained[inds]
phy$hemi[inds]= "S"

#crosses latitude ## BETTER WAY?
inds= which(phy$UpperLat>0 & phy$LowerLat<0  )
con= phy$Sconstrained + phy$Nconstrained
con[which(con==2)]=1

phy$Cconstrained[inds]= con[inds] 
phy$Wconstrained[inds]= 1
phy$hemi[inds]= "B"

#Specify edges to predict
phy$predC=0
phy$predW=0
phy$predC[which( !is.na(phy$Tamb_lowSS) & phy$Cconstrained==0 )]=1
phy$predW[which( !is.na(phy$Tamb_upSS) & phy$Wconstrained==0 )]=1

#Limit to species with data
phy= phy[phy$predC==1 | phy$predW==1,]

#-----------------------------------------
#Estimate range limits in current and future environments

rlim= matrix(NA, nrow=nrow(phy), ncol=13) #ymin past, ymax past, ymin future, ymax future, area past, area future, ymin med past, ymax med past, ymin med future, ymax med future

#plot range maps
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
pdf("MammalRangePred.pdf", height = 6, width = 10)

par(mfrow=c(4,4), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

#change directory back for shapefiles
setwd(paste(mydir,"Data\\Shapefiles\\TERRESTRIAL_MAMMALS\\", sep=""))

#LOOP SPECIES
for(spec in 1:nrow(phy) ){

	#LOAD SHAPEFILE AND EXTRACT EXT
	shape= shapefile(paste(phy[spec,"Spec.syn"],".shp",sep=""))  
	extent2= extent(shape)
  #load(file = paste(phy[spec,"Spec.syn"],".rda",sep="") )
	
for(clim in 1:2){ #present, future  

  if(clim==1) clim.pmin1= clim.pmin; clim.pmax1= clim.pmax
  if(clim==2) clim.pmin1= clim.fmin; clim.pmax1= clim.fmax

## CURRENT CLIMATE
#Subset to above Tmin and below Tmax
#account for species with only one constraint
if(phy$predC[spec]==1) clim.pmin1[clim.pmin1< phy$Tamb_lowSS[spec]]<- NA
if(phy$predW[spec]==1) clim.pmax1[clim.pmax1> phy$Tamb_upSS[spec]]<- NA

#specify prediction
if(phy$predC[spec]==1 & phy$predW[spec]==1) clim.spec= clim.pmin1 + clim.pmax1
if(phy$predC[spec]==1 & phy$predW[spec]==0) clim.spec= clim.pmin1
if(phy$predC[spec]==0 & phy$predW[spec]==1) clim.spec= clim.pmax1
  
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
df <- data.frame(xy, clump_id, is_clump = sc[] %in% freq(sc, useNA = 'no')[,1])
df= df[df$is_clump == T, ]

clumps= unique(df[df$y<extent2@ymax & df$y>extent2@ymin,"clump_id"])
rclump= df[df$"clump_id" %in% clumps, ]

clim.spec.crop[!sc %in% clumps] <- NA
#-----------------------------

# calculate frequency of each clump/patch
sc.freq <- as.data.frame(freq(sc))
sc.freq<- sc.freq[!is.na(sc.freq$value),]

#get rid of small clumps
## store clump ID's with small size
rmID <- sc.freq$value[which(sc.freq$count > sum(sc.freq$count)/100)]
##restrict to biggest clump
#rmID <- sc.freq$value[which.max(sc.freq$count)]
clim.spec.crop[!sc %in% rmID] <- NA
clim.spec.crop[sc %in% rmID] <- 1 

clim.spec.crop= trim(clim.spec.crop)

##Cut off tails, need to account for Central America, etc
#rsum= rowSums(as.matrix(clim.spec.crop), na.rm=TRUE) 
#clim.spec.crop[ which(rsum< (max(rsum)/10)), ]<- NA #cut off tails<10% row sums

#-----------------------------
##restrict to biggest clump again
###find clumps
#sc= clump(clim.spec.crop, directions=4)

## calculate frequency of each clump/patch
#sc.freq <- as.data.frame(freq(sc))
#sc.freq<- sc.freq[!is.na(sc.freq$value),]

## store clump ID's with small size
##rmID <- sc.freq$value[which(sc.freq$count <100)]
##restrict to biggest clump
#rmID <- sc.freq$value[which.max(sc.freq$count)]
#clim.spec.crop[!sc %in% rmID] <- NA
#clim.spec.crop[sc %in% rmID] <- 1 

range.s= extent(clim.spec.crop)

if(clim==1) range.p= clim.spec.crop
if(clim==2) range.f= clim.spec.crop

#Control for constraints
lat.range= c(range.s@ymin, range.s@ymax) 
if(phy$predC[spec]==0 & phy$hemi[spec]=="N")lat.range[2]=NA
if(phy$predC[spec]==0 & phy$hemi[spec]=="S")lat.range[1]=NA
if(phy$predC[spec]==0 & phy$hemi[spec]=="B"){lat.range[1]=NA; lat.range[2]=NA}
if(phy$predW[spec]==0 & phy$hemi[spec]=="N")lat.range[1]=NA
if(phy$predW[spec]==0 & phy$hemi[spec]=="S")lat.range[2]=NA
  
#find range
if(clim==1) rlim[spec,1:2]= lat.range 
if(clim==2) rlim[spec,3:4]= lat.range

#find lat range by column
yFromRow(clim.spec.crop, row=1:nrow(clim.spec.crop))

pred1=clim.spec.crop
pred1[!is.na(pred1)] <- 1

yr= raster(clim.spec.crop)
yr <- init(yr, 'y')
yr= yr*pred1
yr= as.matrix(yr)

ymax=  median(apply(yr, MARGIN=2, max, na.rm=TRUE),na.rm=TRUE)
ymin=  median(apply(yr, MARGIN=2, min, na.rm=TRUE),na.rm=TRUE)

if(clim==1) rlim[spec,10:11]= c(ymin,ymax) 
if(clim==2) rlim[spec,12:13]= c(ymin,ymax)

} #end clim loop
#----------------
#plot

#set extent
extent1= union(extent(range.p),extent(range.f))
extent2= extent(shape)
cext= union(extent1, extent2)

image(range.p, col="blue", main=phy[spec,"Species"], xlim=c(cext@xmin,cext@xmax), ylim=c(cext@ymin,cext@ymax))
plot(wrld_simpl, add=TRUE, border="gray")
#plot(range.p, col="blue", legend=FALSE, main=phy[spec,"Species"], xlim=c(cext@xmin,cext@xmax), ylim=c(cext@ymin,cext@ymax))
#plot(range.p, col="blue", legend=FALSE, main=phy[spec,"Species"], ext=cext)
plot(range.f, col=rainbow(1, alpha=0.5),add=TRUE, legend=FALSE)
plot(shape, add=TRUE)

#-------------------
#Extract stats
rlim[spec,5]=cellStats(range.p, stat='sum', na.rm=TRUE)
rlim[spec,6]=cellStats(range.f, stat='sum', na.rm=TRUE)

range.p2=range.p*-2

#cell stats to get area change N/S, use sum -2,1?
range.b=mosaic(range.p2, range.f, fun=mean, tolerance=0.05)
freqs= raster::freq(range.b, useNA='no')

freq1= freqs[freqs[,"value"]==-2, "count"]
if(length(freq1>0)) rlim[spec,7]= freq1

freq1= freqs[freqs[,"value"]==-1, "count"]
if(length(freq1>0)) rlim[spec,8]= freq1

freq1= freqs[freqs[,"value"]==1, "count"]
if(length(freq1>0)) rlim[spec,9]= freq1

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
colnames(rlim)= c('pmin','pmax','fmin','fmax', 'parea','farea','ponly','pandf','fonly')
rlim$dmax= rlim$fmax - rlim$pmax
rlim$dmin= rlim$fmin - rlim$pmin
rlim$index= 1:nrow(rlim)

#write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
rlim.out= cbind(phy$Species, rlim)
write.csv(rlim.out,"MammalRlim.csv")

#convert to poleward / equatorial shift
#! need to fix for species spanning equator
rlim$poleward= rlim$dmax
rlim$poleward[which(rlim$pmax<0 & rlim$pmin<0)]= abs(rlim$dmin[which(rlim$pmax<0 & rlim$pmin<0)])

rlim$eq= abs(rlim$dmin)
rlim$eq[which(rlim$pmax<0 & rlim$pmin<0)]= rlim$dmax[which(rlim$pmax<0 & rlim$pmin<0)]

#plot(rlim$index, rlim$fmax, col="red")
#points(rlim$index, rlim$pmax, col="blue")

move= na.omit(rlim$poleward)
move= move[is.finite(move) ]
hist(move, breaks=50)

move= na.omit(rlim$eq)
move= move[is.finite(move) ]
hist(move, breaks=20)

#area
d.area= rlim$farea/rlim$parea
hist(d.area, breaks=20)

hist(rlim$ponly/rlim$parea)
hist(rlim$fonly/rlim$parea)
# 'ponly','pandf','fonly'