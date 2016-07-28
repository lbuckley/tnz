 ## set directory for data files

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#READ PHYS DATA
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy= read.csv("Phy_all.csv")

#---------------------------------------

#ADD Tmin / Tmax
phy$Tmin= NA; phy$Tmedian.min= NA; phy$T5q.min= NA; phy$T10q.min= NA; phy$Tsd.min= NA;
phy$Tmax= NA; phy$Tmedian.max= NA; phy$T5q.max= NA; phy$T10q.max= NA; phy$Tsd.max= NA;

phy$UpperLat= NA; phy$LowerLat=NA; phy$ShapeName=NA; phy$Area=NA;

#MAMMALS
#Add Tmin Tmix data
setwd(paste(mydir,"Out\\", sep=""))
TminTmax= read.csv("MammalTminTmax.csv", na.strings = c("NA","Inf","-Inf"))

#Match species
match1= match(as.character(phy$Spec.syn), as.character(TminTmax$species) )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy$Tmin[matched]<- TminTmax$Tmin[match1[matched]]
phy$Tmedian.min[matched]<- TminTmax$Tmedian.min[match1[matched]]
phy$T5q.min[matched]<- TminTmax$T5q.min[match1[matched]]
phy$T10q.min[matched]<- TminTmax$T10q.min[match1[matched]]
phy$Tsd.min[matched]<- TminTmax$Tsd.min[match1[matched]]

phy$Tmax[matched]<- TminTmax$Tmax[match1[matched]]
phy$Tmedian.max[matched]<- TminTmax$Tmedian.max[match1[matched]]
phy$T5q.max[matched]<- TminTmax$T5q.max[match1[matched]]
phy$T10q.max[matched]<- TminTmax$T10q.max[match1[matched]]
phy$Tsd.max[matched]<- TminTmax$Tsd.max[match1[matched]]

phy$UpperLat[matched]<- TminTmax$UpperLat[match1[matched]]
phy$LowerLat[matched]<- TminTmax$LowerLat[match1[matched]]
phy$Area[matched]<- TminTmax$NumberGrids[match1[matched]]

#BIRDS
#Add Tmin Tmix data
setwd(paste(mydir,"Out\\", sep=""))
TminTmax= read.csv("BirdTminTmax.csv", na.strings = c("NA","Inf","-Inf"))

specgen<- gsub("_", " ", TminTmax$species)
TminTmax$genspec<- substr(specgen, 1, nchar(specgen)-9)

#Match species
match1= match(as.character(phy$Spec.syn), as.character(TminTmax$genspec) )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy$Tmin[matched]<- TminTmax$Tmin[match1[matched]]
phy$Tmedian.min[matched]<- TminTmax$Tmedian.min[match1[matched]]
phy$T5q.min[matched]<- TminTmax$T5q.min[match1[matched]]
phy$T10q.min[matched]<- TminTmax$T10q.min[match1[matched]]
phy$Tsd.min[matched]<- TminTmax$Tsd.min[match1[matched]]

phy$Tmax[matched]<- TminTmax$Tmax[match1[matched]]
phy$Tmedian.max[matched]<- TminTmax$Tmedian.max[match1[matched]]
phy$T5q.max[matched]<- TminTmax$T5q.max[match1[matched]]
phy$T10q.max[matched]<- TminTmax$T10q.max[match1[matched]]
phy$Tsd.max[matched]<- TminTmax$Tsd.max[match1[matched]]

phy$UpperLat[matched]<- TminTmax$UpperLat[match1[matched]]
phy$LowerLat[matched]<- TminTmax$LowerLat[match1[matched]]
phy$Area[matched]<- TminTmax$NumberGrids[match1[matched]]
phy$ShapeName[matched]<- as.character(TminTmax$shapename[match1[matched]])

#Calculate ambient prediction
#Calculate MR elevation

Tmin= phy$T10q.min
Tmax= phy$T10q.max
#Tmin= phy$T5q.min
#Tmax= phy$T5q.max
#================================
#METABOLIC CALCULATIONS

#Calculate CMIN at upper TNZ
phy$Cmin= abs((0-phy$BMR_mlO2_h)/(phy$Tb-phy$Tlc) )
plot(phy$CMIN_mlO2_hC, phy$Cmin )
abline(a=0,b=1)

#Calculate MR elevation
NBMR= abs(phy$Tlc- Tmin)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev<- NBMR / phy$BMR_mlO2_h

NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev.hot<- NBMR / phy$BMR_mlO2_h

#----------------------
## CALCULATE SPECIES SPECIFIC Tamb

#Add temp prediction
phy$Tamb_lowSS= phy$Tlc-(phy$MetElev-1)* phy$BMR_mlO2_h / phy$Cmin

phy$Tamb_upSS= phy$Tuc+(phy$MetElev.hot-1)* phy$BMR_mlO2_h / phy$Cmin

#========================================
#recode constrained
phy$Cconstrained= NA 
phy$Wconstrained= NA
#N
inds= which(phy$UpperLat>0 & phy$LowerLat>0  )
phy$Cconstrained[inds]= phy$Nconstrained[inds] 
phy$Wconstrained[inds]= phy$Sconstrained[inds]

#S
inds= which(phy$UpperLat<0 & phy$LowerLat<0  )
phy$Cconstrained[inds]= phy$Sconstrained[inds] 
phy$Wconstrained[inds]= phy$Nconstrained[inds]

#crosses latitude ## BETTER WAY?
inds= which(phy$UpperLat>0 & phy$LowerLat<0  )
con= phy$Sconstrained + phy$Nconstrained
con[which(con==2)]=1

phy$Cconstrained[inds]= con[inds] 
phy$Wconstrained[inds]= 1
#-------

phy= phy[-which(is.na(phy$Species)),]

#write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy, "MRelevation_all.csv", row.names = FALSE)
#--------------------------------------
