#Plot range shifts

library(ggplot2)
#for mapping
library(ggmap)
library(maps)
library(mapdata)
library(colorRamps)     # for matlab.like(...)
#--------------------
#pick model

mod="CC"
#mod="HD"

#--------------------
mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#rad data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))

rlim.b.cc=read.csv("BirdRlim_CC.csv")
rlim.b.hd=read.csv("BirdRlim_HD.csv")
rlim.m.cc=read.csv("MammalRlim_CC.csv")
rlim.m.hd=read.csv("MammalRlim_HD.csv")

#add field for projection
rlim.b.cc$proj="CC"
rlim.b.hd$proj="HD"
rlim.m.cc$proj="CC"
rlim.m.hd$proj="HD"

rlim=rbind(rlim.b.cc[,1:16], rlim.b.hd[,1:16], rlim.m.cc[,1:16], rlim.m.hd[,1:16])

#match to phy
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRexpansibility_Buckleyetal.csv")

#subset to model
rlim_mod= rlim[rlim$proj==mod,]

#Match species
match1= match(as.character(rlim_mod$phy.Species), as.character(phy$Species) )
pall= cbind(phy[match1,],rlim_mod)

#-----------------------------
#fix infinite medians

pall$pmin.median[!is.finite(pall$pmin.median)]=NA
pall$pmax.median[!is.finite(pall$pmax.median)]=NA
pall$fmin.median[!is.finite(pall$fmin.median)]=NA
pall$fmax.median[!is.finite(pall$fmax.median)]=NA

#===================================================
#Calculate conductance (CMIN)
pall$Cmin= abs((0-pall$BMR_mlO2_h)/(pall$Tb-pall$Tlc) )

#assign temperature metric
Tmin= pall$Tmedian.min
Tmax= pall$Tmedian.max

#Calculate MR elevation
NBMR= abs(pall$Tlc- Tmin)*pall$Cmin +pall$BMR_mlO2_h
pall$MetElev<- NBMR / pall$BMR_mlO2_h

NBMR= abs(pall$Tuc - Tmax)*pall$Cmin +pall$BMR_mlO2_h
pall$MetElev.hot<- NBMR / pall$BMR_mlO2_h

#Add temp prediction
pall$Tamb_lowSS= pall$Tlc-(pall$MetElev-1)* pall$BMR_mlO2_h / pall$Cmin
pall$Tamb_upSS= pall$Tuc+(pall$MetElev.hot-1)* pall$BMR_mlO2_h / pall$Cmin

#===================================================
#Add edges to predict
#Recode North / South range constraints to Cold / Warm boundaries
#C=1: cold constraint; W=1: warm constraint

pall$Cconstrained= NA 
pall$Wconstrained= NA

#Northern hemisphere
inds= which(pall$UpperLat>0 & pall$LowerLat>0  )
pall$Cconstrained[inds]= pall$Nconstrained[inds] 
pall$Wconstrained[inds]= pall$Sconstrained[inds]

#Southern hemisphere
inds= which(pall$UpperLat<0 & pall$LowerLat<0  )
pall$Cconstrained[inds]= pall$Sconstrained[inds] 
pall$Wconstrained[inds]= pall$Nconstrained[inds]

#Crosses hemispheres 
#Change to both cold boundaries, consider constrained if either constraint
inds= which(pall$UpperLat>0 & pall$LowerLat<0  )
con= pall$Sconstrained + pall$Nconstrained
con[which(con==2)]=1

pall$Cconstrained[inds]= con[inds] 
pall$Wconstrained[inds]= 1

#----------
#Specify edges to predict
pall$predC=0
pall$predW=0
pall$predC[which( !is.na(pall$Tamb_lowSS) & pall$Cconstrained==0 )]=1
pall$predW[which( !is.na(pall$Tamb_upSS) & pall$Wconstrained==0 )]=1

#add hemisphere
pall$hemi=NA
inds= which(pall$LowerLat>0)
pall$hemi[inds]="N"
inds= which(pall$UpperLat<0)
pall$hemi[inds]="S"
inds= which(pall$LowerLat<0 & pall$UpperLat>0)
pall$hemi[inds]="B"

#Limit to species with data
pall= pall[pall$predC==1 | pall$predW==1,]

#=======================

#Control for edges to predict

#compare cold and warm edges for unconstrained
pall$cp= NA
pall$cp.med= NA
pall$wp= NA
pall$wp.med= NA

pall$cf= NA
pall$cf.med= NA
pall$wf= NA
pall$wf.med= NA

#North
inds= which( pall$predC==1 & pall$hemi=="N")
pall$cp[inds]= pall$pmax[inds]
pall$cp.med[inds]= pall$pmax.median[inds]
pall$cf[inds]= pall$fmax[inds]
pall$cf.med[inds]= pall$fmax.median[inds]

inds= which( pall$predW==1 & pall$hemi=="N")
pall$wp[inds]= pall$pmin[inds]
pall$wp.med[inds]= pall$pmin.median[inds]
pall$wf[inds]= pall$fmin[inds]
pall$wf.med[inds]= pall$fmin.median[inds]

#South
inds= which( pall$predC==1 & pall$hemi=="S")
pall$cp[inds]= pall$pmin[inds]
pall$cp.med[inds]= pall$pmin.median[inds]
pall$cf[inds]= pall$fmin[inds]
pall$cf.med[inds]= pall$fmin.median[inds]

inds= which( pall$predW==1 & pall$hemi=="S")
pall$wp[inds]= pall$pmax[inds]
pall$wp.med[inds]= pall$pmax.median[inds]
pall$wf[inds]= pall$fmax[inds]
pall$wf.med[inds]= pall$fmax.median[inds]

#Both
inds= which( pall$predC==1 & pall$hemi=="B")
imin= which(abs(pall$UpperLat[inds])>abs(pall$LowerLat[inds])  )  
pall$cp[inds[imin]]= pall$pmax[inds[imin]]
pall$cp.med[inds[imin]]= pall$pmax.median[inds[imin]]
pall$cf[inds[imin]]= pall$fmax[inds[imin]]
pall$cf.med[inds[imin]]= pall$fmax.median[inds[imin]]

imin= which(abs(pall$UpperLat[inds])<abs(pall$LowerLat[inds])  )  
pall$cp[inds[imin]]= pall$pmin[inds[imin]]
pall$cp.med[inds[imin]]= pall$pmin.median[inds[imin]]
pall$cf[inds[imin]]= pall$fmin[inds[imin]]
pall$cf.med[inds[imin]]= pall$fmin.median[inds[imin]]

#--------------
#RANGE SHIFT METRICS 
pall$c.shift= pall$cf - pall$cp 
pall$cmed.shift= pall$cf.med - pall$cp.med 
pall$w.shift= pall$wf - pall$wp 
pall$wmed.shift= pall$wf.med - pall$wp.med 

#plot of median shift
dl=ggplot(pall, aes(cmed.shift, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Latitude shift at cold range boundary (°)")+ scale_fill_manual(values = c("darkgreen","blue"))+theme_bw()

#hist(pall$cmed.shift, breaks=40)

#-----------------
#RANGE SIZE CHANGE
rall= pall[which(pall$predC==1 & pall$predW==1),]

#area
rall$area.shift= rall$farea/rall$parea

da=ggplot(rall, aes(area.shift, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Latitude shift at cold range boundary (°)")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(0,4))

plot(da)

#-----------
#Stats
mshift=(pall[which(pall$Taxa=="Mammal"),"cmed.shift"])
bshift=(pall[which(pall$Taxa=="Bird"),"cmed.shift"])

mshift[!is.finite(mshift)]=NA
bshift[!is.finite(bshift)]=NA

mean(mshift,na.rm=TRUE)
mean(bshift,na.rm=TRUE)
median(mshift,na.rm=TRUE)
median(bshift,na.rm=TRUE)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#MAMMALS
#########################################

library(dismo)
library(raster)
library(maptools)
data(wrld_simpl)

#-------------
#pick model
#mod="HD"
mod="CC"

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
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_all.csv")
#phy=na.omit(phy[,1:30])

#Calculate ambient prediction
#Calculate MR elevation

#TRY MEDIAN
Tmin= phy$Tmedian.min
Tmax= phy$Tmedian.max
#Tmin= phy$T10q.min
#Tmax= phy$T10q.max

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

#----------
#Specify edges to predict
phy$predC=0
phy$predW=0
phy$predC[which( !is.na(phy$Tamb_lowSS) & phy$Cconstrained==0 )]=1
phy$predW[which( !is.na(phy$Tamb_upSS) & phy$Wconstrained==0 )]=1

#add hemisphere
phy$hemi=NA
inds= which(phy$LowerLat>0)
phy$hemi[inds]="N"
inds= which(phy$UpperLat<0)
phy$hemi[inds]="S"
inds= which(phy$LowerLat<0 & phy$UpperLat>0)
phy$hemi[inds]="B"

#Limit to species with data
phy= phy[phy$predC==1 | phy$predW==1,]

#-----------------------------------------
#Estimate range limits in current and future environments
library(gridBase)

##FIG 4
setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf(paste("Fig4",mod,".pdf", sep=""), height=8, width=8)
par(mfrow=c(2,2), cex=1.1, mar=c(3, 3, 1, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

#change directory back for shapefiles
setwd(paste(mydir,"Data\\Shapefiles\\TERRESTRIAL_MAMMALS\\", sep=""))

#speci= rev(which(phy$Species %in% c("Spermophilus beecheyi","Microtus montanus","Myodes gapperi") ))
speci= rev(which(phy$Species %in% c("Tamiasciurus hudsonicus", "Neotoma lepida","Peromyscus californicus" ) ))
#Tamiasciurus hudsonicus, Red Squirrel 
#Neotoma lepida, desert woodrat
#Peromyscus californicus, California Mouse

#LOOP SPECIES
for(spec in speci ){
  
  #LOAD SHAPEFILE AND EXTRACT EXT
  shape= readRDS(paste(phy[spec,"Spec.syn"],".rds",sep=""))
  
  extent2= extent(shape)
  
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
    
    if(!is.na(sum(as.matrix(clim.spec.crop))) )clim.spec.crop= trim(clim.spec.crop)
    
    if(clim==1) range.p= clim.spec.crop
    if(clim==2) range.f= clim.spec.crop
    
  } #end clim loop
 
  #-----------------------------
  #RANGE MAPS
  xlims=c(-130,-75)
  ylims=c(25,80)
  
  image(range.p, col="blue", main=phy[spec,"Species"], xlim=xlims, ylim=ylims, ylab="Latitude (°)", xlab="Longitude (°)", cex=1.2, cex.lab=1.2)
  plot(wrld_simpl, add=TRUE, border="gray")
  plot(range.f, col=rainbow(1, alpha=0.5),add=TRUE, legend=FALSE)
  plot(shape, add=TRUE)
 
  if(spec==169) legend("topleft", legend= c("1960-1990","2070") , fill=c("blue", rainbow(1, alpha=0.5)), bty="n")
  
} #end species loop  
  #=======================================
  
#add ggplot
 
  plot.new()              
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp1 <-plotViewport(c(0,0,0,0))  #(c(1.8,1,0,1)) ## create new vp with margin
  print(dl,vp = vp1)        
  
  dev.off()
  