library(ggplot2)
library(stringr)

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

## READ DATA
#-----------------------------
#Read physiology data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_all.csv")
#phy=read.csv("MRelevation_allMASTER.csv")

#----------------------------
#read Msum data
setwd(paste(mydir,"Data\\Msum\\", sep=""))
birdMsum=read.csv("StagerBirdMsum.csv")
mammMsum=read.csv("LovegroveMammMsum.csv")
birdMsum2=read.csv("SwansonBirdMsum.csv")
mammMsum2=read.csv("RezendeMammMsum.csv")
#----------------------------
#add Msum data
phy$Msum_mlO2_h=NA

#add birds
match1= match(phy$Spec.syn, birdMsum$Species) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[matched]<- birdMsum$Msum_mlO2_h[match1[matched]]

birdMsum2$Species= str_trim( as.character(birdMsum2$Species),side="both" )
match1= match(phy$Spec.syn, birdMsum2$Species) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[matched]<- birdMsum2$Msum_mlO2_h[match1[matched]]

#add mammals
mammMsum$Species= paste(mammMsum$gen, mammMsum$species, sep=" ")

match1= match(phy$Spec.syn, mammMsum$Species) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[matched]<- mammMsum$Msum_mlO2_h[match1[matched]]

match1= match(phy$Spec.syn, mammMsum2$Species) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[matched]<- mammMsum2$Msum_mlO2_h[match1[matched]]

#----------------------------
#Proportion Msum at range edge
Tmin= phy$T10q.min
Tmax= phy$T10q.max

NBMR= abs(phy$Tlc- Tmin)*phy$Cmin +phy$BMR_mlO2_h
phy$pMsum<- phy$Msum_mlO2_h / NBMR

NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
phy$pMsum.hot<- phy$Msum_mlO2_h / NBMR

#--------------
plot(density(na.omit(phy$pMsum)))

phy2= phy[which(!is.na(phy$pMsum)),]

#Write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy, "MRelevation_Msum.csv")
