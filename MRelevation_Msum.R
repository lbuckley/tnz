library(ggplot2)
library(stringr)
library(grid)

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

## READ DATA
#-----------------------------
#Read physiology data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
#phy=read.csv("MRexpansibility_Buckleyetal.csv")
#phy=read.csv("MRelevation_all.csv")

#include constrained species
#phy=read.csv("MRelevation_all_wConstrained.csv")
phy=read.csv("MRexpansibility_Buckleyetal_wQual_noUCTdrop.csv")

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
phy$BMR_mlO2_h_Msum=phy$BMR_mlO2_h

#add birds
match1= match(birdMsum$Species, phy$Spec.syn) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[match1[matched]]<- birdMsum$Msum_mlO2_h[matched]
phy$BMR_mlO2_h_Msum[match1[matched]]<- birdMsum$BMR_mlO2_h[matched]

birdMsum2$Species= str_trim( as.character(birdMsum2$Species),side="both" )
match1= match(birdMsum2$Species, phy$Spec.syn) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[match1[matched]]<- birdMsum2$Msum_mlO2_h[matched]

#add mammals
mammMsum$Species= paste(mammMsum$gen, mammMsum$species, sep=" ")

match1= match(mammMsum$Species, phy$Spec.syn) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[match1[matched]]<- mammMsum$Msum_mlO2_h[matched]
phy$BMR_mlO2_h_Msum[match1[matched]]<- mammMsum$BMR__mlO2_h[matched]

match1= match(mammMsum2$Species, phy$Spec.syn) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[match1[matched]]<- mammMsum2$Msum_mlO2_h[matched]
phy$BMR_mlO2_h_Msum[match1[matched]]<- mammMsum2$BMR__mlO2_h[matched]

# drop species with Msum~ BMR
phy[which(phy$Species=="Sorex cinereus"), "Msum_mlO2_h"]= NA
phy[which(phy$Species=="Sorex cinereus"), "BMR_mlO2_h_Msum"]= NA

#count of birds and mammals
table(phy2$Taxa)

#----------------------------
#Proportion Msum at range edge
Tmin= phy$Tmedian.min
Tmax= phy$Tmedian.max

NBMR= abs(phy$Tlc- Tmin)*phy$Cmin +phy$BMR_mlO2_h
phy$pMsum<- NBMR / phy$Msum_mlO2_h

NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
phy$pMsum.hot<- NBMR / phy$Msum_mlO2_h 

phy2= phy[which(!is.na(phy$Msum_mlO2_h)),]

#find peak
d.msum= density(phy2$pMsum)
peak.msum=d.msum$x[which.max(d.msum$y)] 

#median=0.88, mean=0.96
median(phy$pMsum, na.rm=TRUE)
mean(phy$pMsum, na.rm=TRUE)
MsumE= peak.msum #peak=0.74

sd(phy$pMsum, na.rm=TRUE)
Mdif= phy$Msum_mlO2_h- phy$BMR_mlO2_h

phy$Tamb_low_Msum= phy$Tlc- MsumE*(phy$Msum_mlO2_h- phy$BMR_mlO2_h_Msum) / phy$Cmin

phymsum= phy[!is.na(phy$pMsum),]

#--------------

#Write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy, "MRelevation_Msum.csv")

write.csv(phymsum, "MRelevation_Msum_nonoa.csv")

#===============================
#PLOT

hist(na.omit(phy$pMsum),breaks=8)

#-----------

hl=ggplot(phy, aes(pMsum, fill = Taxa)) + 
  geom_histogram(binwidth = 0.2)+labs(x=expression(MR[CRB] / Msum))+ scale_fill_manual(values = c("darkgreen","blue"))+theme_bw() +xlim(c(0.25,2.25))+theme(axis.title=element_text(size=rel(1.3)))

#Plot TRAITS

#lower
xyrange= range(c(phy$Tamb_low_Msum, phy$Tmin.use), na.rm=TRUE)
xyrange[1]= -50

phy$Torpor= as.factor(phy$torpor)
phy$Mass= phy$Mass_g

p <- ggplot(data = phy, aes(x = Tamb_low_Msum, y = T10q.min, shape=Taxa, color=Torpor, size= log(Mass))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pl= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.3)))
#+ facet_wrap(~Taxa)

#------------------------
setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("Fig1_Msum.pdf", height=10, width=5)

#plot
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(hl,vp=vplayout(1,1))
print(pl,vp=vplayout(2,1))

dev.off()

#--------------------------
#msum factor
phymsum$Msum_fact=phymsum$Msum_mlO2_h/phymsum$BMR_mlO2_h

plot(phymsum$Msum_fact, phymsum$pMsum)
summary(phymsum$Msum_fact)
