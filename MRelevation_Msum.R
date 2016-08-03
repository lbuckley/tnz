library(ggplot2)
library(stringr)
library(grid)

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
phy$BMR_mlO2_h_Msum=phy$BMR_mlO2_h

#add birds
match1= match(phy$Spec.syn, birdMsum$Species) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[matched]<- birdMsum$Msum_mlO2_h[match1[matched]]
phy$BMR_mlO2_h_Msum[matched]<- birdMsum$BMR_mlO2_h[match1[matched]]

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
phy$BMR_mlO2_h_Msum[matched]<- mammMsum$BMR__mlO2_h[match1[matched]]

match1= match(phy$Spec.syn, mammMsum2$Species) 
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
phy$Msum_mlO2_h[matched]<- mammMsum2$Msum_mlO2_h[match1[matched]]
phy$BMR_mlO2_h_Msum[matched]<- mammMsum2$BMR__mlO2_h[match1[matched]]

#----------------------------
#Proportion Msum at range edge
Tmin= phy$T10q.min
Tmax= phy$T10q.max

NBMR= abs(phy$Tlc- Tmin)*phy$Cmin +phy$BMR_mlO2_h
phy$pMsum<- phy$Msum_mlO2_h / NBMR

NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
phy$pMsum.hot<- phy$Msum_mlO2_h / NBMR

#use peak of 1: median=0.95, mean=1.0
median(phy$pMsum, na.rm=TRUE)
mean(phy$pMsum, na.rm=TRUE)

sd(phy$pMsum, na.rm=TRUE)
Mdif= phy$Msum_mlO2_h- phy$BMR_mlO2_h

phy$Tamb_low_Msum= phy$Tlc- (phy$Msum_mlO2_h- phy$BMR_mlO2_h_Msum) / phy$Cmin

phymsum= phy[!is.na(phy$pMsum),]
#--------------

#Write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy, "MRelevation_Msum.csv")

#===============================
#PLOT

hist(na.omit(phy$pMsum),breaks=8)

#-----------

hl=ggplot(phy, aes(pMsum, fill = Taxa)) + 
  geom_histogram(binwidth = 0.2)+labs(x=expression(Msum / MR[CRB]))+ scale_fill_manual(values = c("darkgreen","blue"))+theme_bw() +xlim(c(0,2.25))+theme(axis.title=element_text(size=rel(1.3)))

#Plot TRAITS


#lower
xyrange= range(c(phy$Tamb_low_Msum, phy$Tmin.use), na.rm=TRUE)
#xyrange[1]= -60
p <- ggplot(data = phy, aes(x = Tamb_low_Msum, y = T10q.min, shape=Taxa, color=as.factor(torpor), size= log(Mass_g))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
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

