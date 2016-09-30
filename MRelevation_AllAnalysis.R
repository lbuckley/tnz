##RUN AFTER EXTRACT TMIN SCOPE

library(geiger)
library(plyr)
library(caper)
library (vegan)
library(ggplot2)
library(car)
library(picante)
library(phytools)
library(ape)
library(grid)
library(MuMIn) #for model averaging

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#define functions
count=function(x) length(na.omit(x))

#=======================================
## READ DATASET
#-----------------------------
#Read physiology data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRexpansibility_Buckleyetal.csv")

#-----------------------
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

#========================================
#Recode North / South range constraints to Cold / Warm boundaries
#C=1: cold constraint; W=1: warm constraint

phy$Cconstrained= NA 
phy$Wconstrained= NA

#Northern hemisphere
inds= which(phy$UpperLat>0 & phy$LowerLat>0  )
phy$Cconstrained[inds]= phy$Nconstrained[inds] 
phy$Wconstrained[inds]= phy$Sconstrained[inds]

#Southern hemisphere
inds= which(phy$UpperLat<0 & phy$LowerLat<0  )
phy$Cconstrained[inds]= phy$Sconstrained[inds] 
phy$Wconstrained[inds]= phy$Nconstrained[inds]

#Crosses hemispheres 
#Change to both cold boundaries, consider constrained if either constraint
inds= which(phy$UpperLat>0 & phy$LowerLat<0  )
con= phy$Sconstrained + phy$Nconstrained
con[which(con==2)]=1

phy$Cconstrained[inds]= con[inds] 
phy$Wconstrained[inds]= 1

#======================================
#ADD METABOLIC CONSTRAINTS

#Calculate conductance (CMIN)
phy$Cmin= abs((0-phy$BMR_mlO2_h)/(phy$Tb-phy$Tlc) )

#specify temperature metrics (Tmin and Tmax) to use
phy$Tmin.use=phy$Tmedian.min
phy$Tmax.use=phy$Tmedian.max

plot(phy$Tlc, phy$Tmin.use)
abline(a=0, b=1)
plot(phy$Tuc, phy$Tmax.use)
abline(a=0, b=1)

#Restrict to species with Tmin below the lower extent of the thermoneutral zone (Tmin <Tlc)
phy= phy[which((phy$Tlc - phy$Tmin.use)>0),]

#-------------------------
#CALCULATE METABOLIC EXPANSIBILITY (MetElev)

#Cold boundary
NBMR= abs(phy$Tlc- phy$Tmin.use)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev= NBMR / phy$BMR_mlO2_h

#Warm boundary
NBMR= abs(phy$Tuc - phy$Tmax.use)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev.hot= NBMR / phy$BMR_mlO2_h

#-------------------------
## PREDICT PHYSIOLOGICAL TEMPERATURE LIMITS
#Species specific environmental temperature than we predict the species can tolerate 

#Cold boundary
phy$Tamb_lowSS= phy$Tlc-(phy$MetElev-1)* phy$BMR_mlO2_h / phy$Cmin

#Warm boundary
phy$Tamb_upSS= phy$Tuc+(phy$MetElev.hot-1)* phy$BMR_mlO2_h / phy$Cmin

#change species with Tmax cooler than the upper end of thermal neutral zone (Tuc>Tmax) to NA
inds= which(phy$Tuc> phy$Tmax.use)
phy[inds,"MetElev.hot"]=NA
phy[inds,"Tamb_upSS"]=NA

#=======================================
#ANALYSIS

#remove constrained range boundaries from analysis
phy$MetElev[which(phy$Cconstrained==1)] = NA
phy$MetElev.hot[which(phy$Wconstrained==1)] = NA

#check outliers
phy[which(phy$MetElev>100),]
#Myotis lucifugus: Hibernate in caves isolated from environments, http://animaldiversity.org/accounts/Myotis_lucifugus/
#Ammodramus savannarum: migratory, http://animaldiversity.org/accounts/Ammodramus_savannarum/
#drop outlier for reasons above
phy=phy[-which(phy$MetElev>100),]

phy[order(phy$MetElev,decreasing=TRUE)[1:10],]
#Planigale tenuirostris, ADW: lives in soil cavities

#--------------------------
#Distribution of metabolic expanisbility at cold range boundary
#Plot density of metabolic expansibility
MetElevs.b= na.omit(phy$MetElev[which(phy$Taxa=="Bird")])
MetElevs.m= na.omit(phy$MetElev[which(phy$Taxa=="Mammal")])

#numbers of species
length(MetElevs.b)
length(MetElevs.m)

#calcualte density
db= density(MetElevs.b,n=50)
dm= density(MetElevs.m)
 
#find peak of density distribution
#peak.b=db$x[which.max(db$y)]
#find mean of two edges of peak
peak.b=mean(db$x[order(db$y, decreasing=TRUE)[1:2]])

peak.m=dm$x[which.max(dm$y)] 

#calculate medians, mean
median(MetElevs.b); mean(MetElevs.b);sd(MetElevs.b)
median(MetElevs.m); mean(MetElevs.m);sd(MetElevs.m)
# se
se=function(x) sd(x)/sqrt(sum(!is.na(x)))
se(MetElevs.m)
se(MetElevs.b)

#-----------
#Distribution of metabolic expanisbility at warm range boundary
MetElevs.b= na.omit(phy$MetElev.hot[which(phy$Taxa=="Bird")])
MetElevs.m= na.omit(phy$MetElev.hot[which(phy$Taxa=="Mammal")])

#numbers of species
length(MetElevs.b)
length(MetElevs.m)

#calculate density
dbh= density(MetElevs.b)
dmh= density(MetElevs.m)

#find peak of density distribution
peak.bh=dbh$x[which.max(dbh$y)]
peak.mh=dmh$x[which.max(dmh$y)]

#calculate medians, mean
median(MetElevs.b); mean(MetElevs.b);sd(MetElevs.b)
median(MetElevs.m); mean(MetElevs.m);sd(MetElevs.m)

#------------------------------
#PLOTS PREDICTED AND OBSERVED THERMAL CONSTRAINTS
# TAMB= TCRIT - (FACTOR-1)BMR/COND

#Predict physiological constraints
#Add temp prediction based on metabolic expansibility across species
phy$Tamb_low=NA
phy$Tamb_up=NA

sel= which(phy$Taxa=="Bird")
phy[sel,"Tamb_low"]= phy$Tlc[sel]-(peak.b-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]
phy[sel,"Tamb_up"]= phy$Tuc[sel]+(peak.bh-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]

sel= which(phy$Taxa=="Mammal")
phy[sel,"Tamb_low"]= phy$Tlc[sel]-(peak.m-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]
phy[sel,"Tamb_up"]= phy$Tuc[sel]+(peak.mh-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]

#account for constrained range edges
phy$Tamb_low[which(phy$Cconstrained==1)] = NA
phy$Tamb_up[which(phy$Wconstrained==1)] = NA

#account for species within themonetural zone (Tmax<Tuc)
inds= which(phy$Tuc> phy$Tmax.use)
phy[inds,"Tamb_up"]=NA

birds= phy[which(phy$Taxa=="Bird"),]
mammals= phy[which(phy$Taxa=="Mammal"),]

#=============================
# PLOTS

#FIGURE 1: PLOT SCOPE AND RESIDUALS

#DENSITY PLOTS
dl=ggplot(phy, aes(MetElev, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at cold range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.2))) + theme(legend.position="none")

du=  ggplot(phy, aes(MetElev.hot, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at warm range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.2)))

#---------------------
#Plot TRAITS
#OPTIONS: color=Order, color=log(Mass..g.), diet, ForStrat, Activity.Nocturnal, Activity.Crepuscular,Activity.Diurnal, Nocturnal         
#Food, Climate, Habitat, Torpor, Clutch_Size, Incubation 

phy$Torpor= as.factor(phy$torpor)
phy$Mass= phy$Mass_g

#lower
xyrange= range(c(phy$Tamb_low, phy$Tmin.use), na.rm=TRUE)
xyrange[1]= -60
p <- ggplot(data = phy, aes(x = Tamb_low, y = Tmin.use, shape=Taxa, color=Torpor, size= log(Mass))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pl= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.2))) + theme(legend.position="none")
#+ facet_wrap(~Taxa)

#upper
xyrange= range(c(phy$Tamb_up, phy$Tmax.use[!is.na(phy$Tamb_up)]), na.rm=TRUE)
p <- ggplot(data = phy, aes(x = Tamb_up, y = Tmax.use, shape=Taxa, color=Torpor, size= log(Mass)))+ xlim(xyrange)+ylim(xyrange)+xlab("Physiological temperature limit (°C)")+ylab("Warm range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pu= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.2)))

#----------
setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("Fig1.pdf", height=10, width=12)

#plot
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(dl,vp=vplayout(1,1))
print(du,vp=vplayout(1,2))
print(pl,vp=vplayout(2,1))
print(pu,vp=vplayout(2,2))

dev.off()

#=======================================
#Models
bird= phy[which(phy$Taxa=="Bird"),]
mamm= phy[which(phy$Taxa=="Mammal"),]

#bird
mod1= lm(bird$MetElev~ log(bird$Mass_g) + bird$diet + bird$Nocturnal)
mod1= lm(bird$MetElev.hot~ log(bird$Mass_g) + bird$diet + bird$Nocturnal +bird$torpor)
summary(mod1)
AIC(mod1)

#mammal
mod1= lm(mamm$MetElev~ log(mamm$Mass_g) + mamm$diet + mamm$Nocturnal + mamm$torpor)
mod1= lm(mamm$MetElev.hot~ log(mamm$Mass_g) + mamm$diet + mamm$Nocturnal +mamm$torpor)

summary(mod1)
AIC(mod1)

#residual plots
crPlots(mod1)
#--------------------------
#BIRDS
bird1= na.omit(bird[,c("MetElev","Mass_g","diet","Nocturnal","torpor")])
mod1= lm(bird1$MetElev~ log(bird1$Mass_g) +bird1$diet + bird1$Nocturnal + bird1$torpor +log(bird1$Mass_g):bird1$diet +log(bird1$Mass_g):bird1$Nocturnal +log(bird1$Mass_g):bird1$torpor+bird1$diet:bird1$Nocturnal +bird1$diet:bird1$torpor + bird1$Nocturnal:bird1$torpor, na.action = "na.fail")

#MODEL SELECTION
d_mod=dredge(mod1)

#extract best 5 models and weights
best.mods= d_mod[1:5,]

ma_mod= model.avg(d_mod)
summary(ma_mod)

##KEEP IMPORTANCE >0.5
mod1= lm(bird$MetElev~ log(bird$Mass_g) +bird$diet + bird$Nocturnal + bird$torpor)
anova(mod1)

#-----------
#MAMMALS
mamm1= na.omit(mamm[,c("MetElev","Mass_g","diet","Nocturnal","torpor")])
mod1= lm(mamm1$MetElev~ log(mamm1$Mass_g) +mamm1$diet + mamm1$Nocturnal + mamm1$torpor +log(mamm1$Mass_g):mamm1$diet +log(mamm1$Mass_g):mamm1$Nocturnal +log(mamm1$Mass_g):mamm1$torpor+mamm1$diet:mamm1$Nocturnal +mamm1$diet:mamm1$torpor + mamm1$Nocturnal:mamm1$torpor, na.action = "na.fail")

#MODEL SELECTION
d_mod=dredge(mod1)

#extract best 5 models and weights
best.mods= d_mod[1:5,]

ma_mod= model.avg(d_mod)
summary(ma_mod)

##KEEP IMPORTANCE >0.5
mod1= lm(mamm$MetElev~ log(mamm$Mass_g) +mamm$diet + mamm$Nocturnal + mamm$torpor +mamm$diet:mamm$torpor +mamm$diet:mamm$Nocturnal +log(mamm$Mass_g):mamm$torpor)
anova(mod1)

#===============================================
#phylogenetic analysis

#READ PHYLOGENY
setwd(paste(mydir,"MRelevation\\Data\\Phylo\\", sep=""))

#mammals
speTree1<-read.nexus("SFritz.tre")

#birds
BirdTree<-read.nexus("my_tree.tre")

#matching
## CHANGE TO ACCOUNT FOR SYNONYMS
bird.l= na.omit(phy[which(phy$Taxa=="Bird"),c("gen_spec","MetElev","Mass_g")])
bird.u= na.omit(phy[which(phy$Taxa=="Bird"),c("gen_spec","MetElev.hot")])

mamm.l= na.omit(phy[which(phy$Taxa=="Mammal"),c("gen_spec","MetElev","Mass_g")])
mamm.u= na.omit(phy[which(phy$Taxa=="Mammal"),c("gen_spec", "MetElev.hot")])

MammTree.l<-drop.tip(speTree1, setdiff(speTree1$tip.label, mamm.l$gen_spec));
match1= match(MammTree.l$tip.label, mamm.l$gen_spec)
mamm.l <- mamm.l[match1,]
rownames(mamm.l)= mamm.l$gen_spec

MammTree.u<-drop.tip(speTree1, setdiff(speTree1$tip.label, mamm.u$gen_spec));
match1= match(MammTree.u$tip.label, mamm.u$gen_spec)
mamm.u <- mamm.u[match1,]
rownames(mamm.u)= mamm.u$gen_spec

BirdTree.l<-drop.tip(BirdTree, setdiff(BirdTree$tip.label, bird.l$gen_spec));
match1= match(BirdTree.l$tip.label, bird.l$gen_spec)
bird.l <- bird.l[match1,]
rownames(bird.l)= bird.l$gen_spec

BirdTree.u<-drop.tip(BirdTree, setdiff(BirdTree$tip.label, bird.u$gen_spec));
match1= match(BirdTree.u$tip.label, bird.u$gen_spec)
bird.u <- bird.u[match1,]
rownames(bird.u)= bird.u$gen_spec

#-----------------------------
#MAMMAL PHYLOGENY PLOT
#Conservatism

setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("FigSphy.pdf", height=8, width=5)

par(mfcol=c(3,2)) 

#mammals
me= log(mamm.l[match(MammTree.l$tip.label,as.character(mamm.l$gen_spec)),"MetElev"])
names(me)= mamm.l$gen_spec
obj<-contMap(MammTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

me= log(mamm.u[match(MammTree.u$tip.label,as.character(mamm.u$gen_spec)),"MetElev.hot"])
names(me)= mamm.u$gen_spec
obj<-contMap(MammTree.u,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")
## FIX TRAIT NAs #,method="anc.ML"

##MASS
me= log(mamm.l[match(MammTree.l$tip.label,as.character(mamm.l$gen_spec)),"Mass_g"])
names(me)= mamm.l$gen_spec
obj<-contMap(MammTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

#-------------
#birds
me= log(bird.l[match(BirdTree.l$tip.label,as.character(bird.l$gen_spec)),"MetElev"])
names(me)= bird.l$gen_spec
obj<-contMap(BirdTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

me= log(bird.u[match(BirdTree.u$tip.label,as.character(bird.u$gen_spec)),"MetElev.hot"])
names(me)= bird.u$gen_spec
obj<-contMap(BirdTree.u,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

##MASS
me= log(bird.l[match(BirdTree.l$tip.label,as.character(bird.l$gen_spec)),"Mass_g"])
names(me)= bird.l$gen_spec
obj<-contMap(BirdTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

dev.off()

## OTHER METRICS
#http://rfunctions.blogspot.com/2014/02/measuring-phylogenetic-signal-in-r.html
#http://webpages.sdsmt.edu/~dbapst/Analyzing_Trait_Evolution_in_R.pdf

#------------------------------------------------------
#PGLS
bird= phy[which(phy$Taxa=="Bird"),]
mamm= phy[which(phy$Taxa=="Mammal"),]

#match
match1= match(speTree1$tip.label, mamm$gen_spec)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
MammTree<-drop.tip(speTree1,speTree1$tip.label[not.matched])

match1= match(BirdTree$tip.label, bird$gen_spec)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
BirdTree<-drop.tip(BirdTree,BirdTree$tip.label[not.matched])

#---------------------------
bird$gen_spec= as.character(bird$gen_spec)
birdc <- comparative.data(BirdTree, bird[,c("gen_spec","MetElev","MetElev.hot", "Mass_g","diet","Nocturnal")], names.col="gen_spec", vcv=TRUE)
mod <- pgls(MetElev ~ log(Mass_g) + diet+ Nocturnal, birdc)
#mod <- pgls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal, birdc)
modn <- pgls(MetElev ~ 1, birdc)
anova.pgls(mod, modn)

mamm$gen_spec= as.character(mamm$gen_spec)
mammc <- comparative.data(MammTree, mamm[,c("gen_spec","MetElev", "Mass_g","diet","Nocturnal","torpor")], names.col="gen_spec", vcv=TRUE) #"MetElev.hot",
mod <- pgls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, mammc)
#mod <- pgls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal+torpor, mammc)

#========================================
#Phylogenetic supplement

#http://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/
#http://www.mpcm-evolution.org/OPM/Chapter5_OPM/OPM_chap5.pdf

#corPagel,corGrafen, corMartins and corBlomberg

## MAMMAL COLD
mamm.l= na.omit(phy[which(phy$Taxa=="Mammal"),c("gen_spec","MetElev","Mass_g","Nocturnal","torpor","diet")])

MammTree.l<-drop.tip(speTree1, setdiff(speTree1$tip.label, mamm.l$gen_spec));
match1= match(MammTree.l$tip.label, mamm.l$gen_spec)
mamm.l <- mamm.l[match1,]
rownames(mamm.l)= mamm.l$gen_spec

#--------------------
#CHECK SIGNAL IN PREDICTORS

#Predictors
phylosignal(mamm.l$Mass_g, MammTree.l)
phylosignal(mamm.l$Nocturnal, MammTree.l)
phylosignal(mamm.l$torpor, MammTree.l)

phylosig(MammTree.l, mamm.l$Mass_g, method="lambda")
phylosig(MammTree.l, mamm.l$Nocturnal, method="lambda")
phylosig(MammTree.l, mamm.l$torpor, method="lambda")

#---------------------
#PHYLO SIGNAL IN RESPONSE
phylosignal(mamm.l$MetElev, MammTree.l)
phylosig(MammTree.l, mamm.l$MetElev, method="lambda")

#PHYLOSIGNAL IN RESIDUALS
#non-phylo models
pglsModel= lm(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, data=mamm.l)
plot(pglsModel)
#------------
summary(pglsModel)
anova(pglsModel)

phylosignal( as.vector(pglsModel$residuals), MammTree.l)
phylosig(MammTree.l, as.vector(pglsModel$residuals), method="lambda")

plot(pglsModel, resid(., type="n")~fitted(.), main="Normalized Residuals v Fitted Values",
     abline=c(0,0))
res <- resid(pglsModel, type="n")
qqnorm(res)
qqline(res)

#------------------

pglsModel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corBrownian(phy = MammTree.l), data = mamm.l, method = "ML")

pglsModel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(0.2, phy = MammTree.l, fixed= FALSE), data = mamm.l, method = "ML")

pglsModel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corMartins(1, phy = MammTree.l), data = mamm.l, method = "ML")

summary(pglsModel)
#------------
#Fit lambda simultaneous
fitPagel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation=corPagel(value=0.8, phy=MammTree.l), data=mamm.l)
intervals(fitPagel, which="var-cov")
summary(fitPagel)

fitPagel0 <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(value = 0, phy = MammTree.l, fixed = TRUE), data = mamm.l) # independence
fitPagel1 <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(value = 1, phy = MammTree.l, fixed = TRUE), data = mamm.l) # Brownian motion

anova(fitPagel, fitPagel0) #almost reject independence
anova(fitPagel, fitPagel1) #reject brownian

#likelihood of lambda values
lambda <- seq(0, 1, length.out = 100)
lik <- sapply(lambda, function(lambda) logLik(gls(MetElev ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(value = lambda, phy = MammTree.l, fixed = TRUE), data = mamm.l)))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
abline(v = fitPagel$modelStruct, col = "red")

#-----------------------------------
## MAMMAL WARM
mamm.l= na.omit(phy[which(phy$Taxa=="Mammal"),c("gen_spec","MetElev.hot","Mass_g","Nocturnal","torpor","diet")])

MammTree.l<-drop.tip(speTree1, setdiff(speTree1$tip.label, mamm.l$gen_spec));
match1= match(MammTree.l$tip.label, mamm.l$gen_spec)
mamm.l <- mamm.l[match1,]
rownames(mamm.l)= mamm.l$gen_spec

#---------------------
#PHYLO SIGNAL IN RESPONSE
phylosignal(mamm.l$MetElev.hot, MammTree.l)
phylosig(MammTree.l, mamm.l$MetElev.hot, method="lambda")

#PHYLOSIGNAL IN RESIDUALS
#non-phylo models
pglsModel= lm(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, data=mamm.l)
plot(pglsModel)
#------------
summary(pglsModel)
anova(pglsModel)

phylosignal( as.vector(pglsModel$residuals), MammTree.l)
phylosig(MammTree.l, as.vector(pglsModel$residuals), method="lambda")

#------------------

pglsModel <- gls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corBrownian(phy = MammTree.l), data = mamm.l, method = "ML")

pglsModel <- gls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(0.2, phy = MammTree.l, fixed= FALSE), data = mamm.l, method = "ML")

pglsModel <- gls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corMartins(1, phy = MammTree.l), data = mamm.l, method = "ML")

summary(pglsModel)
#------------
#Fit lambda simultaneous
fitPagel <- gls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation=corPagel(value=0.8, phy=MammTree.l), data=mamm.l)
intervals(fitPagel, which="var-cov")
summary(fitPagel)

fitPagel0 <- gls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(value = 0, phy = MammTree.l, fixed = TRUE), data = mamm.l) # independence
fitPagel1 <- gls(MetElev.hot ~ log(Mass_g) + diet+ Nocturnal +torpor, correlation = corPagel(value = 1, phy = MammTree.l, fixed = TRUE), data = mamm.l) # Brownian motion

anova(fitPagel, fitPagel0) #almost reject independence
anova(fitPagel, fitPagel1) #reject brownian

#-----------------------------------
## BIRD COLD
mamm.l= na.omit(phy[which(phy$Taxa=="Bird"),c("gen_spec","MetElev","Mass_g","Nocturnal","diet")])

MammTree.l<-drop.tip(BirdTree, setdiff(BirdTree$tip.label, mamm.l$gen_spec));
match1= match(MammTree.l$tip.label, mamm.l$gen_spec)
mamm.l <- mamm.l[match1,]
rownames(mamm.l)= mamm.l$gen_spec

#--------------------
#CHECK SIGNAL IN PREDICTORS

#Predictors
phylosignal(mamm.l$Mass_g, MammTree.l)
phylosignal(mamm.l$Nocturnal, MammTree.l)

phylosig(MammTree.l, mamm.l$Mass_g, method="lambda")
phylosig(MammTree.l, mamm.l$Nocturnal, method="lambda")

#---------------------
#PHYLO SIGNAL IN RESPONSE
phylosignal(mamm.l$MetElev, MammTree.l)
phylosig(MammTree.l, mamm.l$MetElev, method="lambda")

#PHYLOSIGNAL IN RESIDUALS
#non-phylo models
pglsModel= lm(MetElev ~ log(Mass_g) + diet+ Nocturnal, data=mamm.l)
plot(pglsModel)
#------------
summary(pglsModel)
anova(pglsModel)

phylosignal( as.vector(pglsModel$residuals), MammTree.l)
phylosig(MammTree.l, as.vector(pglsModel$residuals), method="lambda")

#------------------

pglsModel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal, correlation = corBrownian(phy = MammTree.l), data = mamm.l, method = "ML")

pglsModel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal, correlation = corPagel(0.2, phy = MammTree.l, fixed= FALSE), data = mamm.l, method = "ML")

pglsModel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal, correlation = corMartins(1, phy = MammTree.l), data = mamm.l, method = "ML")
summary(pglsModel)
#------------
#Fit lambda simultaneous
fitPagel <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal, correlation=corPagel(value=0.8, phy=MammTree.l), data=mamm.l)
intervals(fitPagel, which="var-cov")
summary(fitPagel)

fitPagel0 <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal, correlation = corPagel(value = 0, phy = MammTree.l, fixed = TRUE), data = mamm.l) # independence
fitPagel1 <- gls(MetElev ~ log(Mass_g) + diet+ Nocturnal, correlation = corPagel(value = 1, phy = MammTree.l, fixed = TRUE), data = mamm.l) # Brownian motion

anova(fitPagel, fitPagel0) #almost reject independence
anova(fitPagel, fitPagel1) #reject brownian

#------------------
setwd(paste(mydir,"MRelevation\\Data\\", sep=""))
write.csv(phy, "MRelevation_wTraits.csv")

#=========================================================
#PLOT RELATIONSHIP BETWEEN TEMPERATURE METRICS
setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("FigTemp.pdf", height=10, width=12)

par(mfrow=c(1,2), mar=c(4,4,2,0), oma=c(0,0,0,0), bty="l", lty="solid", cex=1.6, mgp=c(2, 1, 0))

#MIN
plot(mammals$Tmedian.min, mammals$T5q.min, type="p", xlab="Tmin median (°C)", ylab="Tmin metric (°C)")
points(mammals$Tmedian.min, mammals$Tmin, type="p", col="blue")
points(mammals$median.min, mammals$T10q.min, type="p", col="red")
points(mammals$Tmedian.min, mammals$T5q.min, type="p")

points(birds$Tmedian.min, birds$Tmin, type="p", col="blue", pch="*")
points(birds$Tmedian.min, birds$T10q.min, type="p", col="red", pch="*")
points(birds$Tmedian.min, birds$T5q.min, type="p", pch="*")

abline(a=0,b=1, lwd=2)
legend("topleft",pch=1, legend=c("minimum","5th quantile","10th quantile"), col=c("blue","black","red"),bty="n")

#plot(phy$Tmin, phy$Tmedian.min, type="p")
#abline(a=0,b=1, lwd=2)

summary(phy$Tmedian.min-phy$Tmin)
summary(phy$Tsd.min)

#------------------------
#MAX
plot(mammals$Tmedian.max, mammals$T5q.max, type="p", xlab="Tmax median  (°C)", ylab="Tmax metric (°C)")
points(mammals$Tmedian.max, mammals$Tmax, type="p", col="red")
points(mammals$Tmedian.max, mammals$T10q.max, type="p", col="blue")
points(mammals$Tmedian.max, mammals$T5q.max, type="p")

points(birds$Tmedian.max, birds$Tmax, type="p", col="red", pch="*")
points(birds$Tmedian.max, birds$T10q.max, type="p", col="blue", pch="*")
points(birds$Tmedian.max, birds$T5q.max, type="p", pch="*")

abline(a=0,b=1, lwd=2)
legend("topleft",pch=1, legend=c("90th quantile","95th quantile","maximum"), col=c("blue","black","red"),bty="n")

#plot(phy$Tmax, phy$Tmedian.max, type="p")
#abline(a=0,b=1, lwd=2)

summary(phy$Tmedian.max-phy$Tmax)

summary(phy$Tsd.max)

dev.off()

#----------------------
#TEMP ISOCLINEs
summary(mammals$Tsd.min)
summary(mammals$Tmad.min)
summary(birds$Tsd.min)
summary(birds$Tmad.min)

summary(mammals$Tsd.max)
summary(mammals$Tmad.max)
summary(birds$Tsd.max)
summary(birds$Tmad.max)

#--------------------------------------
#SENSITIVITY PLOTs

#METABOLIC EXPANSIBILITY BY MASS
p <- ggplot(data = phy, aes(x = log(Mass_g), y = log(MetElev), shape=Taxa, color=as.factor(torpor) ))+ scale_shape_manual(values = c(1,19))
            
pu= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.2)))
plot(pu)

#-------------------
#Randomize Temp
library(Rmisc)

MEmed= rep(NA,1000)
MEmean = rep(NA,1000)

for(i in 1:1000){
  #phy1=phy[phy$Taxa=="Mammal",]
  phy1=phy[phy$Taxa=="Bird",]
  Trand= sample(phy1$Tmedian.min)
  NBMR= abs(phy$Tlc- Trand)*phy1$Cmin +phy1$BMR_mlO2_h
  MErand<- NBMR / phy1$BMR_mlO2_h
  MEmed[i]= median(MErand)
  MEmean[i]= mean(MErand)
  CI(MErand, ci=0.95)
  #hist(MErand)
}
CI(MEmed, ci=0.95)
CI(MEmean, ci=0.95)

hist(MEmed)
hist(MEmean)

mean(MEmed)
sd(MEmed)

dr= density(phy1$T10q.min) 
dr= density(phy1$Tmedian.min)
plot(dr )

#upper     mean    lower 
#3.982882 3.976690 3.970497

#NBMR= abs(phy$Tuc - Tmax)*phy$Cmin +phy$BMR_mlO2_h
#phy$MetElev.hot<- NBMR / phy$BMR_mlO2_h

MEmed= rep(NA,1000)

for(i in 1:1000){
  Trand= sample(phy$T10q.max)
  NBMR= abs(phy$Tlc- Trand)*phy$Cmin +phy$BMR_mlO2_h
  MErand<- NBMR / phy$BMR_mlO2_h
  MEmed[i]= median(MErand)
}
CI(MEmed, ci=0.95)

dr= density(phy1$T10q.max)
plot(dr )

#-------------------
#Analyze isotherms
#Temp to scange MR_CRB by 10%

#Temp change result in changing MR
fact=0.2
phy$Tdif= -fact*(phy$T10q.min -phy$Tlc - phy$BMR_mlO2_h / phy$Cmin)
hist(phy$Tdif)

plot(phy$Tsd.min, phy$Tdif, xlim=range(0,20), ylim=range(0,20) )
abline(a=0,b=1)

#MR change at Tsd
phy$MRfact= -phy$Tsd.min/(phy$Tmedian.min -phy$Tlc - phy$BMR_mlO2_h / phy$Cmin)
phy$MRfact_max= -phy$Tsd.max/(phy$Tmedian.max - phy$Tuc - phy$BMR_mlO2_h / phy$Cmin)

phy1=phy[phy$Taxa=="Mammal",]
phy1=phy[phy$Taxa=="Bird",]

summary(phy1$MRfact)
summary(na.omit(phy1$MRfact_max[which(phy1$MRfact_max>0)]))

sd(phy1$MRfact)
sd(na.omit(phy1$MRfact_max[which(phy1$MRfact_max>0)]))


