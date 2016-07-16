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

## READ DATA
#-----------------------------
#Read physiology data
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
phy=read.csv("MRelevation_all.csv")
#phy=read.csv("MRelevation_allMASTER.csv")

#=======================================
#CHECK DATA

#check data quality
count=function(x) length(na.omit(x))

aggregate(phy[,c("MetElev")], list(phy$Taxa), FUN=count)

#---------------------
#specify Tmin and Tmax
phy$Tmin.use=phy$T10q.min
phy$Tmax.use=phy$T10q.max

plot(phy$Tlc, phy$Tmin.use)
abline(a=0, b=1)

#Restrict to species with Tmin <Tlc
phy= phy[which((phy$Tlc - phy$Tmin.use)>0),]

#-------------------------
#Update scope and T ambient
NBMR= abs(phy$Tlc- phy$Tmin.use)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev= NBMR / phy$BMR_mlO2_h

NBMR= abs(phy$Tuc - phy$Tmax.use)*phy$Cmin +phy$BMR_mlO2_h
phy$MetElev.hot= NBMR / phy$BMR_mlO2_h

phy$Tamb_lowSS= phy$Tlc-(phy$MetElev-1)* phy$BMR_mlO2_h / phy$Cmin

phy$Tamb_upSS= phy$Tuc+(phy$MetElev.hot-1)* phy$BMR_mlO2_h / phy$Cmin

#change species with Tuc>Tmax to NA
inds= which(phy$Tuc> phy$Tmax.use)
phy[inds,"MetElev.hot"]=NA
phy[inds,"Tamb_upSS"]=NA

#=======================================
#ANALYSIS

phy$scope= phy$MetElev
phy$scope.hot= phy$MetElev.hot

#remove constrained range boundaries
phy$scope[which(phy$Cconstrained==1)] = NA
phy$scope.hot[which(phy$Wconstrained==1)] = NA
# sum(!is.na(phy$scope))
# sum(!is.na(phy$scope.hot))

#ME density
#NEED TO CHECK SCOPES >10 # & phy$scope<10
scopes.b= na.omit(phy$scope[which(phy$Taxa=="Bird")])
scopes.m= na.omit(phy$scope[which(phy$Taxa=="Mammal")])

db= density(scopes.b)
dm= density(scopes.m)

peak.b=db$x[which.max(db$y)]
peak.m=dm$x[which.max(dm$y)]

## HOT SCOPE
scopes.b= na.omit(phy$scope.hot[which(phy$Taxa=="Bird")])
scopes.m= na.omit(phy$scope.hot[which(phy$Taxa=="Mammal")])

dbh= density(scopes.b)
dmh= density(scopes.m)

peak.bh=dbh$x[which.max(dbh$y)]
peak.mh=dmh$x[which.max(dmh$y)]

#------------------------------
#TEMP PLOTS

# TAMB= TCRIT - (FACTOR-1)BMR/COND
#Use peak values

#Add temp prediction
phy$Tamb_low=NA
phy$Tamb_up=NA

sel= which(phy$Taxa=="Bird")
phy[sel,"Tamb_low"]= phy$Tlc[sel]-(peak.b-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]
phy[sel,"Tamb_up"]= phy$Tuc[sel]+(peak.bh-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]

sel= which(phy$Taxa=="Mammal")
phy[sel,"Tamb_low"]= phy$Tlc[sel]-(peak.m-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]
phy[sel,"Tamb_up"]= phy$Tuc[sel]+(peak.mh-1)* phy$BMR_mlO2_h[sel] / phy$Cmin[sel]

#account for constrained
phy$Tamb_low[which(phy$Cconstrained==1)] = NA
phy$Tamb_up[which(phy$Wconstrained==1)] = NA

#account for species with Tmax<Tuc
inds= which(phy$Tuc> phy$Tmax.use)
phy[inds,"Tamb_up"]=NA

birds= phy[which(phy$Taxa=="Bird"),]
mammals= phy[which(phy$Taxa=="Mammal"),]

#check outliers
#drop one outlier
phy=phy[-which(phy$scope>100),]

#=============================
#FIGURE 1: PLOT SCOPE AND RESIDUALS

#DENSITY PLOTS
dl=ggplot(phy, aes(scope, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Factorial scope at cold range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()

du=  ggplot(phy, aes(scope.hot, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Factorial scope at warm range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()

#---------------------
#Plot TRAITS
#OPTIONS: color=Order, color=log(Mass..g.), diet, ForStrat, Activity.Nocturnal, Activity.Crepuscular,Activity.Diurnal, Nocturnal         
#Food, Climate, Habitat, Torpor, Clutch_Size, Incubation 

#lower
xyrange= range(c(phy$Tamb_low, phy$Tmin.use), na.rm=TRUE)
xyrange[1]= -60
p <- ggplot(data = phy, aes(x = Tamb_low, y = Tmin.use, shape=Taxa, color=as.factor(torpor), size= log(Mass_g))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pl= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()
#+ facet_wrap(~Taxa)

#upper
xyrange= range(c(phy$Tamb_up, phy$Tmax.use[!is.na(phy$Tamb_up)]), na.rm=TRUE)
p <- ggplot(data = phy, aes(x = Tamb_up, y = Tmax.use, shape=Taxa, color=as.factor(torpor), size= log(Mass_g)))+ xlim(xyrange)+ylim(xyrange)+xlab("Physiological temperature limit (°C)")+ylab("Warm range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pu= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()

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

#-----------------------------------
#Models
bird= phy[which(phy$Taxa=="Bird"),]
mamm= phy[which(phy$Taxa=="Mammal"),]
  
#bird
mod1= lm(bird$scope~ log(bird$Mass_g) + bird$diet + bird$Nocturnal +bird$torpor)
mod1= lm(bird$scope.hot~ log(bird$Mass_g) + bird$diet + bird$Nocturnal +bird$torpor)
summary(mod1)

#mammal
mod1= lm(mamm$scope~ log(mamm$Mass_g) + mamm$diet + mamm$Nocturnal + mamm$torpor)
mod1= lm(mamm$scope.hot~ log(mamm$Mass_g) + mamm$diet + mamm$Nocturnal +mamm$torpor)

#mod1= lm(mamm$scope~ log(mamm$Mass_g) +mamm$diet + mamm$Nocturnal + mamm$torpor +log(mamm$Mass_g):mamm$diet +log(mamm$Mass_g):mamm$Nocturnal +log(mamm$Mass_g):mamm$torpor+mamm$diet:mamm$Nocturnal +mamm$diet:mamm$torpor + mamm$Nocturnal:mamm$torpor)
summary(mod1)

#residual plots
crPlots(mod1)
#--------------------------
#BIRDS
bird1= na.omit(bird[,c("scope","Mass_g","diet","Nocturnal","torpor")])
mod1= lm(bird1$scope~ log(bird1$Mass_g) +bird1$diet + bird1$Nocturnal + bird1$torpor +log(bird1$Mass_g):bird1$diet +log(bird1$Mass_g):bird1$Nocturnal +log(bird1$Mass_g):bird1$torpor+bird1$diet:bird1$Nocturnal +bird1$diet:bird1$torpor + bird1$Nocturnal:bird1$torpor, na.action = "na.fail")

#MODEL SELECTION
d_mod=dredge(mod1)

#extract best 5 models and weights
best.mods= d_mod[1:5,]

ma_mod= model.avg(d_mod)
summary(ma_mod)

##KEEP IMPORTANCE >0.5
mod1= lm(bird$scope~ log(bird$Mass_g) +bird$diet + bird$Nocturnal + bird$torpor)
anova(mod1)

#-----------
#MAMMALS
mamm1= na.omit(mamm[,c("scope","Mass_g","diet","Nocturnal","torpor")])
mod1= lm(mamm1$scope~ log(mamm1$Mass_g) +mamm1$diet + mamm1$Nocturnal + mamm1$torpor +log(mamm1$Mass_g):mamm1$diet +log(mamm1$Mass_g):mamm1$Nocturnal +log(mamm1$Mass_g):mamm1$torpor+mamm1$diet:mamm1$Nocturnal +mamm1$diet:mamm1$torpor + mamm1$Nocturnal:mamm1$torpor, na.action = "na.fail")

#MODEL SELECTION
d_mod=dredge(mod1)

#extract best 5 models and weights
best.mods= d_mod[1:5,]

ma_mod= model.avg(d_mod)
summary(ma_mod)

##KEEP IMPORTANCE >0.5
mod1= lm(mamm$scope~ log(mamm$Mass_g) +mamm$diet + mamm$Nocturnal + mamm$torpor +mamm$diet:mamm$torpor +mamm$diet:mamm$Nocturnal +log(mamm$Mass_g):mamm$torpor)
anova(mod1)

#===============================================
#phylogenetic analysis

#READ PHYLOGENY
setwd(paste(mydir,"MRelevation\\Data\\Phylo\\", sep=""))

#mammals
# if you want to use a consensus tree
speTree1<-read.nexus("SFritz.tre")

#birds
tree_bird<-read.nexus("my_tree.tre")
#tree_bird2<-read.nexus("Jetz_et_al.tre")

#matching
## CHANGE TO ACCOUNT FOR SYNONYMS
bird.l= na.omit(phy[which(phy$Taxa=="Bird"),c("gen_spec","scope")])
bird.u= na.omit(phy[which(phy$Taxa=="Bird"),c("gen_spec","scope.hot")])

mamm.l= na.omit(phy[which(phy$Taxa=="Mammal"),c("gen_spec","scope","Mass_g")])
mamm.u= na.omit(phy[which(phy$Taxa=="Mammal"),c("gen_spec", "scope.hot")])

MammTree.l<-drop.tip(speTree1, setdiff(speTree1$tip.label, mamm.l$gen_spec));
match1= match(mamm.l$gen_spec, speTree1$tip.label)
matched= which(!is.na(match1))
mamm.l= mamm.l[matched,]

MammTree.u<-drop.tip(speTree1, setdiff(speTree1$tip.label, mamm.u$gen_spec));
match1= match(mamm.u$gen_spec, speTree1$tip.label)
matched= which(!is.na(match1))
mamm.u= mamm.u[matched,]

BirdTree.l<-drop.tip(tree_bird, setdiff(tree_bird$tip.label, bird.l$gen_spec));
match1= match(bird.l$gen_spec, tree_bird$tip.label)
matched= which(!is.na(match1))
bird.l= bird.l[matched,]

BirdTree.u<-drop.tip(tree_bird, setdiff(tree_bird$tip.label, bird.u$gen_spec));
match1= match(bird.u$gen_spec, tree_bird$tip.label)
matched= which(!is.na(match1))
bird.u= bird.u[matched,]

#-----------------------------
#FIG 2: MAMMAL PHYLOGENY
#Conservatism
#PLOT: http://lukejharmon.github.io/ilhabela/instruction/2015/07/05/plotting-methods/

## ADD BARS
#obj<-contMap(anole.tree,exp(svl),plot=FALSE)
#plotTree.wBars(obj$tree,exp(svl),method="plotSimmap",
#               tip.labels=TRUE,fsize=0.7,colors=obj$cols,type="fan",scale=0.002)
#add.color.bar(1.0,obj$cols,title="trait value",lims=obj$lims,prompt=FALSE,
#              x=0.9*par()$usr[1],y=0.9*par()$usr[3])

setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("Fig2.pdf", height=8, width=4)

par(mfcol=c(2,1)) 

#mammals
me= log(mamm.l[match(MammTree.l$tip.label,as.character(mamm.l$gen_spec)),"scope"])
names(me)= mamm.l$gen_spec
obj<-contMap(MammTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

#me= log(mamm.u[match(MammTree.u$tip.label,as.character(mamm.u$gen_spec)),"scope.hot"])
#names(me)= mamm.u$gen_spec
#obj<-contMap(MammTree.u,me,fsize=c(0.1,0.6),outline=FALSE, type="fan")
## FIX TRAIT NAs #,method="anc.ML"

##MASS
me= log(mamm.l[match(MammTree.l$tip.label,as.character(mamm.l$gen_spec)),"Mass_g"])
names(me)= mamm.l$gen_spec
obj<-contMap(MammTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan",ftype="off")

##TORPOP

dev.off()
#-------------
#birds
#me= log(bird.l[match(BirdTree.l$tip.label,as.character(bird.l$gen_spec)),"scope"])
#names(me)= bird.l$gen_spec
#obj<-contMap(BirdTree.l,me,fsize=c(0.1,0.6),outline=FALSE, type="fan")

#me= log(bird.u[match(BirdTree.u$tip.label,as.character(bird.u$gen_spec)),"scope.hot"])
#names(me)= bird.u$gen_spec
#obj<-contMap(BirdTree.u,me,fsize=c(0.1,0.6),outline=FALSE, type="fan")

---------------------
#Phylosignal

phylosignal(mamm.l$scope, MammTree.l)
#phylosignal(mamm.u$scope.hot, MammTree.u)
phylosignal(bird.l$scope, BirdTree.l)
#phylosignal(bird.u$scope.hot, BirdTree.u)

phylosig(MammTree.l, mamm.l$scope, method="lambda") ### SIGNAL
#phylosig(MammTree.u, mamm.u$scope.hot, method="lambda")
phylosig(BirdTree.l, bird.l$scope, method="lambda")
#phylosig(BirdTree.u, bird.u$scope.hot, method="lambda")

## OTHER METRICS
#http://rfunctions.blogspot.com/2014/02/measuring-phylogenetic-signal-in-r.html

#------------------------------------------------------
#PGLS
bird= phy[which(phy$Taxa=="Bird"),]
mamm= phy[which(phy$Taxa=="Mammal"),]

#match
match1= match(speTree1$tip.label, mamm$gen_spec)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
MammTree<-drop.tip(speTree1,speTree1$tip.label[not.matched])

match1= match(tree_bird$tip.label, bird$gen_spec)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))
BirdTree<-drop.tip(tree_bird,tree_bird$tip.label[not.matched])

#---------------------------
bird$gen_spec= as.character(bird$gen_spec)
birdc <- comparative.data(BirdTree, bird[,c("gen_spec","scope","scope.hot", "Mass_g","diet","Nocturnal")], names.col="gen_spec", vcv=TRUE)
mod <- pgls(scope ~ log(Mass_g) + diet+ Nocturnal, birdc)
#mod <- pgls(scope.hot ~ log(Mass_g) + diet+ Nocturnal, birdc)
modn <- pgls(scope ~ 1, birdc)
anova.pgls(mod, modn)

mamm$gen_spec= as.character(mamm$gen_spec)
mammc <- comparative.data(MammTree, mamm[,c("gen_spec","scope", "Mass_g","diet","Nocturnal","torpor")], names.col="gen_spec", vcv=TRUE) #"scope.hot",
mod <- pgls(scope ~ log(Mass_g) + diet+ Nocturnal +torpor, mammc)
#mod <- pgls(scope.hot ~ log(Mass_g) + diet+ Nocturnal+torpor, mammc)

#------------------
setwd(paste(mydir,"MRelevation\\Data\\", sep=""))
write.csv(phy, "MRelevation_wTraits.csv")

#=========================================================
#PLOT RELATIONSHIP BETWEEN TEMPERATURE METRICS

#MIN
plot(phy$T10q.min, phy$T5q.min, type="p")
points(phy$T10q.min, phy$Tmin, type="p", col="blue")
points(phy$T10q.min, phy$Tmedian.min, type="p", col="red")
abline(a=0,b=1, lwd=2)

plot(phy$Tmin, phy$Tmedian.min, type="p")
abline(a=0,b=1, lwd=2)

summary(phy$Tmedian.min-phy$Tmin)

summary(phy$Tsd.min)

#------------------------
#MAX
plot(phy$T10q.max, phy$T5q.max, type="p")
points(phy$T10q.max, phy$Tmax, type="p", col="blue")
points(phy$T10q.max, phy$Tmedian.max, type="p", col="red")
abline(a=0,b=1, lwd=2)

plot(phy$Tmax, phy$Tmedian.max, type="p")
abline(a=0,b=1, lwd=2)

summary(phy$Tmedian.max-phy$Tmax)

summary(phy$Tsd.max)