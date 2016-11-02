 
#=============================
# PLOTS

#FIGURE 1: PLOT SCOPE AND RESIDUALS

#DENSITY PLOTS
dl=ggplot(phy, aes(MetElev, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at cold range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.2))) + theme(legend.position="none")

du=  ggplot(phy, aes(MetElev.hot, fill = Taxa)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at warm range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4)))

#split birds and mammals
phy.mamm= phy[phy$Taxa=="Mammal",]
phy.bird= phy[phy$Taxa=="Bird",]

dl.mamm=ggplot(phy.mamm, aes(MetElev)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at cold range boundary")+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4)))+ theme(legend.position="none")

dl.bird=ggplot(phy.bird, aes(MetElev)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at cold range boundary")+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4))) + theme(legend.position="none")

#facet
phy$Taxa_f = factor(phy$Taxa, levels=c('Mammal','Bird'))

dl=ggplot(phy, aes(MetElev)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)+xlab("Metabolic elevation at cold range boundary")+ scale_fill_manual(values = c("darkgreen","blue"))+xlim(c(1,10))+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4))) + theme(legend.position="none") +  theme(strip.text.x = element_text(size = rel(1.4)))

dl2= dl+ facet_grid(. ~ Taxa_f)

# +   annotate("text", x = 4, y=13000, label = "ship")

#---------------------
#Plot TRAITS
#OPTIONS: color=Order, color=log(Mass..g.), diet, ForStrat, Activity.Nocturnal, Activity.Crepuscular,Activity.Diurnal, Nocturnal         
#Food, Climate, Habitat, Torpor, Clutch_Size, Incubation 

phy$Torpor= as.factor(phy$torpor)
phy$Mass= phy$Mass_g

#split birds and mammals
phy.mamm= phy[phy$Taxa=="Mammal",]
phy.bird= phy[phy$Taxa=="Bird",]

#lower
xyrange= range(c(phy$Tamb_low, phy$Tmin.use), na.rm=TRUE)
xyrange[1]= -60
p <- ggplot(data = phy, aes(x = Tamb_low, y = Tmin.use, shape=Taxa, color=Torpor, size= log(Mass))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pl= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.2))) + theme(legend.position="none")
#+ facet_wrap(~Taxa)
#-------------
#facet

p <- ggplot(data = phy, aes(x = Tamb_low, y = Tmin.use, color=Torpor, size= log(Mass))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")
pl= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4))) + theme(legend.position="bottom")
pl2= pl + facet_grid(. ~ Taxa_f)

#-------------
#split lower by taxa
p <- ggplot(data = phy.mamm, aes(x = Tamb_low, y = Tmin.use, color=Torpor, size= log(Mass))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pl.mamm= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4))) + theme(legend.position = c(.2, .7))+ theme(legend.text = element_text(size = rel(1.4)), legend.title = element_text(size = rel(1.4)))

p <- ggplot(data = phy.bird, aes(x = Tamb_low, y = Tmin.use, color=Torpor, size= log(Mass))) + xlim(xyrange)+ylim(xyrange) +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pl.bird= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4))) + theme(legend.position="none") +theme(axis.title.x=element_blank(),
axis.title.y=element_blank() )

#upper
xyrange= range(c(phy$Tamb_up, phy$Tmax.use[!is.na(phy$Tamb_up)]), na.rm=TRUE)
p <- ggplot(data = phy, aes(x = Tamb_up, y = Tmax.use, shape=Taxa, color=Torpor, size= log(Mass)))+ xlim(xyrange)+ylim(xyrange)+xlab("Physiological temperature limit (°C)")+ylab("Warm range boundary temperature (°C)")+ scale_shape_manual(values = c(1,19))
#+ scale_color_manual(values = c("gray","darkgreen","purple"))
pu= p + geom_point() + geom_abline(intercept=0, slope=1)+theme_bw()+theme(axis.title=element_text(size=rel(1.4)),axis.text=element_text(size=rel(1.4)))

#----------
#mammal and bird overlaid

#setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
#pdf("Fig1.pdf", height=10, width=12)

##plot
#grid.newpage()
#pushViewport(viewport(layout=grid.layout(2,2)))
#vplayout<-function(x,y)
#  viewport(layout.pos.row=x,layout.pos.col=y)
#print(dl,vp=vplayout(1,1))
#print(du,vp=vplayout(1,2))
#print(pl,vp=vplayout(2,1))
#print(pu,vp=vplayout(2,2))

#dev.off()

#----------------------
#Upper limits in supplement

setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("Fig1.pdf", height=10, width=10)

#plot
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(dl2,vp=vplayout(1,1))
print(pl2,vp=vplayout(2,1))

dev.off()

#=======================================
#--------------------------
#Extract legend

library(grid)
library(gridExtra)

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend<-g_legend(pl.mamm)

#----------------------------
setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("Fig1.pdf", height=10,width=10)

lheight <- sum(legend$height)
p <- arrangeGrob(dl.mamm, pl.mamm, dl.bird, pl.bird, ncol=2, left=textGrob("Absorptivity", rot = 90, gp=gpar(fontsize=20)))

theight <- unit(20, "points")
p <- arrangeGrob(p, textGrob("Physiological Temperature Limit (°C)", gp=gpar(fontsize=20)), heights=unit.c(unit(1, "npc") - theight, theight))
p <- arrangeGrob(p, legend, heights=unit.c(unit(1, "npc") - lheight, lheight), ncol=1)
print(p)

dev.off()

#=========================================

setwd(paste(mydir,"MRelevation\\Figures\\", sep=""))
pdf("FigS1_MRub.pdf", height=10, width=6)

#plot
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(du,vp=vplayout(1,1))
print(pu,vp=vplayout(2,1))

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


