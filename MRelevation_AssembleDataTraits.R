## ASSEMBLE DATA AND ADD TRAITS

library(stringr)

## set directory for data files
mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\TNZ\\"

#==================================================
#AssEMBLE PHYS DATASETS
setwd(paste(mydir,"\\Data\\", sep=""))

#READ CLIMATE / GEO DATA
Tbird= read.table("Birds_Max_Min_NEW.txt", header=TRUE)
Tmammal= read.table("Mammals_Max_Min_NEW.txt", header=TRUE)
Tall= rbind(Tbird, Tmammal)

#READ PHYS DATA
phy= read.csv("Khaliq_etal_commented_ProcB.csv")
phy$UCT...C.=as.numeric(as.character(tnz$UCT...C.))
phy$gen_spec= paste(phy$Genus,"_",phy$Species,sep="")
names(phy)[6:9]= c("Mass_g","Tlc","Tuc","TNZ")

  #Limit to high quality
  #remove poor quality ##LOOSE 36 SPECIEs
  phy= phy[-which(phy$Quality %in% c("B","C","D")), ]
  
frist= read.csv("FristoeData.csv")

#Match species
#add TNZ  ## Check several UCT 1?
phy$Tb= NA
phy$BMR_mlO2_h= NA
phy$CMIN_mlO2_hC= NA

match1= match(as.character(phy$gen_spec), frist$gen_spec)
matched= which(!is.na(match1))

phy$Tb[matched]= frist$Tb...C.[match1[matched]]
phy$BMR_mlO2_h[matched]= frist$BMR..mlO2.hour.[match1[matched]]
phy$CMIN_mlO2_hC[matched]= frist$CMIN..ml02.hour..C.[match1[matched]]

#---------------------------------------
#ADD DATASETS

#Add species
setwd(paste(mydir,"Data\\AddTNZ\\", sep=""))
riek=read.csv("Rieketal2013.csv")
bozi=read.csv("DataSet TNZ Bozinovic et al.csv")
cant=read.csv("Canterbury2002.csv", na.strings = c("NA","---") )
cant$TRB= gsub("\\*","",cant$TRB)

#Add BMR Data
bmr=read.csv("anage_data.csv")
bmr$genspec= paste(bmr$Genus,bmr$Species, sep=" ")

setwd(paste(mydir,"Data\\", sep=""))
bmr.bird=read.csv("McNab2009_birds.csv")
bmr.mamm=read.csv("McNab2008_mammals.csv")
noct.mamm=read.csv("Bennie_TimePartitioning.csv")

#--------------------------
#MATCH DATA
phy$Species= paste(phy$Genus,phy$Species, sep=" " )
phy.all= phy[,c("Species","Order","Taxa","Mass_g","Tb","Tlc","Tuc", "BMR_mlO2_h","CMIN_mlO2_hC")]
phy.all$BMR_W=NA
phy.all$Source="KhaliqHof"

#riek
riek$Species=str_trim(riek$Species, side = c("both"))
riek$Taxa= "Mammal"
match1= match(riek$Species, phy.all$Species)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy.add= as.data.frame( matrix(NA,nrow=length(not.matched), ncol=ncol(phy.all)) )
names(phy.add)= names(phy.all)
phy.add[, c("Species","Order","Taxa","Mass_g","Tb","Tlc","Tuc") ]<- riek[not.matched,c("Species","Tax","Taxa","BM..g.","Tb..C.","Tlc..C.","Tuc..C.")]
phy.add[, "Source" ]="Riek" 
phy.all= rbind(phy.all,phy.add)

#try to match Tb
#phy.tb= phy.all[is.na(phy.all$Tb) ,]
#match1= match(riek$Species, phy.tb$Species)
#----------
#bozi

bozi$Taxa= "Mammal"
bozi$Order= "Rodentia"
match1= match(bozi$Species, phy.all$Species)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy.add= as.data.frame( matrix(NA,nrow=length(not.matched), ncol=ncol(phy.all)) )
names(phy.add)= names(phy.all)
phy.add[, c("Species","Order","Taxa","Mass_g","Tlc","Tuc") ]<- bozi[not.matched,c("Species","Order","Taxa","mb","Tlc","Tuc")]
phy.add[, "Source" ]="Bozinovic" 
phy.all= rbind(phy.all,phy.add)

#--------
#Canterbury

cant$Taxa= "Bird"
cant$Order= "Passeriformes"

cant$Species=str_trim(cant$Species, side = c("both"))
match1= match(cant$Species, phy.all$Species)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy.add= as.data.frame( matrix(NA,nrow=length(not.matched), ncol=ncol(phy.all)) )
names(phy.add)= names(phy.all)
phy.add[, c("Species","Order", "Taxa","Mass_g","Tlc","BMR_mlO2_h","CMIN_mlO2_hC","Tb") ]<- cant[not.matched,c("Species","Order", "Taxa","M","TLC", "BMR_mLO2_h", "COND_mLO2_hC","Tb")]
phy.add[, "Source" ]="Canterbury" 
phy.all= rbind(phy.all,phy.add)

#========================
#match up BMR
#Add BMR Data: bmr, bmr.bird,bmr.mamm, noct.mamm

bmr.bird$bmr_W= bmr.bird$BMR..kJ.h./3.6
bmr.bird$Species= str_trim(bmr.bird$Species, side = c("both"))

na.inds= which(is.na(phy.all$BMR_W))

match1= match(phy.all$Species[na.inds], bmr.bird$Species)
matched= which(!is.na(match1))
phy.all$BMR_W[na.inds[matched]]<- bmr.bird$bmr_W[match1[matched]] 

#---------------
bmr.mamm$bmr_W= bmr.mamm$BMR..kJ.h./3.6
bmr.mamm$Species= str_trim(bmr.mamm$Species, side = c("both"))

na.inds= which(is.na(phy.all$BMR_W))

match1= match(phy.all$Species[na.inds], bmr.mamm$Species)
matched= which(!is.na(match1))
phy.all$BMR_W[na.inds[matched]]<- bmr.mamm$bmr_W[match1[matched]] 

## MAtch BMR
bmr$Tb= bmr$Temperature..K. - 273

na.inds= which(is.na(phy.all$BMR_W))

match1= match(phy.all$Species[na.inds], bmr$genspec)
matched= which(!is.na(match1))
phy.all$BMR_W[na.inds[matched]]<- bmr$Metabolic.rate..W.[match1[matched]] 

#add Tb
na.inds= which(is.na(phy.all$Tb))

match1= match(phy.all$Species[na.inds], bmr$genspec)
matched= which(!is.na(match1))
phy.all$Tb[na.inds[matched]]<- bmr$Tb[match1[matched]] 

#--------------------
#Convert W to ml/O2h
na.inds= which(is.na(phy.all$BMR_mlO2_h))
phy.all$BMR_mlO2_h[na.inds] =  phy.all$BMR_W[na.inds]*179

#============================================================
#ADD TRAITS

#READ TRAIT DATA
setwd(paste(mydir,"MRelevation\\Data\\Phylo\\", sep=""))

trait_bird= read.csv("Birds_Max_Min.csv", header=TRUE)

setwd(paste(mydir,"\\Data\\EltonTraits\\", sep=""))

elton_m= read.csv("MamFuncDat.csv")
elton_b= read.csv("BirdFuncDat.csv")

#Classify traits
diets= names(elton_m)[4:13]
diets= gsub("Diet.","", diets)
elton_m$diet= diets[apply(elton_m[,4:13], MARGIN=1, FUN=which.max)]

elton_m$diet[which(apply(elton_m[,4:13], MARGIN=1, FUN=max)<50)]= "Omnivore"
elton_m$diet= gsub("Fruit","FruiNect",elton_m$diet)
elton_m$diet= gsub("Nect","FruiNect",elton_m$diet)
elton_m$diet= gsub("FruiFrui","Frui",elton_m$diet)
elton_m$diet= gsub("Inv","Invertebrate",elton_m$diet)
elton_m$diet= gsub("Seed","PlantSeed",elton_m$diet)
elton_m$diet= gsub("PlantO","PlantSeed",elton_m$diet)
elton_m$diet= gsub("Scav","VertFishScav",elton_m$diet)
elton_m$diet= gsub("Vect","VertFishScav",elton_m$diet)
elton_m$diet= gsub("Vend","VertFishScav",elton_m$diet)
elton_m$diet= gsub("Vfish","VertFishScav",elton_m$diet)
elton_m$diet= gsub("Vunk","VertFishScav",elton_m$diet)

#diets= names(elton_b)[10:20]
#elton_b$diet= diets[unlist(apply(elton_b[,10:20], MARGIN=1, FUN=which.max))]

ForStrat= names(elton_b)[24:30]
ForStrat= gsub("ForStrat.","", ForStrat)
inds= as.numeric(apply(elton_b[,24:30], MARGIN=1, FUN=which.max))
elton_b$ForStrat= ForStrat[inds]

#-------------------------------------
#Match data
phy= phy.all
phy$diet=NA; phy$ForStrat=NA; #phy$Activity.Nocturnal=NA; phy$Activity.Crepuscular=NA; phy$Activity.Diurnal=NA; phy$Nocturnal= NA
#phy$Food=NA; phy$Climate=NA; phy$Habitat=NA; phy$Torpor=NA; phy$Clutch_Size=NA; phy$Incubation

#Match species
#Birds
match1= match(as.character(phy$Species), elton_b$Scientific )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy$diet[matched]<- as.character(elton_b$Diet.5Cat[match1[matched]])
phy$Nocturnal[matched]<- elton_b$Nocturnal[match1[matched]]
#-----
match1= match(as.character(phy$gen_spec), trait_bird$species )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

#phy$Food[matched]<- as.character(trait_bird$Food[match1[matched]])
#phy$Climate[matched]<- as.character(trait_bird$Climate[match1[matched]])
#phy$Habitat[matched]<- as.character(trait_bird$Habitat[match1[matched]])
#phy$Torpor[matched]<- as.character(trait_bird$Torpor[match1[matched]])
#phy$Clutch_Size[matched]<- as.character(trait_bird$Clutch_Size[match1[matched]])
#phy$Incubation[matched]<- as.character(trait_bird$Incubation[match1[matched]])

#---------------------
#mammals
match1= match(as.character(phy$Species), elton_m$Scientific )
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

phy$diet[matched]<- elton_m$diet[match1[matched]]
phy$ForStrat[matched]<- as.character(elton_m$ForStrat.Value[match1[matched]])
#phy$Activity.Nocturnal[matched]<- elton_m$Activity.Nocturnal[match1[matched]]
#phy$Activity.Crepuscular[matched]<- elton_m$Activity.Crepuscular[match1[matched]]
#phy$Activity.Diurnal[matched]<- elton_m$Activity.Diurnal[match1[matched]]
phy$Nocturnal[matched]<- elton_m$Activity.Nocturnal[match1[matched]]

#---------------------
#Add torpor / hibernation

setwd(paste(mydir,"MRelevation\\Data\\", sep=""))
torp= read.csv("RufTorporHib.csv")

match1= match(phy$Species, torp$Species)
matched= which(!is.na(match1))
not.matched= which(is.na(match1))

#phy$TorpHib= NA
#phy$TorpHib[matched]<- as.character(torp$TorporHib[match1[matched]])

#====================
#Add Bennie , McNab

setwd(paste(mydir,"Data\\", sep=""))
bmr.bird=read.csv("McNab2009_birds.csv")
bmr.mamm=read.csv("McNab2008_mammals.csv")
noct.mamm=read.csv("Bennie_TimePartitioning.csv")

bmr.bird$Species= str_trim(bmr.bird$Species, side = c("both"))
bmr.mamm$Species= str_trim(bmr.mamm$Species, side = c("both"))
noct.mamm$GenSpec= paste(noct.mamm$Genus, noct.mamm$Species, sep=" ")
#-------------
phy$Climate_McN=NA; phy$Torpor_McN=NA
#phy$Food_McN=NA; phy$Habitat_McN=NA; 

match1= match(phy$Species, bmr.bird$Species)
matched= which(!is.na(match1))
#phy[matched,c("Food_McN")]<- as.character(bmr.bird[match1[matched],c("Foodc")] )
phy[matched,c("Climate_McN")]<- as.character(bmr.bird[match1[matched],c("Climated")] )
#phy[matched,c("Habitat_McN")]<- as.character(bmr.bird[match1[matched],c("Habitate")] )
phy[matched,c("Torpor_McN")]<- as.character(bmr.bird[match1[matched],c("Torpor.f")] )
  
match1= match(phy$Species, bmr.mamm$Species)
matched= which(!is.na(match1))
#phy[matched,c("Food_McN")]<- as.character(bmr.mamm[match1[matched],c("Fooda")] )
phy[matched,c("Climate_McN")]<- as.character(bmr.mamm[match1[matched],c("Climateb")] )
#phy[matched,c("Habitat_McN")]<- as.character(bmr.mamm[match1[matched],c("Habitatc")] )
phy[matched,c("Torpor_McN")]<- as.character(bmr.mamm[match1[matched],c("Torpor.e")] )

#phy$noct=NA
#match1= match(phy$Species, noct.mamm$GenSpec)
#matched= which(!is.na(match1))
#phy[matched,"noct"]<- as.character(noct.mamm[match1[matched],"Behaviour"] )

#--------------
#recode torpor
phy$torpor= NA
#code as 1 in hibernation or torpor
phy[which(phy$Torpor_McN %in% c("HIB","Y","Y?") ),"torpor"]=1
phy[which(phy$Torpor_McN %in% c("N","N?") ),"torpor"]=0

#========================================
#Add synonyms
phy$Spec.syn= as.character(phy$Species)

#Birds
phy$Spec.syn[which(phy$Spec.syn=="Anas rhynchotis")]="Spatula rhynchotis"
phy$Spec.syn[which(phy$Spec.syn=="Collocalia vanikorensis")]="Aerodramus vanikorensis"
phy$Spec.syn[which(phy$Spec.syn=="Aceros plicatus")]="Rhyticeros plicatus"
phy$Spec.syn[which(phy$Spec.syn=="Coturnix chinensis")]="Synoicus chinensis"
phy$Spec.syn[which(phy$Spec.syn=="Aramides cajanea")]="Aramides cajaneus"
phy$Spec.syn[which(phy$Spec.syn=="Gallinula mortierii")]="Tribonyx mortierii"
phy$Spec.syn[which(phy$Spec.syn=="Gallinula ventralis")]="Tribonyx ventralis"
phy$Spec.syn[which(phy$Spec.syn=="Gallirallus owstoni")]="Hypotaenidia owstoni"
phy$Spec.syn[which(phy$Spec.syn=="Porzana cinerea")]="Amaurornis cinerea"
phy$Spec.syn[which(phy$Spec.syn=="Cacatua roseicapilla")]="Eolophus roseicapilla"
phy$Spec.syn[which(phy$Spec.syn=="Otus leucotis")]="Ptilopsis leucotis"
phy$Spec.syn[which(phy$Spec.syn=="Poecile atricapillus")]="Parus atricapillus"
phy$Spec.syn[which(phy$Spec.syn=="Junco hyemalis hyemalis")]="Junco hyemalis"
phy$Spec.syn[which(phy$Spec.syn=="Junco hyemalis oreganus")]="Junco hyemalis"
phy$Spec.syn[which(phy$Spec.syn=="Hesperiphona vespertina")]="Coccothraustes vespertinus"

#Mammals
phy$Spec.syn[which(phy$Spec.syn=="Galerella sanguinea")]="Herpestes sanguineus"
phy$Spec.syn[which(phy$Spec.syn=="Cryptomys bocagei")]="Fukomys bocagei"
phy$Spec.syn[which(phy$Spec.syn=="Cryptomys damarensis")]="Fukomys damarensis"
phy$Spec.syn[which(phy$Spec.syn=="Cryptomys mechowi")]="Fukomys mechowi"
phy$Spec.syn[which(phy$Spec.syn=="Marmosa microtarsus")]="Gracilinanus microtarsus"
phy$Spec.syn[which(phy$Spec.syn=="Ningaui yvonnae")]="Ningaui yvonneae"

phy$Spec.syn[which(phy$Spec.syn=="Echymipera kalabu")]="Echymipera kalubu"
phy$Spec.syn[which(phy$Spec.syn=="Eptesicus vulturnus")]="Vespadelus vulturnus"
#phy$Spec.syn[which(phy$Spec.syn=="Miniopterus schreibersii")]=""
phy$Spec.syn[which(phy$Spec.syn=="Anoura caudifera")]="Anoura caudifer"
phy$Spec.syn[which(phy$Spec.syn=="Choeroniscus perspicillata")]="Carollia perspicillata" #?
phy$Spec.syn[which(phy$Spec.syn=="Artibeus literatus")]="Artibeus lituratus"
#phy$Spec.syn[which(phy$Spec.syn=="Pteronotus parnelli")]=""
phy$Spec.syn[which(phy$Spec.syn=="Rhinonicteris aurantius")]="Rhinonicteris aurantia"
phy$Spec.syn[which(phy$Spec.syn=="Clethrionomys rufocanus")]="Myodes rufocanus"
phy$Spec.syn[which(phy$Spec.syn=="Clethrionomys rutilus")]="Myodes rutilus"
phy$Spec.syn[which(phy$Spec.syn=="Cricetulus triton")]="Tscherskia triton"
phy$Spec.syn[which(phy$Spec.syn=="Apodemus hermonensis")]="Apodemus witherbyi"
#phy$Spec.syn[which(phy$Spec.syn=="Mormota flaviventris")]=""
phy$Spec.syn[which(phy$Spec.syn=="Nannospalax leucodon")]="Spalax leucodon"
#phy$Spec.syn[which(phy$Spec.syn=="Cavia porcellus")]=""
phy$Spec.syn[which(phy$Spec.syn=="Procavia johnstoni")]="Procavia capensis"
phy$Spec.syn[which(phy$Spec.syn=="Hemiechinus aethiopicus")]="Paraechinus aethiopicus"
#phy$Spec.syn[which(phy$Spec.syn=="Equus caballus")]=""
phy$Spec.syn[which(phy$Spec.syn=="Ovis orientalis aries")]="Ovis orientalis"
phy$Spec.syn[which(phy$Spec.syn=="Capra hircus")]="Capra aegagrus"
phy$Spec.syn[which(phy$Spec.syn=="Cervus elaphus c.")]="Cervus elaphus"
phy$Spec.syn[which(phy$Spec.syn=="Alopex lagopus")]="Vulpes lagopus"
phy$Spec.syn[which(phy$Spec.syn=="Proteles cristatus")]="Proteles cristata"
phy$Spec.syn[which(phy$Spec.syn=="Callithrix pygmaea")]="Cebuella pygmaea"
phy$Spec.syn[which(phy$Spec.syn=="Dipodillus dasyurus")]="Gerbillus dasyurus"
phy$Spec.syn[which(phy$Spec.syn=="Micaelamys namaquensis")]="Aethomys namaquensis"

phy$gen_spec= gsub(" ","_",phy$Species)

#==========================================
#Restrict to species with required data
phy= phy[which(!is.na(phy$Tb) & !is.na(phy$BMR_mlO2_h) ),]

phy$Taxa= gsub("Birds","Bird",phy$Taxa)
phy$Taxa= gsub("Mammals","Mammal",phy$Taxa)
phy$Taxa= as.factor(phy$Taxa)

#========================================
#Add range constraint information
phy$Nconstrained=NA
phy$Sconstrained=NA

setwd(paste(mydir,"Data\\", sep=""))
const=read.csv("RangeConstraints.csv")
const$gen_spec= gsub(" ","_", const$Species)

match1= match(as.character(phy$gen_spec), const$gen_spec)
matched= which(!is.na(match1))

phy$Nconstrained[matched]= const$Nconstrained[match1[matched]]
phy$Sconstrained[matched]= const$Sconstrained[match1[matched]]

#restrict to species with at least one unconstrained edge
phy=phy[phy$Nconstrained==0 | phy$Sconstrained==0,]

#========================================

#write out
setwd(paste(mydir,"MRelevation\\Out\\", sep=""))
write.csv(phy, "Phy_all.csv", row.names = FALSE)
#--------------------------------------
