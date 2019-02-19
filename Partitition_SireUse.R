
################
#Funkcija summarySE
#####################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[colnames(datac) == "mean"] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################
################
#Funkcija summarySE . za SUM
#####################
summarySE1 <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     sum = sum   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[colnames(datac) == "sum"] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################
partSex <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartSex_18022019_gen40.csv")


library(reshape)

#preuredi, da dobiš vse skupine (sex, cat) v en stolpec
partSex <- melt(partSex, measure.vars = c("M", "F"), id.vars = c("Generation", "scenario", "strategy", "rep", "way"))
#primer ene replike
class0 <- partSex[partSex$strategy == "SU55" & partSex$scenario == "Class" & partSex$rep == 0 & partSex$way == "abs",]
library(ggplot2)
ggplot(class0[class0$Generation != 41,], aes(x=Generation, y=value, group = variable, colour=variable)) + geom_line()

#dobi povprečja replik
Means <- summarySE(data = partSex, measurevar = "value", groupvars = c("Generation", "scenario", "strategy", "variable", "way"))
MeanTotal <- summarySE1(data = partSex, measurevar = "value", groupvars = c("Generation", "scenario", "strategy", "way", "rep"))
MeanTotal <- summarySE(data = MeanTotal, measurevar = "value", groupvars = c("Generation", "scenario", "strategy", "way"))
MeanTotal$variable <- "Total"
Means <- rbind(Means, MeanTotal)
Means$PlotGroup <- paste0(Means$scenario, Means$variable)

#plot po spolu
Csex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "abs" & Means$scenario == "Class",], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("PT") + theme(legend.position = "none") + ylim(c(0,5))
PSsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "abs" & Means$scenario == "GenSLO",], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS-PS") + theme(legend.position = "none") + ylim(c(0,5))
GSCsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "abs" & Means$scenario == "OtherCowsGen",], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS-C") + theme(legend.position = "none") + ylim(c(0,5))
GSBDsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "abs" & Means$scenario == "BmGen",], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS-BD") + theme(legend.position = "none") + ylim(c(0,5))
GSsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "abs" & Means$scenario == "Gen",], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS") + ylim(c(0,5))

multiplot(Csex, PSsex, GSCsex, GSBDsex, GSsex, cols=5)

#plot po spolu - relativno
Csex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "rel" & Means$scenario == "Class" & Means$Generation != 40,], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("PT") + theme(legend.position = "none")  + ylim(c(0,1))
PSsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "rel" & Means$scenario == "GenSLO" & Means$Generation != 40,], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS-PS") + theme(legend.position = "none")+ ylim(c(0,1))
GSCsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "rel" & Means$scenario == "OtherCowsGen" & Means$Generation != 40,], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS-C") + theme(legend.position = "none") + ylim(c(0,1))
GSBDsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "rel" & Means$scenario == "BmGen" & Means$Generation != 40,], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS-BD") + theme(legend.position = "none") + ylim(c(0,1))
GSsex <- ggplot(data = Means[Means$strategy == "SU55" & Means$way == "rel" & Means$scenario == "Gen" & Means$Generation != 40,], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path() + ggtitle("GS") + ylim(c(0,1))

multiplot(Csex, PSsex, GSCsex, GSBDsex, GSsex, cols=5)


#CATEGORIJE
class <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_Classgen40_18022019.csv")
class <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_Classgen20_18022019.csv")
cats <- colnames(class)[!(colnames(class) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
class <- melt(class, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))
# class$Action <- NA
# class$Action[class$variable %in% c("cak", "mladi", "pb", "vhlevljeni")] <- "selM"
# class$Action[class$variable %in% c("bik12", "telM", "", "pripust1", "pripust2")] <- "neSelM"
# class$Action[class$variable %in% c("bm", "pBM")] <- "selF"
# class$Action[class$variable %in% c("telF", "pt", "k")] <- "neSelF"
# class$Action[class$variable %in% c("nr", "potomciNP")] <- "NR"
# table(class$Generation[class$variable == "izl"])
# class$Action[class$variable %in% c("izl")] <- "izl"


genSLO <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_GenSLOgen40_18022019.csv")
genSLO <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_GenSLOgen20_18022019.csv")
cats <- colnames(genSLO)[!(colnames(genSLO) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
genSLO <- melt(genSLO, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))

ocowsGen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_OtherCowsGengen40_18022019.csv")
ocowsGen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_OtherCowsGengen20_18022019.csv")
cats <- colnames(ocowsGen)[!(colnames(ocowsGen) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
ocowsGen <- melt(ocowsGen, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))
table(ocowsGen$variable)



BmGen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_BmGengen40_18022019.csv")
BmGen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_BmGengen20_18022019.csv")
cats <- colnames(BmGen)[!(colnames(BmGen) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
BmGen <- melt(BmGen, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))

#class$Action[class$variable %in% c("izl")] <- "izl"


Gen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_Gengen40_18022019.csv")
Gen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_Gengen20_18022019.csv")
cats <- colnames(Gen)[!(colnames(Gen) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
Gen <- melt(Gen, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))


catPart <- rbind(class, genSLO)
catPart <- rbind(catPart, ocowsGen)
catPart <- rbind(catPart, BmGen)
catPart <- rbind(catPart, Gen)

catPart$Action <- NA
catPart$Action[catPart$variable %in% c("cak", "mladi", "vhlevljeni")] <- "selM_PTesting"
catPart$Action[catPart$variable %in% c("pb")] <- "selM_PT"
catPart$Action[catPart$variable %in% c("genTest", "gpb")] <- "selM_GS"
catPart$Action[catPart$variable %in% c("bik12", "telM", "", "pripust1", "pripust2")] <- "neSelM"
catPart$Action[catPart$variable %in% c("bm", "pBM")] <- "selF"
catPart$Action[catPart$variable %in% c("telF", "pt", "k")] <- "neSelF"
catPart$Action[catPart$variable %in% c("nr", "potomciNP")] <- "NR"
table(catPart$Generation[catPart$variable == "izl"])
catPart$Action[catPart$variable %in% c("izl")] <- "izl"
catPart$Action <- as.factor(catPart$Action)
summary(catPart$Action)

#najprej sum po action po replikah
head(catPart)
nrow(catPart)
catPart$scenario <- as.factor(catPart$scenario)
catPart$strategy <- as.factor(catPart$strategy)
catPart$Action <- as.factor(catPart$Action)
table(catPart$Action)
head(catPart)
SumsCat <- aggregate(catPart$value ~ catPart$Generation + catPart$scenario + catPart$strategy + catPart$Action + catPart$way + catPart$rep, FUN="sum")
SumsTotal <- aggregate(catPart$value ~ catPart$Generation + catPart$scenario + catPart$strategy + catPart$way + catPart$rep, FUN="sum")
head(SumsTotal)
head(SumsCat)
nrow(SumsCat)
colnames(SumsCat) <- c("Generation", "scenario", "strategy", "Action", "way", "Rep", "value")
head(SumsCat)
colnames(SumsTotal) <- c("Generation", "scenario", "strategy", "way", "Rep", "value")
SumsTotal$Action <- "Total"
SumsCat <- rbind(SumsCat, SumsTotal)


#sedaj povprečje po action prek replik
MeansCat <- summarySE(data = SumsCat, measurevar = "value", groupvars = c("Generation", "scenario", "strategy", "Action", "way"))
nrow(MeansCat)
head(MeansCat)
MeansCat$PlotGroup <- paste0(MeansCat$scenario, MeansCat$Action)


#CLASS
#relativno - brez 1. generacije
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "Class" & MeansCat$way == "rel" & MeansCat$Generation != 40,], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path()
#absolutno
absC <- ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "Class" & MeansCat$way == "abs",], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() + ylim(c(0,10)) + ggtitle("PT") + 
  scale_colour_manual(values = c("firebrick", "orange", "forestgreen", "steelblue2", "red", "purple", "magenta", "black"))

#GENSLO
#relativno - brez 1. generacije
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "GenSLO" & MeansCat$way == "rel" & MeansCat$Generation != 40,], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path()
#absolutno
absGS <- ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "GenSLO" & MeansCat$way == "abs",], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() + ylim(c(0,10)) + ggtitle("GS-PS" ) + 
  scale_colour_manual(values = c("firebrick", "orange", "forestgreen", "steelblue2","red",  "grey", "purple", "magenta", "black"))

#OTHECOWS GEN
#relativno - brez 1. generacije
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "OtherCowsGen" & MeansCat$way == "rel" & MeansCat$Generation != 40,], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() 
#absolutno
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "OtherCowsGen" & MeansCat$way == "abs",], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() #GEN

#BMGEN
#relativno - brez 1. generacije
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "BmGen" & MeansCat$way == "rel" & MeansCat$Generation != 40,], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() 
#absolutno
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "BmGen" & MeansCat$way == "abs",], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() #GEN

#GEN
#relativno - brez 1. generacije
ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "Gen" & MeansCat$way == "rel" & MeansCat$Generation != 40,], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() +
  scale_colour_manual(values = c("firebrick", "orange", "forestgreen", "steelblue2", "red", "purple", "magenta", "grey", "black"))
  
#absolutno
absG <- ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "Gen" & MeansCat$way == "abs",], 
       aes(x=Generation, y=value, group=PlotGroup, colour=Action)) + geom_path() + ylim(c(0,10)) + ggtitle("GS") + 
  scale_colour_manual(values = c("firebrick", "orange", "forestgreen", "steelblue2", "red", "grey", "purple", "magenta", "black"))


library(Rmisc)
multiplot(absC, absGS, absG, cols=3)

table(class$Action)
class$Action <- as.factor(class$Action)
summary(class$Action)
table(class$variable[is.na(class$Action)])
classMean <-  summarySE(class, measurevar = "value", groupvars = c("Generation", "scenario", "Action", "way"))
table(classMean$Action, classMean$Generation)
ggplot(data = classMean[classMean$scenario == "Class" & classMean$way == "abs",], 
       aes(x=Generation, y=value, group=Action, colour=Action)) + geom_path()
