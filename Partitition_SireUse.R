
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
partSex$Total <- partSex$F + partSex$M
#naredi razliko class in genSlo
head(partSex)
class55 <- partSex[partSex$scenario == "Class" & partSex$strategy == "SU55" & partSex$way == "abs", c("Generation", "rep", "F", "M", "Total")]
head(class55)
colnames(class55) <- c("Generation", "Rep", "F_PT", "M_PT", "PT_Total")
head(class55)
gsps55 <- partSex[partSex$scenario == "GenSLO" & partSex$strategy == "SU55" & partSex$way == "abs", c("Generation", "rep", "F", "M", "Total" )]
colnames(gsps55) <- c("Generation", "Rep", "F_GSPS", "M_GSPS", "GSPS_Total")
head(gsps55)
gs55 <- partSex[partSex$scenario == "Gen" & partSex$strategy == "SU55" & partSex$way == "abs", c("Generation", "rep", "F", "M", "Total" )]
colnames(gs55) <- c("Generation", "Rep", "F_GS", "M_GS", "GS_Total")
head(gs55)
dif1 <- merge(class55, gsps55, by=c("Generation", "Rep"))
dif1 <- merge(dif1, gs55, by=c("Generation", "Rep"))
head(dif1)

dif1$DIFF_F_PT_GSPS <- dif1$F_GSPS - dif1$F_PT
dif1$DIFF_M_PT_GSPS <- dif1$M_GSPS - dif1$M_PT
dif1$DIFF_F_PT_GS <- dif1$F_GS - dif1$F_PT
dif1$DIFF_M_PT_GS <- dif1$M_GS - dif1$M_PT
head(dif1)
ggplot(data=dif1, aes(x=PT_Total, y=M_PT)) + geom_line()
summary(glm(dif1$M_PT ~dif1$PT_Total + dif1$Generation))
summary(glm(dif1$M_GSPS ~dif1$GSPS_Total + dif1$Generation))
summary(glm(dif1$M_GS ~dif1$GS_Total + dif1$Generation))
ggplot(data=dif1, aes(x=GS_Total, y=M_GS)) + geom_line()
dif1$Rep <- as.factor(dif1$Rep)

ggplot(data=dif1, aes(x=Generation, y=DIFF_M, group=Rep, colour=Rep)) + geom_path()
ggplot(data=dif1, aes(x=Generation, y=DIFF_F, group=Rep, colour=Rep)) + geom_path()
dif1M <- melt(dif1, measure.vars = c("M_PT", "M_GSPS", "DIFF_M"), id.vars = c("Generation", "Rep"))
ggplot(data=dif1M, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path()

#average
dif1$Generation <- as.factor(dif1$Generation)
dif1Avga <- summarySE(data=dif1, measurevar = c("M_PT"), groupvars = "Generation")[,c(1,3,4)]
colnames(dif1Avga)[3] <- "M_PTsd"
dif1Avgb <- summarySE(data=dif1, measurevar = c("M_GSPS"), groupvars = "Generation")[,c(1,3,4)]
colnames(dif1Avgb)[3] <- "M_GSPSsd"
dif1Avgc <- summarySE(data=dif1, measurevar = c("M_GS"), groupvars = "Generation")[,c(1,3,4)]
colnames(dif1Avgc)[3] <- "M_GSsd"
dif1Avgd <- summarySE(data=dif1, measurevar = c("DIFF_M_PT_GSPS"), groupvars = "Generation")[,c(1,3,4)]
colnames(dif1Avgd)[3] <- "DIFF_M_PT_GSPSsd"
dif1Avge <- summarySE(data=dif1, measurevar = c("DIFF_M_PT_GS"), groupvars = "Generation")[,c(1,3,4)]
colnames(dif1Avge)[3] <- "DIFF_M_PT_GSsd"

diffAvg <- merge(dif1Avga, dif1Avgb, by="Generation")
diffAvg <- merge(diffAvg, dif1Avgc, by="Generation")
diffAvg <- merge(diffAvg, dif1Avgd, by="Generation")
diffAvg <- merge(diffAvg, dif1Avge, by="Generation")
head(diffAvg)

diffSD <- melt(diffAvg[,c(1,3,5,7,9, 11)], id.vars = c("Generation"))
diffAvg <- melt(diffAvg[,c(1,2,4,6,8, 10)], id.vars = c("Generation"))
colnames(diffSD) <- c("Generation", "Variable", "SD")
diffSD$Variable <- as.character(diffSD$Variable)
diffSD$variable <- substr(diffSD$Variable,1,nchar(diffSD$Variable)-2)

diffAvg <- unique(merge(diffAvg, diffSD, by=c("Generation", "variable")))

head(diffSD)
head(diffAvg)
head(diffAvg)
diffAvg$what <- ifelse(diffAvg$variable %in% c("DIFF_M_PT_GSPS", "DIFF_M_PT_GS"), "Difference", "Male Partition")
Mplot <- ggplot(data=diffAvg[diffAvg$Generation != 40,], aes(x=Generation, y=value, group=variable, colour=variable, linetype=what)) +
  geom_ribbon(aes(ymin = value - SD, ymax=value + SD), fill = "grey80", colour=NA) +  geom_line(aes(y=value), size=0.7) + 
  ylab("Genetic mean") + scale_linetype_manual("", values=c("solid", "dashed")) + scale_color_discrete("Scenario") + theme_bw() +
  theme(axis.text=element_text(size=18), 
        axis.title.x=element_text(size=16, vjust=-1),
        axis.title.y=element_text(size=16, vjust=2),
        legend.text=element_text(size=16), legend.title=element_text(size=16), 
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

library(reshape)
#averages for females
fdif1$Generation <- as.factor(dif1$Generation)
fdif1Avga <- summarySE(data=dif1, measurevar = c("F_PT"), groupvars = "Generation")[,c(1,3,4)]
colnames(fdif1Avga)[3] <- "F_PTsd"
fdif1Avgb <- summarySE(data=dif1, measurevar = c("F_GSPS"), groupvars = "Generation")[,c(1,3,4)]
colnames(fdif1Avgb)[3] <- "F_GSPSsd"
fdif1Avgc <- summarySE(data=dif1, measurevar = c("F_GS"), groupvars = "Generation")[,c(1,3,4)]
colnames(fdif1Avgc)[3] <- "F_GSsd"
fdif1Avgd <- summarySE(data=dif1, measurevar = c("DIFF_F_PT_GSPS"), groupvars = "Generation")[,c(1,3,4)]
colnames(fdif1Avgd)[3] <- "DIFF_F_PT_GSPSsd"
fdif1Avge <- summarySE(data=dif1, measurevar = c("DIFF_F_PT_GS"), groupvars = "Generation")[,c(1,3,4)]
colnames(fdif1Avge)[3] <- "DIFF_F_PT_GSsd"

fdiffAvg <- merge(fdif1Avga, fdif1Avgb, by="Generation")
fdiffAvg <- merge(fdiffAvg, fdif1Avgc, by="Generation")
fdiffAvg <- merge(fdiffAvg, fdif1Avgd, by="Generation")
fdiffAvg <- merge(fdiffAvg, fdif1Avge, by="Generation")
head(fdiffAvg)

fdiffSD <- melt(fdiffAvg[,c(1,3,5,7,9, 11)], id.vars = c("Generation"))
fdiffAvg <- melt(fdiffAvg[,c(1,2,4,6,8, 10)], id.vars = c("Generation"))
colnames(fdiffSD) <- c("Generation", "Variable", "SD")
fdiffSD$Variable <- as.character(fdiffSD$Variable)
fdiffSD$variable <- substr(fdiffSD$Variable,1,nchar(fdiffSD$Variable)-2)

fdiffAvg <- unique(merge(fdiffAvg, fdiffSD, by=c("Generation", "variable")))

head(fdiffSD)
head(fdiffAvg)
head(fdiffAvg)
fdiffAvg$what <- ifelse(fdiffAvg$variable %in% c("DIFF_F_PT_GSPS", "DIFF_F_PT_GS"), "Difference", "Female Partition")
Fplot <- ggplot(data=fdiffAvg, aes(x=Generation, y=value, group=variable, colour=variable, linetype=what)) +
  geom_ribbon(aes(ymin = value - SD, ymax=value + SD), fill = "grey80", colour=NA) +  geom_line(aes(y=value), size=0.7) + 
  ylab("Genetic mean") + scale_linetype_manual("", values=c("solid", "dashed")) + scale_color_discrete("Scenario") + theme_bw() +
  theme(axis.text=element_text(size=18), 
        axis.title.x=element_text(size=16, vjust=-1),
        axis.title.y=element_text(size=16, vjust=2),
        legend.text=element_text(size=16), legend.title=element_text(size=16), 
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

multiplot(Mplot, Fplot, cols=2)
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


####naredi ločeno particijo na enem pedigreju
ped <- read.csv("~/SU55_OtherCowsGen0_GenPed_EBV_CATS.txt")[,-1]
ped <- ped[ped$Generation %in% 40:60,]

ped$Cat <-  apply(ped[, paste0("Category", 40:60)], 1, function(x) ifelse(sum(x=="izl", na.rm=TRUE)==21, "izl", tail(x[x != "izl"], n=1) ) )
ped$FirstGen <-  apply(ped[, paste0("Category", 40:60)], 1, function(x) ifelse(any(x=="genTest", na.rm=TRUE), "gen", "clas" ))

ped$GenGen <- apply(ped, 1, function(x) substr(colnames(ped)[which(x == "genTest")], 9, 10)[1])
#tukaj določiš zadnjo generacijo, predno so izločeni
ped$ClassGen <- apply(ped[, paste0("Category", 40:60)], 1, function(x) ifelse(sum(x=="izl", na.rm=TRUE)==21, "izl", tail(which(x != "izl"), n=1) + 39 ) )


sols <- read.csv("~/PedSolutions.txt")         
sols <- sols[sols$Generation %in% 40:60,]
nrow(sols)
nrow(ped)
ped <- merge(ped, sols, by="Indiv")

#določi EBV pri genomskem testiranju
ped$genEBV <- NA
for (row in 1:nrow(ped)) {
  if (ped$FirstGen[row] == "gen") {
    ped$genEBV[row] <-  ped[row, paste0("X", (ped$GenGen[row]))]
  }
}
head(ped)
summary(ped$genEBV)    

#določi zadnjo EBV pred izločenje
ped$classEBV <- NA
for (row in 1:nrow(ped)) {
  if (ped$FirstGen[row] == "clas") {
    ped$classEBV[row] <-  ped[row, paste0("X", (ped$ClassGen[row]))]
  }
}
head(ped)

#združi genomske EBV ni klasične
ped$genPart <- ifelse(ped$FirstGen == "gen", ped$genEBV, ped$classEBV)
summary(ped$classEBV)    

#zdaj pa naredi še, ko so progeno testirani
table(ped$LastCat[ped$FirstGen == "gen"])
ped[ped$FirstGen == "gen" & ped$LastCat == "pb",] 
nrow(ped[ped$FirstGen == "gen" & ped$LastCat == "pb",] )
nrow(ped[ped$LastCat == "pb",] )
ped[ped$LastCat == "pb",]
#zadnje leto, ko so genomsko testirani --> potem gredo v progeno testirane
ped$LastGenGen <- apply(ped[, paste0("Category", 40:60)], 1, function(x) tail(which(x == "gpb"), n=1) + 39)
ped[ped$LastCat == "pb",]

#določi zadnjo EBV preden, da grejo v progeno testirane bike
ped$cakEBV <- NA
for (row in 1:nrow(ped)) {
  if (ped$LastCat[row] == "pb") {
    ped$cakEBV[row] <-  ped[row, paste0("X", (ped$LastGenGen[row]))]
  }
}
ped[ped$LastCat == "pb",]
tail(ped[ped$LastCat == "pb",])

#določi klasične EBV za vse ostale
ped$class1EBV <- NA
for (row in 1:nrow(ped)) {
  if (ped$LastCat[row] != "pb") {
    ped$class1EBV[row] <-  ped[row, paste0("X", (ped$ClassGen[row]))]
  }
}
head(ped)

ped$clasPart <- ifelse(ped$LastCat == "pb", ped$cakEBV, ped$class1EBV)
summary(ped$clasPart)

library(partAGV)
part = data.frame()
PartPed = partAGV(x = as.data.frame(ped), sort = FALSE,
                  colId = "Indiv", colFid = "Father.x", colMid = "Mother.x",
                  colPath = "Sex", colAGV = c("genPart", "clasPart"))
PartPedSummary = summary(object = PartPed, by = "Generation.x")
relTmpG <- PartPedSummary$genPart$rel
relTmpG$way <- "rel"
relTmpG$method <- "gen"
part <- rbind(part, relTmpG)

absTmpG <- PartPedSummary$genPart$abs
absTmpG$way <- "abs"
absTmpG$method <- "gen"

part <- rbind(part, absTmpG)

relTmpC <- PartPedSummary$clasPart$rel
relTmpC$way <- "rel"
relTmpC$method <- "clas"
part <- rbind(part, relTmpC)

absTmpC <- PartPedSummary$clasPart$abs
absTmpC$way <- "abs"
absTmpC$method <- "clas"

part <- rbind(part, absTmpC)
head(part)
part$Total <- part$F +  part$M

partM <- melt(part, measure.vars = c("M", "F", "Total"), id.vars = c("Generation.x", "way", "method"))
partM$group <- paste0(partM$method, partM$variable)
ggplot(data=partM[partM$way == "abs",], aes(x = Generation.x, y = value, group=group, colour=group )) + geom_path() + ylim(0, 12)


#zakaj total ni isti?
aggregate(ped$genPart ~ped$Generation.x, FUN="mean")
aggregate(ped$clasPart ~ped$Generation.x, FUN="mean")
