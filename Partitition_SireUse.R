partSex <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartSex_11022019.csv")


library(reshape)

#preuredi, da dobiš vse skupine (sex, cat) v en stolpec
partSex <- melt(partSex, measure.vars = c("M", "F"), id.vars = c("Generation", "scenario", "strategy", "rep", "way"))
#primer ene replike
class0 <- partSex[partSex$strategy == "SU55" & partSex$scenario == "Class" & partSex$rep == 0 & partSex$way == "abs",]
ggplot(class0[class0$Generation != 41,], aes(x=Generation, y=value, group = variable, colour=variable)) + geom_line()

#dobi povprečja replik
Means <- summarySE(data = partSex, measurevar = "value", groupvars = c("Generation", "scenario", "strategy", "variable", "way"))
Means$PlotGroup <- paste0(Means$scenario, Means$variable)

#plot po spolu
ggplot(data = Means[Means$strategy == "SU55" & Means$way == "rel",], aes(x=Generation, y=value, group=PlotGroup, colour=variable, linetype = scenario)) + geom_path()




#CATEGORIJE
class <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_Class_11022019.csv")
cats <- colnames(class)[!(colnames(class) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
class <- melt(class, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))

genSLO <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_GenSLO_11022019.csv")
cats <- colnames(genSLO)[!(colnames(genSLO) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
genSLO <- melt(genSLO, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))

ocowsGen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_OtherCowsGen_11022019.csv")
cats <- colnames(ocowsGen)[!(colnames(ocowsGen) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
ocowsGen <- melt(ocowsGen, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))

BmGen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_BmGen_11022019.csv")
cats <- colnames(BmGen)[!(colnames(BmGen) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
BmGen <- melt(BmGen, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))

Gen <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PartCat_Gen_11022019.csv")
cats <- colnames(Gen)[!(colnames(Gen) %in% c("Generation", "N", "Sum", "way", "scenario", "strategy", "rep")) ]
Gen <- melt(Gen, measure.vars = cats, id.vars = c("Generation", "scenario", "strategy", "rep", "way"))


catPart <- rbind(class, genSLO)
catPart <- rbind(catPart, ocowsGen)
catPart <- rbind(catPart, BmGen)
catPart <- rbind(catPart, Gen)

MeansCat <- summarySE(data = catPart, measurevar = "value", groupvars = c("Generation", "scenario", "strategy", "variable", "way"))
MeansCat$PlotGroup <- paste0(MeansCat$scenario, MeansCat$variable)

#poglej en scenrija
oneGen <- MeansCat[MeansCat$scenario == "Class" & MeansCat$strategy == "SU55" & MeansCat$way == "abs" & MeansCat$Generation == 60,]
oneGen[order(-oneGen$value),]


ggplot(data = MeansCat[MeansCat$strategy == "SU55" & MeansCat$scenario == "Class" & MeansCat$way == "abs",], aes(x=Generation, y=value, group=PlotGroup, colour=variable)) + geom_path()
