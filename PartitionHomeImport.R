setwd("/home/jana/EddieDir/10K/SU55_import/GenGen1_100/")
setwd("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import//")
popSplit <- read.csv("PopulationSplit.txt")
colnames(popSplit) <- c("group", "Indiv")
ped <- read.table("PedigreeAndGeneticValues_cat.txt", header=TRUE)
#ped$gvNormUnres1 <- (ped$gvNormUnres1 - mean( ped$gvNormUnres1[ped$Generation == 40])) / sd( ped$gvNormUnres1[ped$Generation == 40])

#check the fathers
#gen 40 - 60 
ped <- merge(ped, popSplit, by="Indiv")
ped1 <- ped[ped$Generation %in% 41:60,]
table(ped1$cat[ped1$Indiv %in% ped1$Father[ped1$group == "home"]], ped1$group[ped1$Indiv %in% ped1$Father[ped1$group == "home"]])
ped2 <- ped[ped$Generation %in% 61:80,]
table(ped2$group[ped2$Indiv %in% ped2$Father[ped2$group == "home"]],ped2$sex[ped2$Indiv %in% ped2$Father[ped2$group == "home"]])


pedFathers <- ped[ped$Indiv %in% ped$Father,c("Indiv", "cat", "group")]
colnames(pedFathers) <- c("Father", "catFather", "groupFather")
ped <- merge(ped, pedFathers, by="Father")
pedHome1 <- ped[ped$group == "home" & ped$Generation %in% 41:60,]
pedHome2 <- ped[ped$group == "home" & ped$Generation %in% 61:80,]
pedImport1 <- ped[ped$group == "import" & ped$Generation %in% 41:60,]
pedImport2 <- ped[ped$group == "import" & ped$Generation %in% 61:80,]
table(pedHome1$groupFather)
pedHome1[pedHome1$groupFather == "home",]
table(pedHome1[pedHome1$groupFather == "home","Generation"])
table(pedHome1[pedHome1$groupFather == "import","Generation"])
table(pedHome2[pedHome2$groupFather == "import","Generation"])
table(pedImport1$groupFather)
table(pedImport1$groupFather)
pedPotomci <- ped[ped$cat == "potomciNP" & ped$group == "home",]
table(ped[ped$Indiv %in% pedPotomci$Father, "Generation"])

pedHome <- read.csv("GenPed_EBVhome.txt")
pedImport <- read.csv("GenPed_EBVimport.txt")
genped <- read.csv("GenPed_EBV.txt")

# pedHome <- merge(pedHome, popSplit, by="Indiv", all.x=TRUE)
# length(unique(pedHome$Father[pedHome$Father %in% popSplit$Indiv[popSplit$group == "import"]]))
# importFathers <- unique(pedHome$Father[pedHome$Father %in% popSplit$Indiv[popSplit$group == "import"]])


PED <- merge(ped, popSplit, by="Indiv")
#table(pedHome$group)
PED <- merge(PED, genped[,c("Indiv", "EBV")], by="Indiv")
PED <- PED[PED$Generation > 39,]
PED$Generation <- PED$Generation - 40
#PED$gvNormUnres1 <- (PED$gvNormUnres1 - mean( PED$gvNormUnres1[PED$Generation == 40])) / sd( PED$gvNormUnres1[PED$Generation == 40])

library(AlphaPart)
PED$partGroup <- paste0(PED$group, PED$sex)

Part = AlphaPart(x = PED, sort = FALSE,
                 colId = "Indiv", colFid = "Father", colMid = "Mother",
                 colPath = "partGroup", colBV = "gvNormUnres1")


#PartSummary <- read.table("Documents/PhD/Projects/inProgress/GenomicStrategies_import/PartSummary.txt")
PartSummary = summary(object = Part, by = "Generation")
#plot(PartSummary)

head(Part$gvNormUnres1)
head(PartSummary$gvNormUnres1)




PartHome <- Part
PartHome$gvNormUnres1$Generation <- as.factor(PartHome$gvNormUnres1$Generation)
PartHome$gvNormUnres1 <- PartHome$gvNormUnres1[PartHome$gvNormUnres1$group == "home",]
nrow(PartHome$gvNormUnres1)
table(PartHome$gvNormUnres1$Generation)
PartHomeSummary <- summary(PartHome, by = "Generation")
PartHomeSummary <- summary(PartHome, by = "Generation", FUN="var")
#plot(PartHomeSummary)
head(PartHomeSummary$gvNormUnres1)

PartImport <- Part
PartImport$gvNormUnres1$Generation <- as.factor(PartImport$gvNormUnres1$Generation)
PartImport$gvNormUnres1 <- PartImport$gvNormUnres1[PartImport$gvNormUnres1$group == "import",]
PartImportSummary <- summary(PartImport, by = "Generation")
PartImportSummary <- summary(PartImport, by = "Generation", FUN="var")
covByYear <- PartImport$gvNormUnres1 %>% group_by(Generation) %>% summarise(cov(gvNormUnres1_importM, gvNormUnres1_importF))
VarByYearM <- PartImport$gvNormUnres1 %>% group_by(Generation) %>% summarise(var(gvNormUnres1_importM))
VarByYearF <- PartImport$gvNormUnres1 %>% group_by(Generation) %>% summarise(var(gvNormUnres1_importF))
VarByYearTotal <- PartImport$gvNormUnres1 %>% group_by(Generation) %>% summarise(var(gvNormUnres1))

#plot(PartImportSummary)

write.table(PartImportSummary$gvNormUnres1$abs, "PartImportSummary.txt", quote=FALSE, row.names=FALSE)
write.table(PartHomeSummary$gvNormUnres1$abs, "PartHomeSummary.txt", quote=FALSE, row.names=FALSE)


setwd("/home  ")
PartImport <- read.table("Documents/PhD/Projects/inProgress/GenomicStrategies_import/PartImportSummary.txt", header=TRUE)
PartHome <- read.table("Documents/PhD/Projects/inProgress/GenomicStrategies_import/PartHomeSummary.txt", header=TRUE)

library(reshape)
library(ggplot2)
PartImport_m <- melt(PartImportSummary$gvNormUnres1[,-2], id.vars = "Generation" )
PartImport_m$Population <- "Import"
PartImport_m$PerImport <- 100
PartHome_m <- melt(PartHomeSummary$gvNormUnres1[,-2], id.vars = "Generation" )
PartHome_m$Population <- "Home"
PartHome_m$PerImport <- 100

library(ggplot2)
PartPlot <- rbind(PartImport_m, PartHome_m)

ggplot(data = PartPlot, aes(x = Generation, y = value, group=variable, colour=variable, linetype = variable)) + geom_line(size=1) + facet_grid(cols=vars(Population)) + 
  theme_bw(base_size = 16) + theme(panel.grid.major = element_blank()) + ylab("Genetic gain") + 
  scale_colour_manual("", values = c("black", "#f78ba7", "#68b8de", "#991c3d", "#0e4b7d")) + 
  scale_linetype_manual("", values = c("solid", "solid", "solid", "twodash", "twodash")) + 
  scale_x_discrete(breaks = seq(0, 40, 5))

PartPlotVar <- PartPlot
ggplot(data = PartPlotVar, aes(x = Generation, y = value, group=variable, colour=variable, linetype = variable)) + geom_line(size=1) + facet_grid(cols=vars(Population)) + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank()) + 
  scale_colour_manual("", values = c("black", "#f78ba7", "#68b8de", "#991c3d", "#0e4b7d")) + 
  scale_linetype_manual("", values = c("solid", "solid", "solid", "twodash", "twodash")) + 
  scale_x_discrete(breaks = seq(0, 40, 5))
ggplot(data = PartHome_m, aes(x = Generation, y = value, group=variable, colour=variable)) + geom_line() 
ggplot(data = PartImport_m, aes(x = Generation, y = value, colour=variable)) + geom_line() 


genicVariance <- read.csv("GenicVariance_import.csv")
genicVariance$GenicVariance <- as.numeric(genicVariance$GenicVariance  )
TGVsAll <- read.table("TGVsAll_import80_08062020.csv", header=TRUE)

geneticVar <- ggplot(data = TGVsAll, aes(x=Generation, y = SDSt, colour=Group, group=Group)) + geom_line() + ylab("Genetic standard deviation") + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank()) 
geneticVar <- ggplot(data = TGVsAll, aes(x=Generation, y = SDSt, colour=Group, group=Group)) + geom_line() + ylab("Genetic standard deviation") + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank())
ggplot(data = TGVsAll, aes(x=Generation, y = sd, colour=Group, group=Group)) + geom_line() + ylab("Genetic standard deviation") + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank()) 
genicVar <- ggplot(data = TGVsAll, aes(x=Generation, y = SDGenicSt, colour=Group, group=Group)) + geom_line() + ylab("Genetic standard deviation") + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank())
ggplot(data = TGVsAll, aes(x=Generation, y = SDGenic, colour=Group, group=Group)) + geom_line() + ylab("Genetic standard deviation") + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank()) 
ggplot(data = genicVariance, aes(x=Gen, y = GenVar, colour=Group, group=Group)) + geom_line() + ylab("Genetic standard deviation") + 
  theme_bw(base_size=16) + theme(panel.grid.major = element_blank()) 

library(Rmisc)
multiplot(geneticVar, genicVar)
