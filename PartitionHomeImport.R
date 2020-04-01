setwd("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/")
popSplit <- read.csv("PopulationSplit.txt")
colnames(popSplit) <- c("group", "Indiv")
ped <- read.table("SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
#ped$gvNormUnres1 <- (ped$gvNormUnres1 - mean( ped$gvNormUnres1[ped$Generation == 40])) / sd( ped$gvNormUnres1[ped$Generation == 40])


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
                 colPath = "partGroup", colAGV = "gvNormUnres1")


#PartSummary <- read.table("Documents/PhD/Projects/inProgress/GenomicStrategies_import/PartSummary.txt")
PartSummary = summary(object = Part, by = "Generation")
#plot(PartSummary)


PartHome <- Part
PartHome$gvNormUnres1$Generation <- as.factor(PartHome$gvNormUnres1$Generation)
PartHome$gvNormUnres1 <- PartHome$gvNormUnres1[PartHome$gvNormUnres1$group == "home",]
nrow(PartHome$gvNormUnres1)
table(PartHome$gvNormUnres1$Generation)
PartHomeSummary <- summary(PartHome, by = "Generation")
#plot(PartHomeSummary)


PartImport <- Part
PartImport$gvNormUnres1$Generation <- as.factor(PartImport$gvNormUnres1$Generation)
PartImport$gvNormUnres1 <- PartImport$gvNormUnres1[PartImport$gvNormUnres1$group == "import",]
PartImportSummary <- summary(PartImport, by = "Generation")
#plot(PartImportSummary)

write.table(PartImportSummary$gvNormUnres1$abs, "PartImportSummary.txt", quote=FALSE, row.names=FALSE)
write.table(PartHomeSummary$gvNormUnres1$abs, "PartHomeSummary.txt", quote=FALSE, row.names=FALSE)


setwd("/home  ")
PartImport <- read.table("Documents/PhD/Projects/inProgress/GenomicStrategies_import/PartImportSummary.txt", header=TRUE)
PartHome <- read.table("Documents/PhD/Projects/inProgress/GenomicStrategies_import/PartHomeSummary.txt", header=TRUE)

PartImport_m <- melt(PartImport[,-2], id.vars = "Generation" )
PartHome_m <- melt(PartHome[,-2], id.vars = "Generation" )

library(ggplot2)
ggplot(data = PartHome_m, aes(x = Generation, y = value, colour=variable)) + geom_line() 
ggplot(data = PartImport_m, aes(x = Generation, y = value, colour=variable)) + geom_line() 
