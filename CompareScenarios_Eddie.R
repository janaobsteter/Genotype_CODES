setwd("~/Simulation_Rep1/")
pedC = read.table("PedigreeClass.txt", header=TRUE)
pedB = read.table("PedigreeBmGen.txt", header=TRUE)
pedO = read.table("PedigreeOtherCowsGen.txt", header=TRUE)
pedG = read.table("PedigreeGen.txt", header=TRUE)
pedS = read.table("PedigreeGenSLO.txt", header=TRUE)

library(ggplot2)

pedC2 <- pedC[,c("Generation", "gvNormUnres1")]
pedB2 <- pedB[,c("Generation", "gvNormUnres1")]
pedO2 <- pedO[,c("Generation", "gvNormUnres1")]
pedG2 <- pedG[,c("Generation", "gvNormUnres1")]
pedS2 <- pedS[,c("Generation", "gvNormUnres1")]


pedCA <- aggregate(pedC2$gvNormUnres1 ~pedC2$Generation, FUN="mean")
pedBA <- aggregate(pedB2$gvNormUnres1 ~pedB2$Generation, FUN="mean")
pedOA <- aggregate(pedO2$gvNormUnres1 ~pedO2$Generation, FUN="mean")
pedGA <- aggregate(pedG2$gvNormUnres1 ~pedG2$Generation, FUN="mean")
pedSA <- aggregate(pedS2$gvNormUnres1 ~pedS2$Generation, FUN="mean")

colnames(pedCA) <- c("Generation", "Class")
colnames(pedBA) <- c("Generation", "BmGen")
colnames(pedOA) <- c("Generation", "OtherCowsGen")
colnames(pedGA) <- c("Generation", "Gen")
colnames(pedSA) <- c("Generation", "GenSLO")

PED <- merge(pedCA, pedSA, by="Generation")
PED <- merge(PED, pedOA, by="Generation")
PED <- merge(PED, pedBA, by="Generation")
PED <- merge(PED, pedGA, by="Generation")

library(reshape)
PEDm <- melt(PED, id.vars = "Generation")

library(ggplot2)
ggplot(data=PEDm, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path() + ylab("Mean genetic value")


#standardizacia na generacijo 20
PEDs <- PED[PED$Generation %in% 40:60,]
PEDs$ClassSt <-  (PEDs$Class - PEDs$Class[PEDs$Generation == 40]) / genVar$Class[genVar$Generation==40]
PEDs$GenSLOSt <-  (PEDs$GenSLO - PEDs$GenSLO[PEDs$Generation == 40]) / genVar$GenSLO[genVar$Generation==40]
PEDs$OtherCowsGenSt <-  (PEDs$OtherCowsGen - PEDs$OtherCowsGen[PEDs$Generation == 40]) / genVar$OtherCowsGen[genVar$Generation==40]
PEDs$BmGenSt <- (PEDs$BmGen - PEDs$BmGen[PEDs$Generation == 40]) / genVar$BmGen[genVar$Generation==40]
PEDs$GenSt <- (PEDs$Gen - PEDs$Gen[PEDs$Generation == 40]) / genVar$Gen[genVar$Generation==40]



PEDst <- PEDs[,c(1,7,8,9,10,11)]
PEDstM <- melt(PEDst, id.vars = "Generation")
stGainplot <- ggplot(data=PEDstM, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path() + ylab("Standardised mean genetic value") + scale_color_hue("Scenario", labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"))

######
#Variance
######
genVarC <- read.table("~/VarianceClass.txt", header=TRUE)
genVarB <- read.table("~/VarianceBmGen.txt", header=TRUE)
genVarO <- read.table("~/VarianceOtherCowsGen.txt", header=TRUE)
genVarG <- read.table("~/VarianceGen.txt", header=TRUE)
genVarS <- read.table("~/VarianceGenSLO.txt", header=TRUE)

genVarC <- genVarC[genVarC$QtnModel==1,c("Generation", "AdditGeneticVar1")]
genVarS <- genVarS[genVarS$QtnModel==1,c("Generation", "AdditGeneticVar1")]
genVarO <- genVarO[genVarO$QtnModel==1,c("Generation", "AdditGeneticVar1")]
genVarB <- genVarB[genVarB$QtnModel==1,c("Generation", "AdditGeneticVar1")]
genVarG <- genVarG[genVarG$QtnModel==1,c("Generation", "AdditGeneticVar1")]


colnames(genVarC) <- c("Generation", "Class")
colnames(genVarB) <- c("Generation", "BmGen")
colnames(genVarO) <- c("Generation", "OtherCowsGen")
colnames(genVarG) <- c("Generation", "Gen")
colnames(genVarS) <- c("Generation", "GenSLO")

genVar <- merge(genVarC, genVarS, by="Generation")
genVar <- merge(genVar, genVarO, by="Generation")
genVar <- merge(genVar, genVarB, by="Generation")
genVar <- merge(genVar, genVarG, by="Generation")

genM <- melt(genVar, id.vars = "Generation")
varPlot <- ggplot(data=genM[genM$Generation %in% 40:60,], aes(x=Generation, y=value, group=variable, colour=variable)) + 
  geom_path() + ylab("Genetic variance") +
  scale_color_hue("Scenario", labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"))


library(Rmisc)
multiplot(stGainplot, varPlot)
##################################################################
##################################################################
Gen <- aggregate(ped$gvNormUnres1 ~ped$Generation, FUN="mean")
colnames(Gen) <- c("Generation", "TGV")
GenC <- aggregate(pedC$gvNormUnres1 ~pedC$Generation, FUN="mean")
colnames(GenC) <- c("Generation", "TGV_C")

ggplot(data=GenC, aes(x=Generation, y=TGV)) + geom_path()
ggplot(data=Gen, aes(x=Generation, y=TGV)) + geom_path()


pedO <- read.table("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class2/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ", header=TRUE)
GenO <- aggregate(pedO$gvNormUnres1 ~pedO$Generation, FUN="mean")
colnames(GenO) <- c("Generation", "TGV_O")

PED <- merge(GenC, GenO, by="Generation")
library(reshape)
PEDm <- melt(PED, id.vars = "Generation")

ggplot(data=PEDm, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path()

############################
#compare fill ins with different number of sires
fill4 <- read.table("~/PedigreeFillIn4.txt", header=TRUE)
fill8 <- read.table("~/PedigreeFillIn8.txt", header=TRUE)
fill0 <- read.table("~/PedigreeFillIn0.txt", header=TRUE)
class0 <- read.table("~/PedigreeClass0.txt", header=TRUE)

fill4A <- aggregate(fill4$gvNormUnres1 ~fill4$Generation, FUN="mean")
colnames(fill4A) <- c("Generation", "Fill4")
fill8A <- aggregate(fill8$gvNormUnres1 ~fill8$Generation, FUN="mean")
colnames(fill8A) <- c("Generation", "Fill8")
fill0A <- aggregate(fill0$gvNormUnres1 ~fill0$Generation, FUN="mean")
colnames(fill0A) <- c("Generation", "Fill0")
class0A <- aggregate(class0$gvNormUnres1 ~class0$Generation, FUN="mean")
colnames(class0A) <- c("Generation", "Class0")

FILL <- merge(fill4A, fill8A, by="Generation")
FILL <- merge(FILL, fill0A, by="Generation")
FILL <- merge(FILL, class0A, by="Generation")

FILLm <- melt(FILL, id.vars = "Generation")

ggplot(data=FILLm, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path()
