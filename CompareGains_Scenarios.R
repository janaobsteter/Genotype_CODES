gen <- read.table("~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
genE <- read.table("~/compareSel/PedigreeAndGeneticValues.txt", header=TRUE)
clasE <- read.table("~/compareSel/ClassPedigree.txt", header=TRUE)
bmE <- read.table("~/compareSel/BmGenPedigree.txt", header=TRUE)


genA <- aggregate(gen$gvNormUnres1 ~ gen$Generation, FUN="mean")
colnames(genA) <- c("Generation", "Gen")
genEA <- aggregate(genE$gvNormUnres1 ~ genE$Generation, FUN="mean")
colnames(genEA) <- c("Generation", "GenE")
clasEA <- aggregate(clasE$gvNormUnres1 ~ clasE$Generation, FUN="mean")
colnames(clasEA) <- c("Generation", "ClasE")
bmEA <- aggregate(bmE$gvNormUnres1 ~ bmE$Generation, FUN="mean")
colnames(bmEA) <- c("Generation", "BmGenE")




gain <- merge(genA, genEA, by="Generation")
gain <- merge(gain, clasEA, by="Generation")
gain <- merge(gain, bmEA, by="Generation")

library(reshape)
gainM <- melt(gain, id.vars = "Generation")

ggplot(data=gainM, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_point()


#inspect the genetic variance
genVarG <- read.table("~/TotalGenicAndGeneticVariancesPerGeneration_Gen.txt", header=TRUE)
genVarC <- read.table("~/TotalGenicAndGeneticVariancesPerGeneration_Clas.txt", header=TRUE)
genVarF <- read.table("~/TotalGenicAndGeneticVariancesPerGeneration_FillIn.txt", header=TRUE)
genVarG <- genVarG[genVarG$QtnModel==1,c("Generation", "AdditGeneticVar1")]
genVarC <- genVarC[genVarC$QtnModel==1,c("Generation", "AdditGeneticVar1")]
genVarF <- genVarF[genVarF$QtnModel==1,c("Generation", "AdditGeneticVar1")]
colnames(genVarG) <- c("Generation", "GenVar_G")
colnames(genVarC) <- c("Generation", "GenVar_C")
colnames(genVarF) <- c("Generation", "GenVar_F")
genVar <- merge(genVarG, genVarC, by="Generation")
genVar <- merge(genVar, genVarF, by="Generation", all=TRUE)
genM <- melt(genVar, id.vars = "Generation")
ggplot(data=genM, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path()
