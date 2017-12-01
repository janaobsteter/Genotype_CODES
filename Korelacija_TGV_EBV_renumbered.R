chipInfo <- read.table("~/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/Chip1SnpInformation.txt", header=TRUE)
Map <- chipInfo[,c(2,1,3,3)]
Map$SnpId <- as.character(Map$SnpId)
write.table(Map, "~/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/ChipMap.map", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


cor <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/GenPed_EBV.txt", header=TRUE, sep=",")
sol <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/renumbered_Solutions", sep=" ")
pedO <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", sep=" ", header=TRUE)
ped$EBV <- sol$V3


sol <-  read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/solutions", header=FALSE, skip=1)
colnames(sol) <- c("Trait", "Effect", "RenID", "sol")
ped <-  read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/renadd01.ped", header=FALSE, skip=1)
colnames(ped) <- c("RenID", "F", "M", "a", "b", "c", "d", "e", "f", "Indiv")
SOL <- merge(ped, sol, by="RenID")

PEDSOL <- merge(pedO, SOL, by="Indiv")
cor(PEDSOL$sol, PEDSOL$gvNormUnres1)

PEDSOL30 <- PEDSOL[PEDSOL$Generation==30,]
cor(PEDSOL30$gvNormUnres1, PEDSOL30$sol)
library(plyr)

func <- function(PEDSOL)
{
  return(data.frame(COR = cor(PEDSOL$sol, PEDSOL$gvNormUnres1)))
}

genCor <- ddply(PEDSOL, .(PEDSOL$Generation), func)

library(ggplot2)
ggplot(data=genCor, aes(x=genCor$`PEDSOL$Generation`, y=genCor$COR)) + geom_point()
