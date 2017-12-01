ped <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", sep=" ")
pedDat <- ped[ped$sex=='F',]
pedDat <- pedDat[,c(2,11)]
pedDat$sex <- 2
write.table(pedDat, "/home/jana/Documents/PhD/CompBio/TestingGBLUP/Blupf90.dat", sep=" ", row.names=FALSE, quote=FALSE, col.names=FALSE)

sol <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/renumbered_Solutions", header=FALSE, sep=" ")
#ped <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", sep=" ")
ped <- read.table("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
ped$EBV <- sol$V3
cor(ped$gvNormUnres1, ped$EBV)
library(plyr)
func <- function(ped)
{
  return(data.frame(COR = cor(ped$gvNormUnres1, ped$EBV)))
}

genCor <- ddply(ped, .(ped$Generation), func)
library(ggplot2)
plot1 <- ggplot(genCor, aes(x=genCor$`ped$Generation`, y=genCor$COR)) + geom_point() + xlab("Generation") + ylab("Correlation TGV | EBV") + ggtitle("Renumbered solutions")

#solutiosn wiuthout renumber
noRsol <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/solutions", header=FALSE, skip=1)
ped$nEBV <- noRsol$V4
func <- function(ped)
{
  return(data.frame(COR = cor(ped$gvNormUnres1, ped$nEBV)))
}

genCor <- ddply(ped, .(ped$Generation), func)
plot2 <- ggplot(genCor, aes(x=genCor$`ped$Generation`, y=genCor$COR)) + geom_point()+ xlab("Generation") + ylab("Correlation TGV | EBV") + ggtitle("Non-renumbered solutions")


library(Rmisc)
multiplot(plot1, plot2, cols=2)
