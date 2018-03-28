library(pedigree)
library(optiSel)
#gen <- read.table('//home/jana/bin/AlphaSim1.05Linux/IndForGeno_5gen.txt') #tukaj so krave, pb in potomciNP
#colnames(gen) <- "Indiv"
#gen$Indiv <- as.numeric(gen$Indiv)

ped <- read.table('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep=" ", header=TRUE)
#pedCat <- merge(gen, ped, by="Indiv", all.x=TRUE)
#pedC <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/Ped", sep=" ", header=TRUE)
dataA <- ped$cat != "izl"
dataAa <- ped$Indiv[ped$cat != "izl"]
#trimPed - tukaj izubereš, koliko generacij nazaj
pedT <- ped[trimPed(ped[,c(2,3,4)], data=dataA, ngenback = 5),]
#write.table(pedT[,c(2,3,4)], "/home/jana/Documents/PhD/CompBio/TestingGBLUP/Pedigree.txt", col.names = FALSE, sep=",", row.names=FALSE, quote=FALSE)


makeA(pedT[,c(2,3,4)], which=c(pedT$Indiv %in% dataAa)) #to so unique elementi plus vsak sam s sabo
#A <- read.table("A.txt")




#Ai <- read.table("Ainv.txt")
