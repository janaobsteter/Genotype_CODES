library(pedigree)
library(optiSel)
gen <- read.table('//home/jana/bin/AlphaSim1.05Linux/IndForGeno_5gen.txt') #tukaj so krave, pb in potomciNP
colnames(gen) <- "Indiv"
gen$Indiv <- as.numeric(gen$Indiv)

ped <- read.table('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep=" ", header=TRUE)
pedCat <- merge(gen, ped, by="Indiv", all.x=TRUE)
#pedC <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/Ped", sep=" ", header=TRUE)
dataG <- ped$Indiv %in% gen$Indiv
dataGI <- ped[ped$cat %in% c("k", "potomciNP", "nr"), "Indiv"]
#trimPed - tukaj izubereš, koliko generacij nazaj
pedT <- ped[trimPed(ped[,c(2,3,4)], data=dataG, ngenback = 5),]
write.table(pedT[,c(2,3,4)], "/home/jana/Documents/PhD/CompBio/TestingGBLUP/Pedigree.txt", col.names = FALSE, sep=",", row.names=FALSE, quote=FALSE)


ped$Father[ped$Father==0] <- NA
ped$Mother[ped$Mother==0] <- NA
ped$Indiv <- as.numeric(ped$Indiv)
ped$Mother <- as.numeric(ped$Mother)
ped$Father <- as.numeric(ped$Father)
Pedig <- ped[,c("Indiv", "Father", "Mother")]
A <- makeA(Pedig, keep.only=gen$Indiv)

write.table(A, "A_optiSel.txt", quote=FALSE)


#Ai <- read.table("Ainv.txt")
