library(pedigree)

ped <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", sep=" ", header=TRUE)
pedC <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/Ped", sep=" ", header=TRUE)
dataG <- ped$cat %in% c("k", "potomciNP", "nr")
dataGI <- ped[ped$cat %in% c("k", "potomciNP", "nr"), "Indiv"]
pedT <- ped[trimPed(ped[,c(2,3,4)], data=dataG, ngenback = 5),]

makeA(pedT[,c(2,3,4)], which=c(pedT$Indiv %in% dataGI))
A <- read.table("A.txt")

makeAinv(pedT[,c(2,3,4)])
Ai <- read.table("Ainv.txt")
nInd <- nrow(pedT)
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd