library(pedigree)
ped <- read.table('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep=" ", header=TRUE)

gen <- read.table('//home/jana/bin/AlphaSim1.05Linux/IndForGeno_5gen.txt') #tukaj so krave, pb in potomciNP
colnames(gen) <- "Indiv"
dataG <- ped$Indiv %in% gen$Indiv


pedT <- ped[trimPed(ped[,c(2,3,4)], data=dataG, ngenback = 5),]


makeAinv(pedT[,c(2,3,4)])
Ai <- read.table("Ainv.txt")
nInd <- nrow(pedT)
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd

write.table(Ainv, "Ainv_Work.txt")