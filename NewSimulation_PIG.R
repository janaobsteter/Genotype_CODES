library(AlphaSimR)

setwd("/home/jana/Documents/SimulationAlphaPart/NewSimulation/")
homedir <- "/home/jana/Documents/SimulationAlphaPart/NewSimulation/"

nGNMales = 80
nGNFemales = 80

nPNFemales = 500

#create founder population
founderPop = runMacs(nInd = 2*(nGNMales + nGNFemales),
                     nChr = 10,
                     segSites = 100,
                     nThreads = 4,
                     species = "CATTLE")



#set simulation parameters
SimPar = SimParam$new(founderPop)
SimPar$setGender(gender = "yes_sys")
varAm <- matrix(c(1, 0, 0, 1), nrow=2)
varEm <- matrix(c(3, 0, 0, 3), nrow=2)
varPm <- varAm + varEm
diag(varAm) / diag(varPm)

#add the two traits
SimPar$addTraitA(nQtlPerChr = 100, mean = c(0, 0), var = diag(varAm), corA = cov2cor(varAm))

#select Base Population
Base = newPop(founderPop, simParam=SimPar)
Base = setPheno(pop = Base, varE = varEm, simParam = SimPar)

#select GN males and GN females

GNMales = selectInd(pop = Base, gender = "M", 
                    nInd = nGNMales, 
                    use = "pheno", 
                    trait = function (x) rowMeans(scale(x)), 
                    simParam = SimPar)
GNFemales = selectInd(pop = Base, gender = "F", 
                      nInd = nGNFemales, 
                      use = "pheno", 
                      trait = function (x) rowMeans(scale(x)), 
                      simParam = SimPar)


nGenerationBurnIn = 20
nGenerationSel = 20

#burn in - select only in nucleus
for (gen in 1:nGenerationBurnIn) {
  #cross nucleus males and females
  SelCandGN = randCross2(males = GNMales,
                         females = GNFemales,
                         nCrosses = nGNFemales,
                         nProgeny = 16,
                         simParam = SimPar)
  rm(GNMales, GNFemales)
  
  #set phenotype for selection individuals
  SelCandGN <- setPheno(pop = SelCandGN, 
                        varE = varEm,
                        simParam = SimPar)
  #select new nucleus males
  GNMales = selectInd(pop = SelCandGN,
                      nInd = nGNMales,
                      gender = "M",
                      use = "pheno",
                      trait = function (x) rowMeans(scale(x)),
                      simParam = SimPar)
  #select new nucleus females
  GNFemales = selectInd(pop = SelCandGN,
                        nInd = nGNFemales,
                        gender = "F",
                        use = "pheno",
                        trait = function (x) rowMeans(scale(x)),
                        simParam = SimPar)
  #in the last generation, also choose reproductive females
  if (gen == nGenerationBurnIn) {
    UnselBase <- SelCandGN[!(SelCandGN@id %in% c(GNMales@id, GNFemales@id))]
    PNFemales <- selectInd(pop = UnselBase,
                           gender = "F",
                           nInd = nPNFemales,
                           use = "pheno",
                           trait = 1,
                           simParam = SimPar)
    PNFemales@pheno[,2] <- NA
  }
  rm(SelCandGN)
}

PedEval = data.frame()

#do 20 generations of selection
for (gen in 1:nGenerationSel) {
  #first part is the same - do a round of selection in the nucleus
  #cross nucleus males and females
  SelCandGN = randCross2(males = GNMales,
                         females = GNFemales,
                         nCrosses = nGNFemales,
                         nProgeny = 12,
                         simParam = SimPar)
  
  rm(GNMales, GNFemales)
  
  #set phenotype for selection individuals
  SelCandGN <- setPheno(pop = SelCandGN, 
                        varE = varEm,
                        simParam = SimPar)
  
  #write SelCandGN to PedEval
  PedEval <- rbind(PedEval, 
                   data.frame(
                     ID = SelCandGN@id,
                     MID = SelCandGN@mother,
                     FID = SelCandGN@father,
                     Generation = gen,
                     Program = "GN",
                     Gender = SelCandGN@gender,
                     TBV1 = SelCandGN@gv[,1],
                     TBV2 = SelCandGN@gv[,2],
                     EBV1 = NA,
                     EBV2 = NA,
                     Pheno1 = SelCandGN@pheno[,1],
                     Pheno2 = SelCandGN@pheno[,2]
                   ))
  
  #estiate EBVs
  setwd("GN/")
  write.table(PedEval[,c("ID", "FID", "MID")], "Blupf90.ped", quote=FALSE, row.names=FALSE, sep=" ")
  blupf901 <- PedEval[,c("ID", "Pheno1")]
  blupf901$Mean <- 1
  write.table(blupf901, "Blupf901.dat", quote=FALSE, row.names=FALSE, sep=" ", col.names = FALSE, append = file.exists("Blupf901.dat"))
  blupf902 <- PedEval[,c("ID", "Pheno2")]
  blupf902$Mean <- 1
  write.table(blupf902, "Blupf902.dat", quote=FALSE, row.names=FALSE, sep=" ", col.names = FALSE, append = file.exists("Blupf902.dat"))
  
  system("./renumf90 < renumParam1")
  system("./blupf90 renf90.par")
  system("bash MatchAfterRenum.sh")
  rSol <- read.table("renumbered_Solutions")  
  colnames(rSol) <- c("renID", "ID", "EBV1")
  rSol <- rSol[rSol$ID %in% PedEval$ID,]
  PedEval$EBV1 <- rSol$EBV1[match(rSol$ID, PedEval$ID)]  
 
  system("rm renf90.par solutions renumbered_Solutions")
  system("./renumf90 < renumParam2")
  system("./blupf90 renf90.par")
  system("bash MatchAfterRenum.sh")
  rSol <- read.table("renumbered_Solutions")  
  colnames(rSol) <- c("renID", "ID", "EBV2")
  rSol <- rSol[rSol$ID %in% PedEval$ID,]
  PedEval$EBV2 <- rSol$EBV2[match(rSol$ID, PedEval$ID)]
  setwd(homedir)
  
  
  
  SelCandEBVs <- PedEval[PedEval$ID %in% SelCandGN@id,]
  SelCandGN@ebv <- as.matrix(SelCandEBVs[, c("EBV1", "EBV2")])
    
  
  #select new nucleus males
  GNMales = selectInd(pop = SelCandGN,
                      nInd = nGNMales,
                      gender = "M",
                      use = "ebv",
                      trait = function (x) rowMeans(scale(x)),
                      simParam = SimPar)
  #select new nucleus females
  GNFemales = selectInd(pop = SelCandGN,
                        nInd = nGNFemales,
                        gender = "F",
                        use = "ebv",
                        trait = function (x) rowMeans(scale(x)),
                        simParam = SimPar)
  
  rm(SelCandGN)
  
  
  # ##do a round of selection in nucleus - to select the PNFemales
  SelCandPN = randCross2(males = GNMales,
                         females = PNFemales,
                         nCrosses = nPNFemales,
                         nProgeny = 12,
                         simParam = SimPar)
  rm(PNFemales)

  #set phenotypes to SelCandPN
  SelCandPN <- setPheno(pop = SelCandPN,
                        varE = varEm,
                        simParam = SimPar)
  SelCandPN@pheno[,2] <- NA

  #write SelCandGN to PedEval
  PedEval <- rbind(PedEval,
                   data.frame(
                     ID = SelCandPN@id,
                     MID = SelCandPN@mother,
                     FID = SelCandPN@father,
                     Generation = gen + 1,
                     Program = "PN",
                     Gender = SelCandPN@gender,
                     TBV1 = SelCandPN@gv[,1],
                     TBV2 = SelCandPN@gv[,2],
                     EBV1 = NA,
                     EBV2 = NA,
                     Pheno1 = SelCandPN@pheno[,1],
                     Pheno2 = NA
                   ))

  
  #estiate EBVs
  setwd("PN1/")
  write.table(PedEval[,c("ID", "FID", "MID")], "Blupf90.ped", quote=FALSE, row.names=FALSE, sep=" ")
  blupf901 <- PedEval[,c("ID", "Pheno1")]
  blupf901$Mean <- 1
  write.table(blupf901, "Blupf901.dat", quote=FALSE, row.names=FALSE, sep=" ", col.names = FALSE, append = file.exists("Blupf901.dat"))

  
  system("./renumf90 < renumParam1")
  system("./blupf90 renf90.par")
  system("bash MatchAfterRenum.sh")
  rSol <- read.table("renumbered_Solutions")  
  colnames(rSol) <- c("renID", "ID", "EBV1")
  rSol <- rSol[rSol$ID %in% PedEval$ID,]
  PedEval$EBV1 <- rSol$EBV1[match(rSol$ID, PedEval$ID)]  
  
  system("rm renf90.par solutions renumbered_Solutions")
  setwd(homedir)
  
  
  SelCandEBVs <- PedEval[PedEval$ID %in% SelCandPN@id,]
  SelCandPN@ebv<- as.matrix(SelCandEBVs[, c("EBV1")])
  
  #now select the PNFemales
  PNFemales <- selectInd(pop = SelCandPN,
                         nInd = nPNFemales,
                         gender = "F",
                         use = "ebv",
                         trait = 1,
                         simParam = SimPar)
  rm(SelCandPN)
}

library(AlphaPart)
PedEval$ProgramGender <- paste(PedEval$Program, PedEval$Gender, sep="_")

PedEval$TBV1_s <- (PedEval$TBV1 - mean(PedEval$TBV1[PedEval$Generation == 1])) / sd(PedEval$TBV1[PedEval$Generation == 1])
PedEval$TBV2_s <- (PedEval$TBV2 - mean(PedEval$TBV2[PedEval$Generation == 1])) / sd(PedEval$TBV2[PedEval$Generation == 1])
PedEval$Pheno1_s <- (PedEval$Pheno1 - mean(PedEval$Pheno1[PedEval$Generation == 1])) / sd(PedEval$TBV1[PedEval$Generation == 1])
PedEval$Pheno2_s <- (PedEval$Pheno2 - mean(PedEval$Pheno2[PedEval$Generation == 1])) / sd(PedEval$TBV2[PedEval$Generation == 1])

part1 <- AlphaPart(x = PedEval,
                   colId = "ID",
                   colFid = "FID",
                   colMid = "MID",
                   colPath = "Gender",
                   colAGV = c("TBV1", "EBV1"))


part1Sum <- summary(object = part1, by="Generation")
plot1 <- plot(part1Sum)
plot1


library(pedigree)
who <- PedEval$ID %in% unique(PedEval$FID)
Trimed <- trimPed(ped = PedEval, data = who)
PedEval <- PedEval[Trimed,]


#plot genetic trend
library(reshape)
library(ggplot2)
PedEvalM <- melt(PedEval, measure.vars = c("TBV1", "TBV2", "Pheno1", "Pheno2"))
PedEvalM <- melt(PedEval, measure.vars = c("TBV1_s", "TBV2_s", "Pheno1_s", "Pheno2_s"))
PedEvalMA <- summarySE(data = PedEvalM, measurevar = "value", groupvars = c("variable", "Generation"))
ggplot(data = PedEvalMA, aes(x = Generation, value, group=variable, colour=variable)) + geom_line()
