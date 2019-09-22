library(AlphaSimR)


nGNMales = 10 
nGNFemales = 50

nPNFemales = 500
nPNMales = 10

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
varEm <- matrix(c(0.01, 0, 0, 0.01), nrow=2)
varPm <- varAm + varEm
diag(varAm) / diag(varPm)

#add the two traits
SimPar$addTraitA(nQtlPerChr = 100, mean = c(0, 0), var = diag(varAm), corA = cov2cor(varAm))

#select Base Population
Base = newPop(founderPop, simParam=SimPar)
Base = setPheno(pop = Base, varE = varEm, simParam = SimPar)

#select GN males and GN femals

GNMales = selectInd(pop = Base, gender = "M", 
                    nInd = nGNMales, 
                    use = "pheno", 
                    trait = function (x) rowMeans(x), 
                    simParam = SimPar)
GNFemales = selectInd(pop = Base, gender = "F", 
                      nInd = nGNFemales, 
                      use = "pheno", 
                      trait = function (x) rowMeans(x), 
                      simParam = SimPar)


nGenerationBurnIn = 20
nGenerationSel = 20

#burn in - select only in nucleus
for (gen in 1:nGenerationBurnIn) {
  #cross nucleus males and females
  SelCandGN = randCross2(males = GNMales,
                         females = GNFemales,
                         nCrosses = sum(Base@gender == "F"),
                         nProgeny = 20,
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
                         nCrosses = sum(Base@gender == "F"),
                         nProgeny = 20,
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
                     Pheno1 = SelCandGN@pheno[,1],
                     Pheno2 = SelCandGN@pheno[,2]
                   ))
  
  
  #select new nucleus males
  GNMales = selectInd(pop = SelCandGN,
                      nInd = nGNMales,
                      gender = "M",
                      use = "pheno",
                      trait = 1, #function (x) rowMeans(scale(x)),
                      simParam = SimPar)
  #select new nucleus females
  GNFemales = selectInd(pop = SelCandGN,
                        nInd = nGNFemales,
                        gender = "F",
                        use = "pheno",
                        trait = 1, #function (x) rowMeans(scale(x)),
                        simParam = SimPar)
  
  rm(SelCandGN)
  
  
  # ##do a round of selection in nucleus - to select the PNFemales
  SelCandPN = randCross2(males = GNMales,
                         females = PNFemales,
                         nCrosses = nPNFemales * 12,
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
                     Pheno1 = SelCandPN@pheno[,1],
                     Pheno2 = NA
                   ))

  #now select the PNFemales
  PNFemales <- selectInd(pop = SelCandPN,
                         nInd = nPNFemales,
                         gender = "F",
                         use = "pheno",
                         trait = 1,
                         simParam = SimPar)
  rm(SelCandPN)
}

library(AlphaPart)
?AlphaPart
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
                   colAGV = c("TBV1", "Pheno1"))


part1Sum <- summary(object = part1, by="Generation")
plot1 <- plot(part1Sum)

#plot genetic trend
library(reshape)
library(ggplot2)
PedEvalM <- melt(PedEval, measure.vars = c("TBV1", "TBV2", "Pheno1", "Pheno2"))
PedEvalM <- melt(PedEval, measure.vars = c("TBV1_s", "TBV2_s", "Pheno1_s", "Pheno2_s"))
PedEvalMA <- summarySE(data = PedEvalM, measurevar = "value", groupvars = c("variable", "Generation"))
ggplot(data = PedEvalMA, aes(x = Generation, value, group=variable, colour=variable)) + geom_line()
