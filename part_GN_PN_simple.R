# 
# ---- Environment ----

rm(list = ls())

library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")
library(pedigree)

setwd("/home/jana/Documents/PhD/SimulationAlphaPart")
homedir = "/home/jana/Documents/PhD/SimulationAlphaPart"
# ---- General parameters ----

nGNMales   =  25
nGNFemales = 500

nPNMales       =  100
nPNFemales     = 5000
pPNFemalesPure = 0.15

nGenerationBurn = 20
nGenerationEval = 20

GenMeanCols = c("GenMeanT1", "GenMeanT2", "GenMeanI")
GenVarCols  = c("GenVarT1",  "GenVarT2",  "GenVarI")

# ---- Base population genomes ----

founderPop = runMacs(nInd = nGNMales + nGNFemales,
                     nChr = 10,
                     segSites = 1000,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")

###################################################################################
###################################################################################
# ---- Simulation/Base population parameters ----

SP = SimParam$new(founderPop)
# VarA = matrix(data = c(1.0, 0.1, 0.1, 1.0), nrow = 2); cov2cor(VarA)
VarA = matrix(data = c(1.0, 0.0, 0.0, 1.0), nrow = 2); cov2cor(VarA)
VarE = matrix(data = c(3.0, 0.0, 0.0, 9.0), nrow = 2); cov2cor(VarE)
VarP = VarA + VarE; diag(VarA) / diag(VarP)
SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0), var = diag(VarA), cor = cov2cor(VarA))
# SP$addSnpChip(nSnpPerChr = 1000)
SP$setGender(gender = "yes_rand")

# ---- Base GN population ----

GN = newPop(founderPop)
GNMales   = GN[GN@gender == "M"]
GNFemales = GN[GN@gender == "F"]
rm(GN)

###################################################################################
###################################################################################

# ---- GN burn-in ----

DataBurn = tibble(Generation = 1:nGenerationBurn,
                  GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                  GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
for (Generation in 1:nGenerationBurn) {
  # Generation = 1

  # Mate
  SelCand = randCross2(females = GNFemales, males = GNMales,
                       nCrosses = GNFemales@nInd, nProgeny = 12)

  # Save metrics
  for (gender in c("F", "M")) {
    DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }
  
  
  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)

  # Evaluate (could run BLUP)
  SelCand@ebv = SelCand@pheno

  # Select
  GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                        use = "ebv", trait = function(x) rowMeans(x))
  GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                        use = "ebv", trait = function(x) rowMeans(x))
}

# Plot genetic means
DataBurn %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic mean")

# Plot genetic variances
DataBurn %>% gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic variance")

# ---- Base PN population ----

GNInd = c(GNFemales@id, GNMales@id)
SelCand = SelCand[!(SelCand@id %in% GNInd)]; SelCand@nInd
# Cheat here - consider all animals are females
SelCand@gender[] = "F"
PNFemales = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                      use = "ebv", trait = function(x) rowMeans(x))
# meanG(GNMales); meanG(GNFemales); meanG(PNFemales)


# ---- Evaluation ----

DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                    Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                    GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                    GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
DataEvalGN$Program = "GN"
PedEvalGN = rbind(tibble(Generation = 0,
                         IId        = GNMales@id,
                         FId        = NA,
                         MId        = NA,
                         Gender     = GNMales@gender,
                         Program    = "GN",
                         PhenoT1    = GNMales@pheno[,1],
                         PhenoT2    = GNMales@pheno[,2],
                         EbvT1      = GNMales@ebv[, 1],
                         EbvT2      = GNMales@ebv[, 2],
                         TbvT1      = GNMales@gv[, 1],
                         TbvT2      = GNMales@gv[, 2]),
                  tibble(Generation = 0,
                         IId        = GNFemales@id,
                         FId        = NA,
                         MId        = NA,
                         Gender     = GNFemales@gender,
                         Program    = "GN",
                         PhenoT1    = GNFemales@pheno[,1],
                         PhenoT2    = GNFemales@pheno[,2],
                         EbvT1      = GNFemales@ebv[, 1],
                         EbvT2      = GNFemales@ebv[, 2],
                         TbvT1      = GNFemales@gv[, 1],
                         TbvT2      = GNFemales@gv[, 2]))

DataEvalPN1 = DataEvalGN
DataEvalPN1$Program = "PN1"
PedEvalPN1 = rbind(tibble(Generation = 0,
                          IId        = GNMales@id,
                          FId        = NA,
                          MId        = NA,
                          Gender     = GNMales@gender,
                          Program    = "GN",
                          PhenoT1    = GNMales@pheno[,1],
                          PhenoT2    = GNMales@pheno[,2],
                          EbvT1      = GNMales@ebv[, 1],
                          EbvT2      = GNMales@ebv[, 2],
                          TbvT1      = GNMales@gv[, 1],
                          TbvT2      = GNMales@gv[, 2]),
                   tibble(Generation = 0,
                          IId        = PNFemales@id,
                          FId        = NA,
                          MId        = NA,
                          Gender     = PNFemales@gender,
                          Program    = "PN1",
                          PhenoT1    = PNFemales@pheno[,1],
                          PhenoT2    = PNFemales@pheno[,2],
                          EbvT1      = PNFemales@ebv[, 1],
                          EbvT2      = PNFemales@ebv[, 2],
                          TbvT1      = PNFemales@gv[, 1],
                          TbvT2      = PNFemales@gv[, 2]))

DataEvalPN2 = DataEvalGN
DataEvalPN2$Program = "PN2"
PedEvalPN2 = rbind(tibble(Generation = 0,
                          IId        = GNMales@id,
                          FId        = NA,
                          MId        = NA,
                          Gender     = GNMales@gender,
                          Program    = "GN",
                          PhenoT1    = GNMales@pheno[,1],
                          PhenoT2    = GNMales@pheno[,2],
                          EbvT1      = GNMales@ebv[, 1],
                          EbvT2      = GNMales@ebv[, 2],
                          TbvT1      = GNMales@gv[, 1],
                          TbvT2      = GNMales@gv[, 2]),
                   tibble(Generation = 0,
                          IId        = PNFemales@id,
                          FId        = NA,
                          MId        = NA,
                          Gender     = PNFemales@gender,
                          Program    = "PN2",
                          PhenoT1    = PNFemales@pheno[,1],
                          PhenoT2    = PNFemales@pheno[,2],
                          EbvT1      = PNFemales@ebv[, 1],
                          EbvT2      = PNFemales@ebv[, 2],
                          TbvT1      = PNFemales@gv[, 1],
                          TbvT2      = PNFemales@gv[, 2])
                          )

##############################################################################3
##############################################################################3
# ---- Program 0: GN  ----
for (Generation in 1:nGenerationEval) {

  # Mate
  SelCand = randCross2(females = GNFemales, males = GNMales,
                       nCrosses = GNFemales@nInd, nProgeny = 12)
  # meanG(SelCand)

  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }

  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)

  # Track pedigree
  PedEvalGN = rbind(PedEvalGN,
                    tibble(Generation = Generation,
                           IId        = SelCand@id,
                           FId        = SelCand@father,
                           MId        = SelCand@mother,
                           Gender     = SelCand@gender,
                           Program    = "GN",
                           PhenoT1    = SelCand@pheno[,1],
                           PhenoT2    = SelCand@pheno[,2],
                           EbvT1      = NA,
                           EbvT2      = NA,
                           TbvT1      = SelCand@gv[, 1],
                           TbvT2      = SelCand@gv[, 2]))
  
  # Estimate EBVs with blupf90
  setwd(paste0(homedir, "/GN/"))
  # Create pedigree file
  blupPed <- PedEvalGN[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(PedEvalGN[,c("IId", "FId", "MId")], "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  if (!file.exists("Blupf901.dat")) {
    write.table(PedEvalGN[,c("IId", "PhenoT1")], "Blupf901.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  } else {
    oldblupdat <- read.table("Blupf901.dat")
    colnames(oldblupdat) <- c("IId", "PhenoT1")
    newblupdat <- rbind(oldblupdat, PedEvalGN[,c("IId", "PhenoT1")])
    write.table(newblupdat, "Blupf901.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  }  
  
  # Trait 2
  # Create phenotype file
  if (!file.exists("Blupf902.dat")) {
    write.table(PedEvalGN[,c("IId", "PhenoT2")], "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  } else {
    oldblupdat <- read.table("Blupf902.dat")
    colnames(oldblupdat) <- c("IId", "PhenoT2")
    newblupdat <- rbind(oldblupdat, PedEvalGN[,c("IId", "PhenoT2")])
    write.table(newblupdat, "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  }
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[order(match(sol1$V2, PedEvalGN$IId)),]
  PedEvalGN$EbvT1 <- sol1$V3
  
  # Evaluate Trait 2
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[order(match(sol2$V2, PedEvalGN$IId)),]
  PedEvalGN$EbvT2 <- sol2$V3
  

  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEvalGN[PedEvalGN$IId %in% SelCand@id, c("EbvT1", "EbvT2")])

  # Select
  GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                        use = "ebv", trait = function(x) rowMeans(x))
  GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                        use = "ebv", trait = function(x) rowMeans(x))
  # Clean
  rm(SelCand)



##############################################################################3
##############################################################################3
# ---- Program 1: PN with 100% GNmales, PN does T1 ----

  if (Generation == 1) {
    PNFemales1 = PNFemales
  }

  # Mate
  PNFemales1ForPure = selectInd(pop = PNFemales1, nInd = PNFemales1@nInd * pPNFemalesPure,
                                use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  SelCand = randCross2(females = PNFemales1ForPure, males = GNMales,
                       nCrosses = PNFemales1ForPure@nInd, nProgeny = 14)
  # SelCand@nInd
  # meanG(PNFemales1); meanG(PNFemales1ForPure); meanG(GNMales); meanG(SelCand)

  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalPN1[(DataEvalPN1$Generation == Generation) & (DataEvalPN1$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalPN1[(DataEvalPN1$Generation == Generation) & (DataEvalPN1$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }


  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)
  SelCand@pheno[, 2] = NA


  # Track pedigree
  if (Generation == 1) {
    PedEvalPN1 = rbind(PedEvalPN1,
                       tibble(Generation = Generation,
                              IId        = SelCand@id,
                              FId        = SelCand@father,
                              MId        = SelCand@mother,
                              Gender     = SelCand@gender,
                              Program    = "PN1",
                              PhenoT1    = SelCand@pheno[,1],
                              PhenoT2    = SelCand@pheno[,2],
                              EbvT1      = NA,
                              EbvT2      = NA,
                              TbvT1      = SelCand@gv[, 1],
                              TbvT2      = SelCand@gv[, 2]))
  } else {
    PedEvalPN1 = rbind(PedEvalPN1,
                       tibble(Generation = Generation -1,
                              IId        = GNMales@id,
                              FId        = GNMales@father,
                              MId        = GNMales@mother,
                              Gender     = GNMales@gender,
                              Program    = "GN",
                              PhenoT1    = GNMales@pheno[,1],
                              PhenoT2    = GNMales@pheno[,2],
                              EbvT1      = GNMales@ebv[, 1],
                              EbvT2      = GNMales@ebv[, 2],
                              TbvT1      = GNMales@gv[, 1],
                              TbvT2      = GNMales@gv[, 2]),
                       tibble(Generation = Generation,
                              IId        = SelCand@id,
                              FId        = SelCand@father,
                              MId        = SelCand@mother,
                              Gender     = SelCand@gender,
                              Program    = "PN1",
                              PhenoT1    = SelCand@pheno[,1],
                              PhenoT2    = SelCand@pheno[,2],
                              EbvT1      = NA,
                              EbvT2      = NA,
                              TbvT1      = SelCand@gv[, 1],
                              TbvT2      = SelCand@gv[, 2]))
  }
  
  # Estimate EBVs with blupf90
  setwd(paste0(homedir, "/PN1/"))
  # Create pedigree file
  blupPed <- PedEvalPN1[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(PedEvalPN1[,c("IId", "FId", "MId")], "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  if (!file.exists("Blupf901.dat")) {
    write.table(PedEvalPN1[,c("IId", "PhenoT1")], "Blupf901.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  } else {
    oldblupdat <- read.table("Blupf901.dat")
    colnames(oldblupdat) <- c("IId", "PhenoT1")
    newblupdat <- rbind(oldblupdat, PedEvalPN1[,c("IId", "PhenoT1")])
    write.table(newblupdat, "Blupf901.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  }  
  
  # Trait 2
  # Create phenotype file
  if (!file.exists("Blupf902.dat")) {
    write.table(PedEvalPN1[,c("IId", "PhenoT2")], "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  } else {
    oldblupdat <- read.table("Blupf902.dat")
    colnames(oldblupdat) <- c("IId", "PhenoT2")
    newblupdat <- rbind(oldblupdat, PedEvalPN1[,c("IId", "PhenoT2")])
    write.table(newblupdat, "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  }
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[order(match(sol1$V2, PedEvalPN1$IId)),]
  PedEvalPN1$EbvT1 <- sol1$V3
  

  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEvalPN1[PedEvalPN1$IId %in% SelCand@id, c("EbvT1", "EbvT2")])

  EbvT2Males = merge(x = tibble(FId = unique(SelCand@father)),
                     y = tibble(IId = GNMales@id,
                                EbvT2 = GNMales@ebv[, 2]),
                     by.x = "FId", by.y = "IId")
  rownames(EbvT2Males) = EbvT2Males$FId
  EbvT2Females = merge(x = tibble(MId = unique(SelCand@mother)),
                       y = tibble(IId = PNFemales1ForPure@id,
                                  EbvT2 = PNFemales1ForPure@ebv[, 2]),
                       by.x = "MId", by.y = "IId")
  rownames(EbvT2Females) = EbvT2Females$MId
  SelCand@ebv[, 2] = 0.5 * (EbvT2Males[SelCand@father,   "EbvT2"] +
                              EbvT2Females[SelCand@mother, "EbvT2"])
  
  # Select
  PNMales1   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  PNFemales1 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))


  # Clean
  rm(SelCand)



##############################################################################3
##############################################################################3
# ---- Program 2: PN with 50% GNMales & 50% PNMales, PN does T1 ----
  if (Generation == 1) {
    PNFemales2 = PNFemales
    PNMales2 = GNMales
  }

  # Mate
  PNFemales2ForPure = selectInd(pop = PNFemales2, nInd = PNFemales2@nInd * pPNFemalesPure,
                                use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  SelCand = randCross2(females = PNFemales2ForPure, males = c(GNMales, PNMales2),
                       nCrosses = PNFemales2ForPure@nInd, nProgeny = 14)
  # SelCand@nInd
  # meanG(PNFemales2); meanG(PNFemales2ForPure); meanG(GNMales); meanG(PNMales2); meanG(SelCand)

  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalPN2[(DataEvalPN2$Generation == Generation) & (DataEvalPN2$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalPN2[(DataEvalPN2$Generation == Generation) & (DataEvalPN2$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }


  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)
  SelCand@pheno[, 2] = NA

  EbvT2Males = merge(x = tibble(FId = unique(SelCand@father)),
                     y = tibble(IId = c(GNMales, PNMales2)@id,
                                EbvT2 = c(GNMales, PNMales2)@ebv[, 2]),
                     by.x = "FId", by.y = "IId")
  EbvT2Males = EbvT2Males[!duplicated(EbvT2Males), ]
  rownames(EbvT2Males) = EbvT2Males$FId
  EbvT2Females = merge(x = tibble(MId = unique(SelCand@mother)),
                       y = tibble(IId = PNFemales2ForPure@id,
                                  EbvT2 = PNFemales2ForPure@ebv[, 2]),
                       by.x = "MId", by.y = "IId")
  EbvT2Females = EbvT2Females[!duplicated(EbvT2Females), ]
  rownames(EbvT2Females) = EbvT2Females$MId
  SelCand@ebv[, 2] = 0.5 * (EbvT2Males[SelCand@father,   "EbvT2"] +
                            EbvT2Females[SelCand@mother, "EbvT2"])

  
  # Track pedigree
  PedEvalPN2 = rbind(PedEvalPN2,
                     tibble(Generation = Generation - 1,
                            IId        = GNMales@id,
                            FId        = GNMales@father,
                            MId        = GNMales@mother,
                            Gender     = GNMales@gender,
                            Program    = "GN",
                            EbvT1      = GNMales@ebv[, 1],
                            EbvT2      = GNMales@ebv[, 2],
                            TbvT1      = GNMales@gv[, 1],
                            TbvT2      = GNMales@gv[, 2]),
                     tibble(Generation = Generation,
                            IId        = SelCand@id,
                            FId        = SelCand@father,
                            MId        = SelCand@mother,
                            Gender     = SelCand@gender,
                            Program    = "PN2",
                            EbvT1      = SelCand@ebv[, 1],
                            EbvT2      = SelCand@ebv[, 2],
                            TbvT1      = SelCand@gv[, 1],
                            TbvT2      = SelCand@gv[, 2]))
  
  # Estimate EBVs with blupf90
  setwd(paste0(homedir, "/PN1/"))
  # Create pedigree file
  blupPed <- PedEvalPN2[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(PedEvalPN2[,c("IId", "FId", "MId")], "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  if (!file.exists("Blupf901.dat")) {
    write.table(PedEvalPN2[,c("IId", "PhenoT1")], "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  } else {
    oldblupdat <- read.table("Blupf901.dat")
    colnames(oldblupdat) <- c("IId", "PhenoT1")
    newblupdat <- rbind(oldblupdat, PedEvalPN2[,c("IId", "PhenoT1")])
    write.table(newblupdat, "Blupf90.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  }  
  
  # Trait 2
  # Create phenotype file
  if (!file.exists("Blupf902.dat")) {
    write.table(PedEvalPN2[,c("IId", "PhenoT2")], "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  } else {
    oldblupdat <- read.table("Blupf902.dat")
    colnames(oldblupdat) <- c("IId", "PhenoT1")
    newblupdat <- rbind(oldblupdat, PedEvalPN2[,c("IId", "PhenoT2")])
    write.table(newblupdat, "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  }
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[order(match(sol1$V2, PedEvalPN2$IId)),]
  PedEvalPN2$EbvT1 <- sol1$V3
  
  # Evaluate Trait 2
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[order(match(sol2$V2, PedEvalPN2$IId)),]
  PedEvalPN2$EbvT2 <- sol2$V3
  
  # Select
  PNMales2   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  PNFemales2 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))

  # Clean
  rm(SelCand)
}

# Plot genetic means
ylim = range(rbind(DataEvalGN [, GenMeanCols],
                   DataEvalPN1[, GenMeanCols],
                   DataEvalPN2[, GenMeanCols]))
rbind(DataEvalGN,
      DataEvalPN1,
      DataEvalPN2) %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic mean") +
  ylim(ylim = ylim) +
  facet_wrap(~ Program + Gender, dir="h", nrow=3, ncol=2)

# Plot genetic variances
ylim = range(rbind(DataEvalGN [, GenVarCols],
                   DataEvalPN1[, GenVarCols],
                   DataEvalPN2[, GenVarCols]))
rbind(DataEvalGN,
      DataEvalPN1,
      DataEvalPN2) %>%
  gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic variance") +
  ylim(ylim = ylim) +
  facet_wrap(~ Program+Gender)


# ---- Partitioning the trend ----
PedEvalGN$EbvI = 0.5 * (PedEvalGN$EbvT1 + PedEvalGN$EbvT2)
PedEvalGN$EbvT1_s <- (PedEvalGN$EbvT1 - mean(PedEvalGN$EbvT1[PedEvalGN$Generation == 0])) / sd(PedEvalGN$EbvT1[PedEvalGN$Generation == 0])
PedEvalGN$EbvT2_s <- (PedEvalGN$EbvT2 - mean(PedEvalGN$EbvT2[PedEvalGN$Generation == 0])) / sd(PedEvalGN$EbvT2[PedEvalGN$Generation == 0])
PedEvalGN$EbvI_s <- (PedEvalGN$EbvI - mean(PedEvalGN$EbvI[PedEvalGN$Generation == 0])) / sd(PedEvalGN$EbvI[PedEvalGN$Generation == 0])



if (FALSE) {
  # Sum to zero constraint in the base population
  Sel = is.na(PedEvalGN$FId) & is.na(PedEvalGN$MId)
  Tmp = colMeans(PedEvalGN[Sel, c("EbvT1", "EbvT2", "EbvI")])
  PedEvalGN[, c("EbvT1", "EbvT2", "EbvI")] = PedEvalGN[, c("EbvT1", "EbvT2", "EbvI")] - Tmp
}

PartGN = AlphaPart(x = as.data.frame(PedEvalGN), sort = FALSE,
                 colId = "IId", colFid = "FId", colMid = "MId",
                 colPath = "Gender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
PartGNSummary = summary(object = PartGN, by = "Generation")
print(PartGNSummary)
p1 <- plot(PartGNSummary)
savePlot(p1, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_Gender", type = "jpeg", res=601, width=90, height=70, units="mm")
savePlot(x = p1, filename = "Partition_Gender", type = "jpeg")



PedEvalPN1$EbvI = 0.5 * (PedEvalPN1$EbvT1 + PedEvalPN1$EbvT2)
PedEvalPN1$EbvT1_s <- (PedEvalPN1$EbvT1 - mean(PedEvalPN1$EbvT1[PedEvalPN1$Generation == 0])) / sd(PedEvalPN1$EbvT1[PedEvalPN1$Generation == 0])
PedEvalPN1$EbvT2_s <- (PedEvalPN1$EbvT2 - mean(PedEvalPN1$EbvT2[PedEvalPN1$Generation == 0])) / sd(PedEvalPN1$EbvT2[PedEvalPN1$Generation == 0])
PedEvalPN1$EbvI_s <- (PedEvalPN1$EbvI - mean(PedEvalPN1$EbvI[PedEvalPN1$Generation == 0])) / sd(PedEvalPN1$EbvI[PedEvalPN1$Generation == 0])



PartPN1 = AlphaPart(x = as.data.frame(PedEvalPN1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "Gender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
PartPN1Summary = summary(object = PartPN1, by = "Generation")
print(PartPN1Summary)
p2 <- plot(PartPN1Summary)
savePlot(p2, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_ProgramP1", type = "jpeg", res=601, width=90, height=70, units="mm")





PedEvalPN2$EbvI = 0.5 * (PedEvalPN2$EbvT1 + PedEvalPN2$EbvT2)
PedEvalPN2$EbvT1_s <- (PedEvalPN2$EbvT1 - mean(PedEvalPN2$EbvT1[PedEvalPN2$Generation == 0])) / sd(PedEvalPN2$EbvT1[PedEvalPN2$Generation == 0])
PedEvalPN2$EbvT2_s <- (PedEvalPN2$EbvT2 - mean(PedEvalPN2$EbvT2[PedEvalPN2$Generation == 0])) / sd(PedEvalPN2$EbvT2[PedEvalPN2$Generation == 0])
PedEvalPN2$EbvI_s <- (PedEvalPN2$EbvI - mean(PedEvalPN2$EbvI[PedEvalPN2$Generation == 0])) / sd(PedEvalPN2$EbvI[PedEvalPN2$Generation == 0])



PartPN2 = AlphaPart(x = as.data.frame(PedEvalPN2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "Program", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
PartPN2Summary = summary(object = PartPN2, by = "Generation")
print(PartPN2Summary)
p3 <- plot(PartPN2Summary)
savePlot(p3, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_ProgramP2", type = "jpeg", res=601, width=90, height=70, units="mm")

PedEvalPN2$ProgramGender = paste(PedEvalPN2$Program, PedEvalPN2$Gender, sep = "-")
PartPN2 = AlphaPart(x = as.data.frame(PedEvalPN2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
PartPN2Summary = summary(object = PartPN2, by = "Generation")
print(PartPN2Summary)
p4 <- plot(PartPN2Summary)
savePlot(p4, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_ProgramGenderP2", type = "jpeg", res=601, width=90, height=70, units="mm")

##############################################################
#tgv
# ---- Partitioning the trend ----

PedEvalGN$TbvI = 0.5 * (PedEvalGN$TbvT1 + PedEvalGN$TbvT2)

#standardizirar
PedEvalGN$TbvT1_s <- (PedEvalGN$TbvT1 - mean(PedEvalGN$TbvT1[PedEvalGN$Generation == 0])) / sd(PedEvalGN$TbvT1[PedEvalGN$Generation == 0])
PedEvalGN$TbvT2_s <- (PedEvalGN$TbvT2 - mean(PedEvalGN$TbvT2[PedEvalGN$Generation == 0])) / sd(PedEvalGN$TbvT2[PedEvalGN$Generation == 0])
PedEvalGN$TbvI_s <- (PedEvalGN$TbvI - mean(PedEvalGN$TbvI[PedEvalGN$Generation == 0])) / sd(PedEvalGN$TbvI[PedEvalGN$Generation == 0])

if (FALSE) {
  # Sum to zero constraint in the base population
  Sel = is.na(PedEvalGN$FId) & is.na(PedEvalGN$MId)
  Tmp = colMeans(PedEvalGN[Sel, c("TbvT1", "TbvT2", "TbvI")])
  PedEvalGN[, c("TbvT1", "TbvT2", "TbvI")] = PedEvalGN[, c("TbvT1", "TbvT2", "TbvI")] - Tmp
}
if (FALSE) {
  # Sum to zero constraint in the base population
  Sel = is.na(PedEvalGN$FId) & is.na(PedEvalGN$MId)
  Tmp = colMeans(PedEvalGN[Sel, c("TbvT1", "TbvT2", "TbvI")])
  PedEvalGN[, c("TbvT1_s", "TbvT2_s", "TbvI_s")] = PedEvalGN[, c("TbvT1_s", "TbvT2_s", "TbvI_s")] - Tmp
}

PartGN = AlphaPart(x = as.data.frame(PedEvalGN), sort = FALSE,
                 colId = "IId", colFid = "FId", colMid = "MId",
                 colPath = "Gender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
PartGNSummary = summary(object = PartGN, by = "Generation")

print(PartGNSummary)
p1 <- plot(PartGNSummary)
savePlot(p1, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_Gender", type = "jpeg", res=601, width=90, height=70, units="mm")
avePlot(x = p1, filename = "Partition_Gender", type = "jpeg")


PedEvalPN1$TbvI = 0.5 * (PedEvalPN1$TbvT1 + PedEvalPN1$TbvT2)
PedEvalPN1$TbvT1_s <- (PedEvalPN1$TbvT1 - mean(PedEvalPN1$TbvT1[PedEvalPN1$Generation == 0])) / sd(PedEvalPN1$TbvT1[PedEvalPN1$Generation == 0])
PedEvalPN1$TbvT2_s <- (PedEvalPN1$TbvT2 - mean(PedEvalPN1$TbvT2[PedEvalPN1$Generation == 0])) / sd(PedEvalPN1$TbvT2[PedEvalPN1$Generation == 0])
PedEvalPN1$TbvI_s <- (PedEvalPN1$TbvI - mean(PedEvalPN1$TbvI[PedEvalPN1$Generation == 0])) / sd(PedEvalPN1$TbvI[PedEvalPN1$Generation == 0])

PartPN1 = AlphaPart(x = as.data.frame(PedEvalPN1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "Gender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
PartPN1Summary = summary(object = PartPN1, by = "Generation")
print(PartPN1Summary)
p2 <- plot(PartPN1Summary)
savePlot(p2, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_ProgramP1", type = "jpeg", res=601, width=90, height=70, units="mm")


PedEvalPN2$TbvI = 0.5 * (PedEvalPN2$TbvT1 + PedEvalPN2$TbvT2)
PedEvalPN2$TbvT1_s <- (PedEvalPN2$TbvT1 - mean(PedEvalPN2$TbvT1[PedEvalPN2$Generation == 0])) / sd(PedEvalPN2$TbvT1[PedEvalPN2$Generation == 0])
PedEvalPN2$TbvT2_s <- (PedEvalPN2$TbvT2 - mean(PedEvalPN2$TbvT2[PedEvalPN2$Generation == 0])) / sd(PedEvalPN2$TbvT2[PedEvalPN2$Generation == 0])
PedEvalPN2$TbvI_s <- (PedEvalPN2$TbvI - mean(PedEvalPN2$TbvI[PedEvalPN2$Generation == 0])) / sd(PedEvalPN2$TbvI[PedEvalPN2$Generation == 0])

PartPN2 = AlphaPart(x = as.data.frame(PedEvalPN2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "Program", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
PartPN2Summary = summary(object = PartPN2, by = "Generation")
print(PartPN2Summary)
p3 <- plot(PartPN2Summary)
savePlot(p3, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_ProgramP2", type = "jpeg", res=601, width=90, height=70, units="mm")

PedEvalPN2$ProgramGender = paste(PedEvalPN2$Program, PedEvalPN2$Gender, sep = "-")
PartPN2 = AlphaPart(x = as.data.frame(PedEvalPN2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
PartPN2Summary = summary(object = PartPN2, by = "Generation")
print(PartPN2Summary)
p4 <- plot(PartPN2Summary)
savePlot(p4, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_ProgramGenderP2", type = "jpeg", res=601, width=90, height=70, units="mm")


#Plot true / estimated breeding values by sex
