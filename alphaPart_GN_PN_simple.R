# 
# ---- Environment ----
p1Interaction <- p1
p1gInteraction <- p1g

rm(list = ls())

library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")
library(pedigree)

setwd("/home/jana/Documents/SimulationAlphaPart/")

homedir = "/home/jana/Documents/SimulationAlphaPart/"
# ---- General parameters ----

nGNMales   =  10
nGNFemales = 50

nPNMales       =  10
nPNFemales     = 500
pPNFemalesPure = 0.15

nGenerationBurn = 20
nGenerationEval = 20

GenMeanCols = c("GenMeanT1", "GenMeanT2", "GenMeanI")
GenVarCols  = c("GenVarT1",  "GenVarT2",  "GenVarI")

# ---- Base population genomes ----

founderPop = runMacs(nInd = (nGNMales + nGNFemales)*2,
                     nChr = 10,
                     segSites = 100,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")

###################################################################################
###################################################################################
# ---- Simulation/Base population parameters ----

SP = SimParam$new(founderPop)
# VarA = matrix(data = c(1.0, 0.1, 0.1, 1.0), nrow = 2); cov2cor(VarA)
VarA = matrix(data = c(1.0, 0.0, 0.0, 1.0), nrow = 2); cov2cor(VarA)
VarE = matrix(data = c(3.0, 0.0, 0.0, 3.0), nrow = 2); cov2cor(VarE)

VarP = VarA + VarE; diag(VarA) / diag(VarP)
SP$addTraitA(nQtlPerChr = 100, mean = c(0, 0), var = diag(VarA), cor = cov2cor(VarA))
# SP$addSnpChip(nSnpPerChr = 1000)
SP$setGender(gender = "yes_sys")

# ---- Base GN population ----

GN = newPop(founderPop)
# GN@gender[1:60] <- "F"
# GN@gender[sample(1:60, 10)] <- "M"
table(GN@gender)
BaseGNMales   = GN[GN@gender == "M"]
BaseGNFemales = GN[GN@gender == "F"]
rm(GN)

###################################################################################
###################################################################################

# ---- GN burn-in ----

DataBurn = tibble(Generation = rep(1:nGenerationBurn, each=2), Gender = rep(c("F", "M"), nGenerationBurn),
                  GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                  GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
for (Generation in 1:nGenerationBurn) {
  # Generation = 1
  
  # Mate
  SelCand = randCross2(females = BaseGNFemales, males = BaseGNMales,
                       nCrosses = BaseGNFemales@nInd, nProgeny = 12)
  
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
  mean(SelCand@ebv[SelCand@gender == "M", 1])
  mean(SelCand@ebv[SelCand@gender == "M", 2])
  mean(SelCand@ebv[SelCand@gender == "F", 1])
  mean(SelCand@ebv[SelCand@gender == "F", 2])
  
  # Select
  BaseGNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                            use = "ebv", trait = function(x) rowMeans(x))
  mean(BaseGNMales@ebv[,1])
  mean(BaseGNMales@ebv[,2])
  BaseGNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                            use = "ebv", trait = function(x) rowMeans(x))
  mean(BaseGNFemales@ebv[,1])
  mean(BaseGNFemales@ebv[,2])
  
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

GNInd = c(BaseGNFemales@id, BaseGNMales@id)
SelCand = SelCand[!(SelCand@id %in% GNInd)]; SelCand@nInd
mean(SelCand@ebv[SelCand@gender == "M", 1])
mean(SelCand@ebv[SelCand@gender == "M", 2])
mean(SelCand@ebv[SelCand@gender == "F", 1])
mean(SelCand@ebv[SelCand@gender == "F", 2])
# Cheat here - consider all animals are females
SelCand@gender[] = "F"
BasePNFemales = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                          use = "ebv", trait = function(x) rowMeans(x))
mean(BasePNFemales@ebv[,1])
mean(BasePNFemales@ebv[,2])
mean(BaseGNFemales@ebv[,1])
mean(BaseGNFemales@ebv[,2])
mean(BaseGNMales@ebv[,1])
mean(BaseGNMales@ebv[,2])
plot(density(BasePNFemales@ebv))
plot(density(BaseGNFemales@ebv))
plot(density(BaseGNMales@ebv))
# meanG(GNMales); meanG(GNFemales); meanG(PNFemales)


# ---- PN1 ----

# PedEval = rbind(tibble(Generation = 0,
#                        IId        = BaseGNMales@id,
#                        FId        = NA,
#                        MId        = NA ,
#                        Gender     = BaseGNMales@gender,
#                        Program    = "GN",
#                        PhenoT1    = NA,
#                        PhenoT2    = NA,
#                        EbvT1      = BaseGNMales@ebv[, 1],
#                        EbvT2      = BaseGNMales@ebv[, 2],
#                        TbvT1      = BaseGNMales@gv[, 1],
#                        TbvT2      = BaseGNMales@gv[, 2]),
#                 tibble(Generation = 0,
#                        IId        = BaseGNFemales@id,
#                        FId        = NA,
#                        MId        = NA,
#                        Gender     = BaseGNFemales@gender,
#                        Program    = "GN",
#                        PhenoT1    = NA,
#                        PhenoT2    = NA,
#                        EbvT1      = BaseGNFemales@ebv[, 1],
#                        EbvT2      = BaseGNFemales@ebv[, 2],
#                        TbvT1      = BaseGNFemales@gv[, 1],
#                        TbvT2      = BaseGNFemales@gv[, 2]),
#                 tibble(Generation = 1,
#                        IId        = BasePNFemales@id,
#                        FId        = NA,
#                        MId        = NA,
#                        Gender     = BasePNFemales@gender,
#                        Program    = "PN1",
#                        PhenoT1    = NA,
#                        PhenoT2    = NA,
#                        EbvT1      = BasePNFemales@ebv[, 1],
#                        EbvT2      = BasePNFemales@ebv[, 2],
#                        TbvT1      = BasePNFemales@gv[, 1],
#                        TbvT2      = BasePNFemales@gv[, 2])
# )


DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                    Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                    GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                    GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
DataEvalGN$Program = "GN"
DataEvalPN1 = DataEvalGN
DataEvalPN1$Program = "PN1"
DataEvalPN2 = DataEvalGN
DataEvalPN2$Program = "PN2"

accuraciesPN1 <- data.frame(Program = NA, Generation = NA, Trait = NA, Cor = NA)
accuraciesPN2 <- data.frame(Program = NA, Generation = NA, Trait = NA, Cor = NA)
##############################################################################3
##############################################################################3
# ---- Program PN1  ----
BasePNFemales@pheno[,2] <- NA
for (Generation in 1:nGenerationEval) {  ##nGenerationEval) {
  if (Generation == 1) {
    GNFemales = BaseGNFemales
    GNMales = BaseGNMales
  }
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
  if (Generation == 1) {
    PedEval = rbind(tibble(Generation = Generation,
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
    
  } else {
  PedEval = rbind(PedEval,
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
  }
  
  # Estimate EBVs with blupf90
  setwd(paste0(homedir, "/GN/"))
  if (Generation == 1) {
    system("rm *ped *dat")
  }
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  blup1dat <- PedEval[,c("IId", "PhenoT1")]
  blup1dat$Mean <- 1
  write.table(blup1dat, "Blupf901.dat",
                          quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                          na = "0", append = file.exists("Blupf901.dat"))
  # write.table(PedEval[PedEval$Program == "GN",c("IId", "PhenoT1", "Program", "Generation")], "Blupf901.dat",
  #             quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
  #             na = "0", append = file.exists("Blupf901.dat"))
  
  
  # Trait 2
  # Create phenotype file
  blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2")]
  blup2dat$Mean <- 1
  write.table(blup2dat, "Blupf902.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf902.dat"))
  # write.table(PedEval[,c("IId", "PhenoT2", "Program", "Generation")], "Blupf902.dat",
  #             quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
  #             na = "0", append = file.exists("Blupf902.dat"))
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId,]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3
  
  system("rm renf90.par solutions")
  # Evaluate Trait 2
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId,]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  accuraciesPN1 <- rbind(accuraciesPN1, c("GN", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN1 <- rbind(accuraciesPN1, c("GN", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
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
    PNFemales1 = BasePNFemales
  }
  PNFemales1@pheno[,2] <- NA
  
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
  PedEval = rbind(PedEval, 
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
  
  
  # Estimate EBVs with blupf90
  setwd(paste0(homedir, "/PN1/"))
  if (Generation == 1) {
    system("rm *ped *dat")
  }
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  blup1dat <- PedEval[,c("IId", "PhenoT1")]
  blup1dat$Mean <- 1
  write.table(blup1dat, "Blupf901.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf901.dat"))
  # write.table(PedEval[PedEval$Program == "GN",c("IId", "PhenoT1", "Program", "Generation")], "Blupf901.dat",
  #             quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
  #             na = "0", append = file.exists("Blupf901.dat"))
  
  
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId, ]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3  
  
  # Trait 2
  # Create phenotype file
  blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2")]
  blup2dat$Mean <- 1
  write.table(blup2dat, "Blupf902.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf902.dat"))
  # write.table(PedEval[,c("IId", "PhenoT2", "Program", "Generation")], "Blupf902.dat",
  #             quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
  #             na = "0", append = file.exists("Blupf902.dat"))
  
  system("rm renf90.par solutions")
  # Evaluate Trait 2
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId, ]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  
  # Select
  # PNMales1   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
  #                        use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  PNFemales1 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  
  accuraciesPN1 <- rbind(accuraciesPN1, c("PN1", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN1 <- rbind(accuraciesPN1, c("PN1", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  
  # Clean
  rm(SelCand)
}

PedEval1 = PedEval
rm(PedEval)
DataEvalGN1 = DataEvalGN
rm(DataEvalGN)
write.table(PedEval1,"PedEval_NewModel_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_Intercept_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_Gender_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_Intercept_missingT1_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_NoMean_missingT1_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_NoMean_ChangedGen_missingT1_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_NoMean_NoBaseInPed_missingT1_0.25.csv", quote=FALSE, row.names=FALSE)
write.table(PedEval1,"PedEval_Intercept_NoBaseInPed_missingT1_0.25.csv", quote=FALSE, row.names=FALSE)

# ---- Partitioning the trend PN1 ----
# 
# PedEval1$EbvI = 0.5 * (PedEval1$EbvT1 + PedEval1$EbvT2)
# PedEval1$EbvT1_s[PedEval1$Gender == "M"] <- (PedEval1$EbvT1[PedEval1$Gender == "M"] - mean(PedEval1$EbvT1[PedEval1$Generation == 2 & PedEval1$Gender == "M"])) / sd(PedEval1$TbvT1[PedEval1$Generation == 2 & PedEval1$Gender == "M"])
# PedEval1$EbvT2_s[PedEval1$Gender == "M"] <- (PedEval1$EbvT2[PedEval1$Gender == "M"] - mean(PedEval1$EbvT2[PedEval1$Generation == 2 & PedEval1$Gender == "M"])) / sd(PedEval1$TbvT2[PedEval1$Generation == 2 & PedEval1$Gender == "M"])
# PedEval1$EbvT1_s[PedEval1$Gender == "F"] <- (PedEval1$EbvT1[PedEval1$Gender == "F"] - mean(PedEval1$EbvT1[PedEval1$Generation == 2 & PedEval1$Gender == "F"])) / sd(PedEval1$TbvT1[PedEval1$Generation == 2 & PedEval1$Gender == "F"])
# PedEval1$EbvT2_s[PedEval1$Gender == "F"] <- (PedEval1$EbvT2[PedEval1$Gender == "F"] - mean(PedEval1$EbvT2[PedEval1$Generation == 2 & PedEval1$Gender == "F"])) / sd(PedEval1$TbvT2[PedEval1$Generation == 2 & PedEval1$Gender == "F"])
# #PedEval1$EbvI_s <- (PedEval1$EbvI - mean(PedEval1$EbvI[PedEval1$Generation == 0])) / sd(PedEval1$EbvI[PedEval1$Generation == 0])
# PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)
# 
# 
# PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
# PedEval1$TbvT1_s[PedEval1$Program == "GN"] <- (PedEval1$TbvT1[PedEval1$Program == "GN"] - mean(PedEval1$TbvT1[PedEval1$Generation == 2 & PedEval1$Program == "GN"])) / sd(PedEval1$TbvT1[PedEval1$Generation == 2 & PedEval1$Program == "GN"])
# PedEval1$TbvT2_s[PedEval1$Program == "GN"] <- (PedEval1$TbvT2[PedEval1$Program == "GN"] - mean(PedEval1$TbvT2[PedEval1$Generation == 2 & PedEval1$Program == "GN"])) / sd(PedEval1$TbvT2[PedEval1$Generation == 2 & PedEval1$Program == "GN"])
# PedEval1$TbvT1_s[PedEval1$Program == "PN1"] <- (PedEval1$TbvT1[PedEval1$Program == "PN1"] - mean(PedEval1$TbvT1[PedEval1$Generation == 2 & PedEval1$Program == "PN1"])) / sd(PedEval1$TbvT1[PedEval1$Generation == 2 & PedEval1$Program == "PN1"])
# PedEval1$TbvT2_s[PedEval1$Program == "PN1"] <- (PedEval1$TbvT2[PedEval1$Program == "PN1"] - mean(PedEval1$TbvT2[PedEval1$Generation == 2 & PedEval1$Program == "PN1"])) / sd(PedEval1$TbvT2[PedEval1$Generation == 2 & PedEval1$Program == "PN1"])

PedEvalCOPY <- PedEval1


Ped <- PedEval1[,c("IId", "MId", "FId")]
library(pedigree)
mothers <- unique(Ped$MId)
fathers <- unique(Ped$FId)
who <- c(Ped$IId %in% c(mothers, fathers))
who <- c(Ped$IId %in% c(PedEval1$IId[PedEval1$Generation == max(PedEval1$Generation)]))
Trimed <- trimPed(ped = Ped, data = who)
PedEval1 <- PedEval1[Trimed,]

PedEval1 <- read.csv("PN1/PedEval_Intercept_0.25.csv", header=TRUE)
PedEval1 <- read.table("PN1/PedEval_NewModel_0.25.csv", header=TRUE)
PedEval1$ProgramGender <- paste0(PedEval1$Program, PedEval1$Gender)
#PedEvalM <- PedEval1[PedEval1$Gender == "M" & PedEval1$Program == "GN",]
InterceptSD <- aggregate(PedEvalM$EbvT1 ~ PedEvalM$Generation, FUN="sd")
InterceptSD <- aggregate(PedEval1$EbvT1 ~ PedEval1$Generation + PedEval1$ProgramGender, FUN="sd")
InterceptSD <- aggregate(PedEval1$TbvT1 ~ PedEval1$Generation + PedEval1$ProgramGender, FUN="sd")
InterceptSD <- aggregate(PedEval1$TbvT1 ~ PedEval1$Generation + PedEval1$ProgramGender, FUN="var")
InterceptSD <- aggregate(PedEval1$TbvT1 ~ PedEval1$Generation + PedEval1$ProgramGender, FUN="mean")
InterceptSD$Model <- "Int"
colnames(InterceptSD) <- c("Generation", "ProgramGender", "TbvT1", "Model")
ggplot(InterceptSD, aes(x=Generation, y=TbvT1, group=ProgramGender, colour=ProgramGender)) + geom_line()


PedEval1 <- read.table("PedEval_NewModel_0.25.csv", header=TRUE)
#PedEvalM <- PedEval1[PedEval1$Gender == "M" & PedEval1$Program == "GN",]
PedEval1$ProgramGender <- paste0(PedEval1$Program, PedEval1$Gender)
NewModelSD <- aggregate(PedEvalM$EbvT1 ~ PedEvalM$Generation, FUN="sd")
NewModelSD <- aggregate(PedEval1$EbvT1 ~ PedEval1$Generation + PedEval1$ProgramGender, FUN="sd")
NewModelSD <- aggregate(PedEval1$TbvT1 ~ PedEval1$Generation + PedEval1$ProgramGender, FUN="sd")
NewModelSD$Model <- "New"

SD <- rbind(InterceptSD, NewModelSD)
colnames(SD) <- c("Generation", "SD", "Model")
colnames(SD) <- c("Generation", "ProgramGender", "SD", "Model")
ggplot(data = SD, aes(x = Generation, y = SD, group=Model, colour=Model)) + geom_line() + facet_grid(.~ProgramGender)

PedEval1 <- read.table("PedEval_Intercept_NoBaseInPed_missingT1_0.25.csv", header=TRUE)
table(PedEval1$Generation)
PedEval1$EbvI = 0.5 * (PedEval1$EbvT1 + PedEval1$EbvT2)
PedEval1$EbvT1_s <- (PedEval1$EbvT1 - mean(PedEval1$EbvT1[PedEval1$Generation == 1])) / sd(PedEval1$TbvT1[PedEval1$Generation == 1])
PedEval1$EbvT2_s <- (PedEval1$EbvT2 - mean(PedEval1$EbvT2[PedEval1$Generation == 1])) / sd(PedEval1$TbvT2[PedEval1$Generation == 1])
#PedEval1$EbvI_s <- (PedEval1$EbvI - mean(PedEval1$EbvI[PedEval1$Generation == 0])) / sd(PedEval1$EbvI[PedEval1$Generation == 0])
PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)


PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
PedEval1$TbvT1_s <- (PedEval1$TbvT1 - mean(PedEval1$TbvT1[PedEval1$Generation == 1])) / sd(PedEval1$TbvT1[PedEval1$Generation == 1])
PedEval1$TbvT2_s <- (PedEval1$TbvT2 - mean(PedEval1$TbvT2[PedEval1$Generation == 1])) / sd(PedEval1$TbvT2[PedEval1$Generation == 1])
#PedEval1$TbvI_s <- (PedEval1$TbvI - mean(PedEval1$TbvI[PedEval1$Generation == 0])) / sd(PedEval1$TbvI[PedEval1$Generation == 0])
PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)


PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1", "EbvT2", "EbvI"))
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1_s", "TbvT1_s"))
Part1g = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                   colId = "IId", colFid = "FId", colMid = "MId",
                   colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
Part1g = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                   colId = "IId", colFid = "FId", colMid = "MId",
                   colPath = "ProgramGender", colAGV = c("TbvT1", "TbvT2", "TbvI"))


Part1Summary = summary(object = Part1, by = "Generation")
Part1gSummary = summary(object = Part1g, by = "Generation")
p1 <- plot(Part1Summary)
p1g <- plot(Part1gSummary)



Part1Summary$EbvT1_s$rel

accuraciesPN1$Cor <- as.numeric(accuraciesPN1$Cor)
ggplot(accuraciesPN1, aes(x=Generation, y=Cor, group=Program, colour=Program, shape=Trait)) +  geom_point()


PedEvalIntercept <- PedEval1

head(PedEval1)
gainAe <- summarySE(data=PedEval1, measurevar = "EbvT1", groupvars = c("Generation", "Program", "Gender"))
gainAe$BV <- "EBV"
colnames(gainAe)[5] <- "value"
gainAt <- summarySE(data=PedEval1, measurevar = "TbvT1", groupvars = c("Generation", "Program", "Gender"))
gainAt$BV <- "TBV"
colnames(gainAt)[5] <- "value"
gainA <- rbind(gainAe, gainAt)

t1 <-ggplot(data = gainA, aes(x=Generation, y=value, group=BV, colour=BV)) + geom_line() + facet_grid(~ Program  + Gender)

head(PedEval1)
gainAe <- summarySE(data=PedEval1, measurevar = "EbvT2", groupvars = c("Generation", "Program", "Gender"))
gainAe$BV <- "EBV"
colnames(gainAe)[5] <- "value"
gainAt <- summarySE(data=PedEval1, measurevar = "TbvT2", groupvars = c("Generation", "Program", "Gender"))
gainAt$BV <- "TBV"
colnames(gainAt)[5] <- "value"
gainA <- rbind(gainAe, gainAt)

t2 <- ggplot(data = gainA, aes(x=Generation, y=value, group=BV, colour=BV)) + geom_line() + facet_grid(~ Program + Gender)
t2 <- ggplot(data = gainA, aes(x=Generation, y=value, group=Gender, colour=Gender)) + geom_line() + facet_grid(~ BV + Program)



write.table(PedEval1, paste0("PedEval1.csv"), quote=FALSE, row.names=FALSE)
write.table(DataEvalGN1, paste0("DataEval1.csv"), quote=FALSE, row.names=FALSE)
write.table(accuraciesPN1, paste0("Accuracies1.csv"), quote=FALSE, row.names=FALSE)


summarySE(Part1$EbvT1_s, measurevar = "EbvT1_s_w", groupvars = c("Gender", "Program"))
summary(Part1$EbvT1_s$EbvT1_s_w[Part1$EbvT1_s$Program == "PN1"])
length(Part1$EbvT1_s$EbvT1_s_w[Part1$EbvT1_s$Program == "PN1"])
length(Part1$EbvT1_s$EbvT1_s_w[Part1$EbvT1_s$Program == "GN"])
ggplot(data = Part1$EbvT1_s[Part1$EbvT1_s$Program == "PN1",], aes(x = EbvT1_s_w)) + geom_density(aes(fill = Gender), alpha = 0.4) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+ 
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))
ggplot(data = Part1$EbvT1_s[Part1$EbvT1_s$Program == "GN",], aes(x = EbvT1_s_w)) + geom_density(aes(fill = Gender), alpha = 0.4) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+ 
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))
##############################################################################3
##############################################################################3
# ---- PN2 ----
DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                    Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                    GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                    GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
DataEvalGN$Program = "GN"

PedEval = rbind(tibble(Generation = 0,
                       IId        = BaseGNMales@id,
                       FId        = NA,
                       MId        = NA,
                       Gender     = BaseGNMales@gender,
                       Program    = "GN",
                       PhenoT1    = BaseGNMales@pheno[,1],
                       PhenoT2    = BaseGNMales@pheno[,2],
                       EbvT1      = BaseGNMales@ebv[, 1],
                       EbvT2      = BaseGNMales@ebv[, 2],
                       TbvT1      = BaseGNMales@gv[, 1],
                       TbvT2      = BaseGNMales@gv[, 2]),
                tibble(Generation = 0,
                       IId        = BaseGNFemales@id,
                       FId        = NA,
                       MId        = NA,
                       Gender     = BaseGNFemales@gender,
                       Program    = "GN",
                       PhenoT1    = BaseGNFemales@pheno[,1],
                       PhenoT2    = BaseGNFemales@pheno[,2],
                       EbvT1      = BaseGNFemales@ebv[, 1],
                       EbvT2      = BaseGNFemales@ebv[, 2],
                       TbvT1      = BaseGNFemales@gv[, 1],
                       TbvT2      = BaseGNFemales@gv[, 2]),
                tibble(Generation = 0,
                       IId        = BasePNFemales@id,
                       FId        = NA,
                       MId        = NA,
                       Gender     = BasePNFemales@gender,
                       Program    = "PN2",
                       PhenoT1    = BasePNFemales@pheno[,1],
                       PhenoT2    = BasePNFemales@pheno[,2],
                       EbvT1      = BasePNFemales@ebv[, 1],
                       EbvT2      = BasePNFemales@ebv[, 2],
                       TbvT1      = BasePNFemales@gv[, 1],
                       TbvT2      = BasePNFemales@gv[, 2])
)

for (Generation in 1:nGenerationEval) {
  if (Generation == 1) {
    GNFemales = BaseGNFemales
    GNMales = BaseGNMales
  }
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
  PedEval = rbind(PedEval,
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
  if (Generation == 1) {
    system("rm *ped *dat")
  }
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  write.table(PedEval[,c("IId", "PhenoT1", "Program", "Generation")], "Blupf901.dat", 
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", 
              na = "0", append = file.exists("Blupf901.dat"))
  
  # Trait 2
  # Create phenotype file
  write.table(PedEval[PedEval$Program == "GN",c("IId", "PhenoT2", "Program", "Generation")], "Blupf902.dat", 
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", 
              na = "0", append = file.exists("Blupf902.dat"))
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId, ]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3
  
  # Evaluate Trait 2
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId, ]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  accuraciesPN2 <- rbind(accuraciesPN2, c("GN", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN2 <- rbind(accuraciesPN2, c("GN", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  # Select
  GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                        use = "ebv", trait = function(x) rowMeans(x))
  GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                        use = "ebv", trait = function(x) rowMeans(x))
  # Clean
  rm(SelCand)
  
  # ---- Program 2: PN with 50% GNMales & 50% PNMales, PN does T1 ----
  if (Generation == 1) {
    PNFemales2 = BasePNFemales
    PNMales2 = BaseGNMales
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
  
  
  # Track pedigree
  PedEval = rbind(PedEval, 
                  tibble(Generation = Generation,
                         IId        = SelCand@id,
                         FId        = SelCand@father,
                         MId        = SelCand@mother,
                         Gender     = SelCand@gender,
                         Program    = "PN2",
                         PhenoT1    = SelCand@pheno[,1],
                         PhenoT2    = SelCand@pheno[,2],
                         EbvT1      = NA,
                         EbvT2      = NA,
                         TbvT1      = SelCand@gv[, 1],
                         TbvT2      = SelCand@gv[, 2]))
  
  # Estimate EBVs with blupf90
  setwd(paste0(homedir, "/PN2/"))
  if (Generation == 1) {
    system("rm *ped *dat")
  }
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file
  # Trait 1
  write.table(PedEval[,c("IId", "PhenoT1", "Program", "Generation")], "Blupf901.dat", 
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", 
              na = "0", append = file.exists("Blupf901.dat"))
  
  # Evaluate Trait 1
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId, ]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3
  
  
  # Create phenotype file
  # Trait 2
  write.table(PedEval[PedEval$Program == "GN",c("IId", "PhenoT2", "Program", "Generation")], "Blupf902.dat", 
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", 
              na = "0", append = file.exists("Blupf902.dat"))
  
  
  # Evaluate Trait 2
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId, ]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  accuraciesPN2 <- rbind(accuraciesPN2, c("PN2", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN2 <- rbind(accuraciesPN2, c("PN2", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  # Select
  PNMales2   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  PNFemales2 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  
  
  # Clean
  rm(SelCand)
}

PedEval2 = PedEval
rm(PedEval)
DataEvalGN2 <- DataEvalGN
rm(DataEvalGN)

accuraciesPN2$Cor <- as.numeric(accuraciesPN2$Cor)
ggplot(accuraciesPN2, aes(x=Generation, y=Cor, group=Program, colour=Program, shape=Trait)) +  geom_point()


write.table(PedEval2, paste0("PedEval2.csv"), quote=FALSE, row.names=FALSE)
write.table(DataEvalGN2, paste0("DataEval2.csv"), quote=FALSE, row.names=FALSE)
write.table(accuraciesPN2, paste0("Accuracies2.csv"), quote=FALSE, row.names=FALSE)
###################################################################################################
####################################################################################################

# Plot genetic means
# PN1
ylim = range(rbind(DataEvalGN1 [, GenMeanCols],
                   DataEvalPN1[, GenMeanCols]))
rbind(DataEvalGN1,
      DataEvalPN1) %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic mean") +
  ylim(ylim = ylim) +
  facet_wrap(~ Program + Gender, dir="h", nrow=3, ncol=2)

# Plot genetic variances
ylim = range(rbind(DataEvalGN1 [, GenVarCols],
                   DataEvalPN1[, GenVarCols]))

rbind(DataEvalGN1,
      DataEvalPN1) %>%
  gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic variance") +
  ylim(ylim = ylim) +
  facet_wrap(~ Program+Gender)

# Plot genetic means
# PN2
ylim = range(rbind(DataEvalGN2 [, GenMeanCols],
                   DataEvalPN2[, GenMeanCols]))
rbind(DataEvalGN2,
      DataEvalPN2) %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic mean") +
  ylim(ylim = ylim) +
  facet_wrap(~ Program + Gender, dir="h", nrow=3, ncol=2)

# Plot genetic variances
ylim = range(rbind(DataEvalGN2 [, GenVarCols],
                   DataEvalPN2[, GenVarCols]))

rbind(DataEvalGN2,
      DataEvalPN2) %>%
  gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic variance") +
  ylim(ylim = ylim) +
  facet_wrap(~ Program+Gender)


PedEval1 <- PedEval1[PedEval1$Generation > 10,]
# ---- Partitioning the trend PN1 ----
PedEval1$EbvI = 0.5 * (PedEval1$EbvT1 + PedEval1$EbvT2)
PedEval1$EbvT1_s <- (PedEval1$EbvT1 - mean(PedEval1$EbvT1[PedEval1$Generation == 5])) / sd(PedEval1$EbvT1[PedEval1$Generation == 5])
PedEval1$EbvT2_s <- (PedEval1$EbvT2 - mean(PedEval1$EbvT2[PedEval1$Generation == 5])) / sd(PedEval1$EbvT2[PedEval1$Generation == 5])
#PedEval1$EbvI_s <- (PedEval1$EbvI - mean(PedEval1$EbvI[PedEval1$Generation == 0])) / sd(PedEval1$EbvI[PedEval1$Generation == 0])
PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)

# ---- Partitioning the trend PN2 ----
PedEval2$EbvI = 0.5 * (PedEval2$EbvT1 + PedEval2$EbvT2)
PedEval2$EbvT1_s <- (PedEval2$EbvT1 - mean(PedEval2$EbvT1[PedEval2$Generation == 0])) / sd(PedEval2$EbvT1[PedEval2$Generation == 0])
PedEval2$EbvT2_s <- (PedEval2$EbvT2 - mean(PedEval2$EbvT2[PedEval2$Generation == 0])) / sd(PedEval2$EbvT2[PedEval2$Generation == 0])
#PedEval2$EbvI_s <- (PedEval2$EbvI - mean(PedEval2$EbvI[PedEval2$Generation == 0])) / sd(PedEval2$EbvI[PedEval2$Generation == 0])
PedEval2$EbvI_s <- 0.5 * (PedEval2$EbvT1_s + PedEval2$EbvT2_s)


PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")
PedEval1$GenerationProgram <- paste(PedEval1$Generation, PedEval1$Program, sep="-")
PedEval2$GenerationProgram <- paste(PedEval2$Generation, PedEval2$Program, sep="-")
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = TRUE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1", "EbvT2", "EbvI"))

Part2 = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))

Part1Summary = summary(object = Part1, by = "Generation")
Part2Summary = summary(object = Part2, by = "Generation")


p1 <- plot(Part1Summary)
p2 <- plot(Part2Summary)
#poreži živali iz nukleusa
Part1PN <- Part1
Part1PN$EbvT1_s <- Part1PN$EbvT1_s[Part1PN$EbvT1_s$Program == "PN1",]
Part1PN$EbvT2_s <- Part1PN$EbvT2_s[Part1PN$EbvT2_s$Program == "PN1",]
Part1PN$EbvI_s <- Part1PN$EbvI_s[Part1PN$EbvI_s$Program == "PN1",]

Part1GN <- Part1
Part1GN$EbvT1_s <- Part1GN$EbvT1_s[Part1GN$EbvT1_s$Program == "GN",]
Part1GN$EbvT2_s <- Part1GN$EbvT2_s[Part1GN$EbvT2_s$Program == "GN",]
Part1GN$EbvI_s <- Part1GN$EbvI_s[Part1GN$EbvI_s$Program == "GN",]

Part1SummaryPN = summary(object = Part1PN, by = "Generation")
p1cPN <- plot(Part1SummaryPN)
Part1SummaryGN = summary(object = Part1GN, by = "Generation")
p1cGN <- plot(Part1SummaryGN)



Part2PN <- Part2
Part2PN$EbvT2_s <- Part2PN$EbvT2_s[Part2PN$EbvT2_s$Program == "PN2",]
Part2PN$EbvT2_s <- Part2PN$EbvT2_s[Part2PN$EbvT2_s$Program == "PN2",]
Part2PN$EbvI_s <- Part2PN$EbvI_s[Part2PN$EbvI_s$Program == "PN2",]

Part2GN <- Part2
Part2GN$EbvT2_s <- Part2GN$EbvT2_s[Part2GN$EbvT2_s$Program == "GN",]
Part2GN$EbvT2_s <- Part2GN$EbvT2_s[Part2GN$EbvT2_s$Program == "GN",]
Part2GN$EbvI_s <- Part2GN$EbvI_s[Part2GN$EbvI_s$Program == "GN",]

Part2SummaryPN = summary(object = Part2PN, by = "Generation")
p2cPN <- plot(Part2SummaryPN)
Part2SummaryGN = summary(object = Part2GN, by = "Generation")
p2cGN <- plot(Part2SummaryGN)

#savePlot(p1, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_Gender", type = "jpeg", res=601, width=90, height=70, units="mm")
#savePlot(x = p1, filename = "Partition_Gender", type = "jpeg")



# ---- Partitioning the trend PN1 ----
PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
PedEval1$TbvT1_s <- (PedEval1$TbvT1 - mean(PedEval1$TbvT1[PedEval1$Generation == 0])) / sd(PedEval1$TbvT1[PedEval1$Generation == 0])
PedEval1$TbvT2_s <- (PedEval1$TbvT2 - mean(PedEval1$TbvT2[PedEval1$Generation == 0])) / sd(PedEval1$TbvT2[PedEval1$Generation == 0])
#PedEval1$TbvI_s <- (PedEval1$TbvI - mean(PedEval1$TbvI[PedEval1$Generation == 0])) / sd(PedEval1$TbvI[PedEval1$Generation == 0])
PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)

# ---- Partitioning the trend PN2 ----
PedEval2$TbvI = 0.5 * (PedEval2$TbvT1 + PedEval2$TbvT2)
PedEval2$TbvT1_s <- (PedEval2$TbvT1 - mean(PedEval2$TbvT1[PedEval2$Generation == 0])) / sd(PedEval2$TbvT1[PedEval2$Generation == 0])
PedEval2$TbvT2_s <- (PedEval2$TbvT2 - mean(PedEval2$TbvT2[PedEval2$Generation == 0])) / sd(PedEval2$TbvT2[PedEval2$Generation == 0])
#PedEval2$TbvI_s <- (PedEval2$TbvI - mean(PedEval2$TbvI[PedEval2$Generation == 0])) / sd(PedEval2$TbvI[PedEval2$Generation == 0])
PedEval2$TbvI_s <- 0.5 * (PedEval2$TbvT1_s + PedEval2$TbvT2_s)


PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = TRUE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("TbvT1", "TbvT2", "TbvI"))
Part2 = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))



#poreži živali iz nukleusa
Part1PNt <- Part1
Part1PNt$EbvT1_s <- Part1PNt$TbvT1_s[Part1PNt$TbvT1_s$Program == "PN1",]
Part1PNt$EbvT2_s <- Part1PNt$TbvT2_s[Part1PNt$TbvT2_s$Program == "PN1",]
Part1PNt$EbvI_s <- Part1PNt$TbvI_s[Part1PNt$TbvI_s$Program == "PN1",]

Part1GNt <- Part1
Part1GNt$TbvT1_s <- Part1GNt$TbvT1_s[Part1GNt$TbvT1_s$Program == "PN1",]
Part1GNt$TbvT2_s <- Part1GNt$TbvT2_s[Part1GNt$TbvT2_s$Program == "PN1",]
Part1GNt$TbvI_s <- Part1GNt$TbvI_s[Part1GNt$TbvI_s$Program == "PN1",]

Part1SummaryPN = summary(object = Part1PNt, by = "Generation")
p1cPNt <- plot(Part1SummaryPN)
Part1SummaryGN = summary(object = Part1GNt, by = "Generation")
p1cGNt <- plot(Part1SummaryGN)


Part2PNt <- Part2
Part2PNt$TbvT2_s <- Part2PNt$TbvT2_s[Part2PNt$TbvT2_s$Program == "PN2",]
Part2PNt$TbvT2_s <- Part2PNt$TbvT2_s[Part2PNt$TbvT2_s$Program == "PN2",]
Part2PNt$TbvI_s <- Part2PNt$TbvI_s[Part2PNt$TbvI_s$Program == "PN2",]

Part2GNt <- Part2
Part2GNt$TbvT2_s <- Part2GNt$TbvT2_s[Part2GNt$TbvT2_s$Program == "PN2",]
Part2GNt$TbvT2_s <- Part2GNt$TbvT2_s[Part2GNt$TbvT2_s$Program == "PN2",]
Part2GNt$TbvI_s <- Part2GNt$TbvI_s[Part2GNt$TbvI_s$Program == "PN2",]

Part2SummaryPN = summary(object = Part2PNt, by = "Generation")
p2cPNt <- plot(Part2SummaryPN)
Part2SummaryGN = summary(object = Part2GNt, by = "Generation")
p2cGNt <- plot(Part2SummaryGN)

#poreži živali iz nukleusa
Part1$TbvT1_s <- Part1$TbvT1_s[Part1$TbvT1_s$Program == "PN1",]
Part1$TbvT2_s <- Part1$TbvT2_s[Part1$TbvT2_s$Program == "PN1",]
Part1$TbvI_s <- Part1$TbvI_s[Part1$TbvI_s$Program == "PN1",]

Part2$TbvT1_s <- Part2$TbvT1_s[Part2$TbvT1_s$Program == "PN2",]
Part2$TbvT2_s <- Part2$TbvT2_s[Part2$TbvT2_s$Program == "PN2",]
Part2$TbvI_s <- Part2$TbvI_s[Part2$TbvI_s$Program == "PN2",]


Part1Summary = summary(object = Part1, by = "Generation")
Part2Summary = summary(object = Part2, by = "Generation")
p1g <- plot(Part1Summary)
p2g <- plot(Part2Summary)
savePlot(p1, filename = "/home/jana/Documents/PhD/Projects/inProgress/AlphaPartition//Figures/Partition_Gender", type = "jpeg", res=601, width=90, height=70, units="mm")
savePlot(x = p1, filename = "Partition_Gender", type = "jpeg")

PedEval1 <- read.table("PedEval_NewModel_0.25.csv", header=TRUE)
detach("package:plyr")
library(dplyr)
PN11 <- PedEval1 %>% 
  group_by(Generation, Program, Gender) %>%
  summarize(COR=cor(TbvT1, EbvT1))

ggplot(as.data.frame(PN11), aes(x=Generation, y=COR, group=Gender, colour=Gender)) + geom_line() + facet_grid(. ~ Program)

PN12 <- PedEval1 %>% 
  group_by(Generation, Program, Gender) %>%
  summarize(COR=cor(TbvT2, EbvT2))
ggplot(as.data.frame(PN12), aes(x=Generation, y=COR, group=Gender, colour=Gender)) + geom_line() + facet_grid(. ~ Program)

PN21 <- PedEval2 %>% 
  group_by(Generation, Program) %>%
  summarize(COR=cor(TbvT1, EbvT1))

PN22 <- PedEval2 %>% 
  group_by(Generation, Program) %>%
  summarize(COR=cor(TbvT2, EbvT2))


#Korelacija PA!
head(PedEval1)
PedEval1 <- PedEval1[,1:13]
# Ebv & Tbv PN1
PN1Fathers <- unique(PedEval1[PedEval1$IId %in% PedEval1$FId, c("IId", "EbvT1", "EbvT2", "TbvT1", "TbvT2")])
colnames(PN1Fathers) <- c("FId", "F_EbvT1", "F_EbvT2", "F_TbvT1", "F_TbvT2" )
PN1Mothers <- unique(PedEval1[PedEval1$IId %in% PedEval1$MId, c("IId", "EbvT1", "EbvT2", "TbvT1", "TbvT2")])
colnames(PN1Mothers) <- c("MId", "M_EbvT1", "M_EbvT2", "M_TbvT1", "M_TbvT2")
length(PN1Fathers$FId %in% PedEval1$FId)
length(PN1Mothers$MId %in% PedEval1$MId)
PedEval1 <- merge(PedEval1, PN1Fathers, by="FId", all.x=TRUE)
PedEval1 <- merge(PedEval1, PN1Mothers, by="MId", all.x=TRUE)
PedEval1$PAEbv1 <- 0.5 * (PedEval1$F_EbvT1 + PedEval1$M_EbvT1)
PedEval1$PAEbv2 <- 0.5 * (PedEval1$F_EbvT2 + PedEval1$M_EbvT2)
PedEval1$PATbv1 <- 0.5 * (PedEval1$F_TbvT1 + PedEval1$M_TbvT1)
PedEval1$PATbv2 <- 0.5 * (PedEval1$F_TbvT2 + PedEval1$M_TbvT2)
PedEval1$MSTEbv1 <- PedEval1$EbvT1  - PedEval1$PAEbv1
PedEval1$MSTEbv2 <- PedEval1$EbvT2  - PedEval1$PAEbv2
PedEval1$MSTTbv1 <- PedEval1$TbvT1  - PedEval1$PATbv1
PedEval1$MSTTbv2 <- PedEval1$TbvT1  - PedEval1$PATbv1

head(PedEval1)
#calculate BIAS of MSTs
PedEval1$biasMST1 <- PedEval1$MSTEbv1 - PedEval1$MSTTbv1
PedEval1$biasMST2 <- PedEval1$MSTEbv2 - PedEval1$MSTTbv2
PedEval1$biasEBV1 <- PedEval1$EbvT1 - PedEval1$TbvT1
PedEval1$biasEBV2 <- PedEval1$EbvT2 - PedEval1$TbvT2

write.csv(PedEval1, "PN1/PedEval_Intercept_0.25.csv", quote=FALSE, row.names=FALSE)
write.csv(PedEval1, "PN1/PedEval_NewModel_0.25.csv", quote=FALSE, row.names=FALSE)
model <- lm(PedEval1$TbvT1 ~ PedEval1$EbvT1 + PedEval1$Generation + PedEval1$ProgramGender)
summary(model)


plot(model)

library(dplyr)
detach("package:plyr")
  #negative means underestimation
PN11 <- as.data.frame(PedEval1 %>% 
  group_by(Generation, ProgramGender) %>% 
  summarise(MEAN=mean(biasEBV1)))

PN11 <- as.data.frame(PedEval1 %>% 
  group_by(Generation, ProgramGender) %>% 
  summarise(MEAN=mean(biasMST1)))

PN11 <- as.data.frame(PedEval1 %>% 
  group_by(Generation, ProgramGender) %>% 
  summarise(MEAN=mean(MSTEbv1)))

PN12 <- as.data.frame(PedEval1 %>% 
  group_by(Generation, ProgramGender) %>% 
  summarise(MEAN=mean(biasMST2)))
PN11 <- summarySE(data = PedEval1, groupvars = c("Generation", "ProgramGender"), measurevar = "biasMST1")

ggplot(data = PN11, aes(x=Generation, y=MEAN, colour=ProgramGender, group=ProgramGender)) + geom_line()
ggplot(data = PedEval1, aes(x=Generation, y=biasMST1, colour=ProgramGender, group=ProgramGender)) + geom_point() + geom_smooth()
ggplot(data = PN11, aes(x=Generation, y=biasMST1, colour=ProgramGender, group=ProgramGender)) + geom_line()
maleGen1 <- PedEval1[PedEval1$Generation == 1 & PedEval1$Gender == "M",]
head(maleGen1[order(maleGen1$biasMST1),])


#MST of the fathers

FATH <- PedEval1[PedEval1$IId %in% unique(PedEval1$FId),]
ggplot(data = FATH, aes(x=Generation, y=biasMST1, colour=Program, group=Program)) + geom_point() + geom_smooth()
ggplot(data = FATH, aes(x=Generation, y=biasMST2, colour=Program, group=Program)) + geom_point() + geom_smooth()
MOTH <- PedEval1[PedEval1$IId %in% unique(PedEval1$MId),]
ggplot(data = MOTH, aes(x=Generation, y=biasMST1)) + geom_point() + geom_smooth() + facet_grid(.~Program)
ggplot(data = MOTH, aes(x=Generation, y=biasMST2)) + geom_point() + geom_smooth() + facet_grid(.~Program)



#which animals have the lowest GN-M partition


cor(PedEval1$PAEbv1, PedEval1$PATbv1, use="pairwise.complete.obs")
cor(PedEval1$PAEbv2, PedEval1$PATbv2, use="pairwise.complete.obs")
cor(PedEval1$MSTEbv1, PedEval1$MSTTbv1, use="pairwise.complete.obs")
cor(PedEval1$MSTEbv2, PedEval1$MSTTbv2, use="pairwise.complete.obs")


PN11 <- PedEval1 %>% group_by(Generation, Gender, Program) %>% summarise(COR=cor(TbvT1, EbvT1))
PN12 <- PedEval1 %>% group_by(Generation, Gender, Program) %>% summarise(COR=cor(TbvT2, EbvT2))
options(dplyr.print_max = 1e9)

qplot(data=PN11, x=Generation, y=COR, group=c("Gender"), colour=Gender, shape=Program,size=Program, geom="point") 
qplot(data=PN12, x=Generation, y=COR, group=c("Gender"), colour=Gender, shape=Program,size=Program, geom="point") 

###################################################################################################
###################################################################################################
####gibbsf90
#PN1 - trait 1
g <- read.csv("../PN1/Gibbs1.csv")[,-1]
colnames(g) <- c("IId", paste0("EBV1_", 1:10))
g <- g[g$IId %in% PedEval1$IId,]
PedEval1g <- merge(PedEval1, g, by="IId")

for (sample in 1:samples) {
  col = paste0("EBV1_", sample)
  col1 = paste0("EBV1_", sample, "_s")
  PedEval1g[[col1]] <- (PedEval1g[[col]] - mean(PedEval1g[[col]][PedEval1g$Generation == 0])) / sd(PedEval1g[[col]][PedEval1g$Generation == 0])
}


Part1gibbs = AlphaPart(x = as.data.frame(PedEval1g), sort = FALSE,
                       colId = "IId", colFid = "FId", colMid = "MId",
                       colPath = "ProgramGender", colAGV = c(paste0("EBV1_", 1:10, "_s")))
Part1Summary_gibbs = summary(object = Part1gibbs, by = "Generation")
print(Part1Summary_gibbs)

samples = 10
gibbsPart <- data.frame()
for (sample in 1:samples) {
  col = paste0("EBV1_", sample, "_s")
  tmp <- Part1Summary_gibbs[[col]]$abs
  gibbsPart <- rbind(gibbsPart, tmp)
}

avgGibbsPart <- data.frame(Generation=0:20)
for (path in c("GN-F", "GN-M", "PN1-F", "PN1-M")) {
  tmp <- summarySE(data = gibbsPart, measurevar = path, groupvars = "Generation")[,c(1,3,4)]
  colnames(tmp)[3] <- paste0("SD_", path)
  avgGibbsPart <- cbind(avgGibbsPart, tmp)
}


avgGibbsPart_m <- melt(avgGibbsPart, id.vars = "Generation")
avgGibbsPart_Mean <- avgGibbsPart_m[avgGibbsPart_m$variable %in% c("GN-F", "GN-M", "PN1-F", "PN1-M"),]
ggplot(data = avgGibbsPart_Mean, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line() 

############################
############################
setwd("/home/jana/Documents/SimulationAlphaPart/")
PedEval1 <- read.table("PedEval1.csv", header=TRUE)
#check the bias
PedEval1$bias1 <- PedEval1$TbvT1 - PedEval1$EbvT1
PedEval1$bias2 <- PedEval1$TbvT2 - PedEval1$EbvT2
summarySE(data = PedEval1, groupvars = "Gender", measurevar = "bias1")
summarySE(data = PedEval1, groupvars = "Gender", measurevar = "bias2")

library(dplyr)
PedEval1 %>% 
  group_by(Generation) %>%
  summarize(COR=cor(TbvT1, EbvT1))


PedEval1$ProgramGender <- paste(PedEval1$Program, PedEval1$Gender)
#PedEval1 <- PedEval1[,1:22]
#PedEval1Copy <- PedEval1


#PN1 - trait 2
library(tidyr)
library(AlphaPart)
for (split in 0:89) {
  print(paste0("SPLIT: ", split))
  #read in a split portion of the solutions
  gR <- read.table("GIBBS2_sort.csv", skip = split * 3360500, nrows = 3360500)[,-1]
  colnames(gR) <- c("ID", "EBV2", "Sample")
  g <- spread(gR, Sample, EBV2)
  
  #rename the spread data frame 
  numSample <- (ncol(g)-1)
  Samples <- unique(gR$Sample)
  colnames(g) <- c("IId", paste0("EBV2_", Samples))
  g <- g[g$IId %in% PedEval1$IId,]
  
  #merge with the original PedEval for the program
  PedEval1 <- merge(PedEval1, g, by="IId")
  
  #standardise all "traits"
  for (sample in Samples) {
    col = paste0("EBV2_", sample)
    col1 = paste0("EBV2_", sample, "_s")
    PedEval1[[col1]] <- (PedEval1[[col]] - mean(PedEval1[[col]][PedEval1$Generation == 0])) / sd(PedEval1[[col]][PedEval1$Generation == 0])
  }
  
  #partition the current number of samples ("traits")
  Part1gibbs2 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                          colId = "IId", colFid = "FId", colMid = "MId",
                          colPath = "ProgramGender", colAGV = c(paste0("EBV2_", Samples, "_s")))
  
  #summarize the partitions
  Part1Summary_gibbs2 = summary(object = Part1gibbs2, by = "Generation")
  #print(Part1Summary_gibbs2)
  
  if (split == 0) {
    gibbsPart2 <- data.frame()
    # }
    
    for (sample in Samples) {
      col = paste0("EBV2_", sample, "_s")
      tmp <- Part1Summary_gibbs2[[col]]$abs
      tmp$Sample <- sample
      gibbsPart2 <- rbind(gibbsPart2, tmp)
    }
    
  }
  
  
  #average the partition sums
  avgGibbsPart2 <- data.frame(Generation=0:20)
  for (path in c("GN-F", "GN-M", "PN1-F", "PN1-M")) {
    tmp <- summarySE(data = gibbsPart2, measurevar = path, groupvars = "Generation")[,c(1,3,4)]
    colnames(tmp)[3] <- paste0("SD2_", path)
    avgGibbsPart2 <- cbind(avgGibbsPart2, tmp)
  } 
  
  avgGibbsPart2_m <- melt(avgGibbsPart2, id.vars = "Generation")
  avgGibbsPart2_Mean <- avgGibbsPart2_m[avgGibbsPart2_m$variable %in% c("GN-F", "GN-M", "PN1-F", "PN1-M"),]
  ggplot(data = avgGibbsPart2_Mean, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line() + ggtitle("Trait 2")
  
  ###########
  #small example
  library(AlphaPart)
  library(AlphaSimR)
  ped = data.frame(ID = c("A", "B", "C", "T", "E", "D", "U", "V"),
                   Mother = c(NA, NA, "A", NA, "C", "C", NA, "E"),
                   Father = c(NA, NA, "B", "B", "D", NA, "D", NA),
                   Country = c("X", "Y", "X", "Y", "X", "Y", "Y", "X"))
  
  #library(AlphaSimR)
  founderPop = runMacs(nInd = 8,
                       nChr = 10,
                       segSites = 1000,
                       nThreads = 4,
                       # species = "GENERIC")
                       species = "CATTLE")
  
  SP = SimParam$new(founderPop)
  # VarA = matrix(data = c(1.0, 0.1, 0.1, 1.0), nrow = 2); cov2cor(VarA)
  VarA = matrix(data = c(1.0, 0.0, 0.0, 1.0), nrow = 2); cov2cor(VarA)
  VarE = matrix(data = c(3.0, 0.0, 0.0, 9.0), nrow = 2); cov2cor(VarE)
  
  VarP = VarA + VarE; diag(VarA) / diag(VarP)
  SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0), var = diag(VarA), cor = cov2cor(VarA))
  # SP$addSnpChip(nSnpPerChr = 1000)
  SP$setGender(gender = "yes_rand")
  
  
  TestPop = newPop(founderPop)
  TestPop = setPheno(pop = TestPop, varE = VarE)
  blupPed = ped
  setwd("/home/jana/Documents/SimulationAlphaPart/Test/")
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, na = "0")
  TestPop@pheno[,1]
  pedTmp <- ped
  pedTmp$pheno1 <- TestPop@pheno[,1]
  pedTmp$pheno2 <- TestPop@pheno[,2]
  write.table(pedTmp[,c("ID", "pheno1")], "Blupf901.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, na = "0")
  write.table(pedTmp[,c("ID", "pheno2")], "Blupf902.dat", quote=FALSE, row.names=FALSE, col.names=FALSE, na = "0")
  
  system(command = "./renumf90 < renumParam1")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1" ))
  
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2" ))
  
  ##check correlation
  sol1 <- read.table("renumbered_Solutions_1")
  colnames(sol1) <- c("renID", "ID", "EBV1")
  sol2 <- read.table("renumbered_Solutions_2")
  colnames(sol2) <- c("renID", "ID", "EBV2")
  ped <- merge(ped, sol1[,2:3], by= "ID")
  ped <- merge(ped, sol2[,2:3], by= "ID")
  ped$TGV1 <- TestPop@gv[,1]
  ped$TGV2 <- TestPop@gv[,2]
  ped$Gender = TestPop@gender
  
  library(dplyr)
  ped$EBV1 <- as.numeric(ped$EBV1)
  ped$EBV2 <- as.numeric(ped$EBV2)
  ped$TGV1 <- as.numeric(ped$TGV1)
  ped$TGV2 <- as.numeric(ped$TGV2)
  ped$Country <- as.factor(ped$Country)
  
  
  library(plyr)
  cor(ped$EBV1[ped$Country == "X"], ped$TGV1[ped$Country == "X"])
  cor(ped$EBV1[ped$Country == "Y"], ped$TGV1[ped$Country == "Y"])
  cor(ped$EBV1[ped$Gender == "M"], ped$TGV1[ped$Gender == "M"])
  cor(ped$EBV1[ped$Gender == "F"], ped$TGV1[ped$Gender == "F"])
  
  ped$bias1 <- ped$TGV1 - ped$EBV1
  ped$bias2 <- ped$TGV2 - ped$EBV2
  
  
  ped %>%
    group_by(Gender) %>%
    summarize(COR=cor(TGV1, EBV1))
  ped %>%
    group_by(Gender) %>%
    summarize(COR=cor(TGV2, EBV2))
  
  testPartG <- AlphaPart( ped, sort = FALSE,
                          colId = "ID", colFid = "Father", colMid = "Mother",
                          colPath = "Country", colAGV = c("TGV1", "TGV2"))
  sumTestPartG <- summary(object = testPartG, by = "Gender")
  plot(sumTestPartG)
  
  testPart <- AlphaPart( ped, sort = FALSE,
                         colId = "ID", colFid = "Father", colMid = "Mother",
                         colPath = "Country", colAGV = c("EBV1", "EBV2"))
  
  sumTestPart <- summary(object = testPart, by = "Gender")
  plot(sumTestPart)
  
  
  ##################
  dat2 <- read.table("Blupf902.dat")
pedF <- PedEval1[PedEval1$IId %in% dat2$V1,]  
table(pedF$Program)
table(pedF$Gender)

as.data.frame(PedEval1 %>% 
  group_by(Generation, Gender, Program) %>%
  summarize(COR=sd(EbvT1)))
