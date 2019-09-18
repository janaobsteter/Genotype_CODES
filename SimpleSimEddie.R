# 
# ---- Environment ----

rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
REP = args[1]

library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")
library(pedigree)


homedir = getwd()
# ---- General parameters ----

nGNMales   =  100
nGNFemales = 100


nGenerationBurn = 20
nGenerationEval = 20

GenMeanCols = c("GenMeanT1")
GenVarCols  = c("GenVarT1")

# ---- Base population genomes ----

founderPop = runMacs(nInd = nGNMales + nGNFemales,
                     nChr = 10,
                     segSites = 1000,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")

MST <- data.frame(Rep = NA, H2 = NA, Generation = NA, Gender = NA, rMST = NA)
Acc <- data.frame(Rep = NA, H2 = NA, Generation = NA, Gender = NA, rEBV = NA)
###################################################################################
###################################################################################
# ---- Simulation/Base population parameters ----


system(paste0("mkdir ", homedir, "/Rep", REP))
RepDir = paste0(homedir, "/Rep", REP, "/")
setwd(RepDir)

for (h2 in c(0.25)) {
  system(paste0('mkdir ', RepDir, '/H_', h2))
  hDir = paste0(RepDir, '/H_', h2, "/")
  setwd(hDir)
  system(paste0('cp ', homedir, '/Essentials/* ', hDir))
  
  SP = SimParam$new(founderPop)
  # VarA = matrix(data = c(1.0, 0.1, 0.1, 1.0), nrow = 2); cov2cor(VarA)
  VarA = matrix(data = c(1.0), nrow = 1); cov2cor(VarA)
  VarE = matrix(data = c(3.0), nrow = 1); cov2cor(VarE)
  VarP = VarA + VarE; diag(VarA) / diag(VarP)
  SP$addTraitA(nQtlPerChr = 1000, mean = c(0), var = diag(VarA), cor = cov2cor(VarA))
  # SP$addSnpChip(nSnpPerChr = 1000)
  SP$setGender(gender = "yes_sys")
  
  # ---- Base GN population ----
  
  GN = newPop(founderPop)
  BaseGNMales   = GN[GN@gender == "M"]
  BaseGNFemales = GN[GN@gender == "F"]
  rm(GN)
  
  ###################################################################################
  ###################################################################################
  
  # ---- GN burn-in ----
  
  DataBurn = tibble(Generation = rep(1:nGenerationBurn, each=2), Gender = rep(c("F", "M"), nGenerationBurn),
                    GenMeanT1 = NA, 
                    GenVarT1  = NA)
  for (Generation in 1:nGenerationBurn) {
    # Generation = 1
    
    # Mate
    SelCand = randCross2(females = BaseGNFemales, males = BaseGNMales,
                         nCrosses = BaseGNFemales@nInd, nProgeny = 12)
    
    # Save metrics
    for (gender in c("F", "M")) {
      DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender), GenMeanCols] =
        mean(SelCand@gv[SelCand@gender == gender,])
      DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender), GenVarCols] =
        var(SelCand@gv[SelCand@gender == gender,])
    }
    
    
    # Phenotype
    SelCand = setPheno(pop = SelCand, varE = VarE)
    
    # Evaluate (could run BLUP)
    SelCand@ebv = SelCand@pheno
    
    # Select
    BaseGNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                              use = "ebv", trait = function(x) rowMeans(x))
    BaseGNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
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
  
  GNInd = c(BaseGNFemales@id, BaseGNMales@id)
  SelCand = SelCand[!(SelCand@id %in% GNInd)]; SelCand@nInd
  # Cheat here - consider all animals are females
  SelCand@gender[] = "F"
  
  
  
  # ---- Evaluation ----
  
  PedEval = rbind(tibble(Generation = 0,
                         IId        = BaseGNMales@id,
                         FId        = NA,
                         MId        = NA,
                         Gender     = BaseGNMales@gender,
                         Program    = "GN",
                         PhenoT1    = BaseGNMales@pheno[,1],
                         EbvT1      = BaseGNMales@ebv[, 1],
                         TbvT1      = BaseGNMales@gv[, 1]),
                  tibble(Generation = 0,
                         IId        = BaseGNFemales@id,
                         FId        = NA,
                         MId        = NA,
                         Gender     = BaseGNFemales@gender,
                         Program    = "GN",
                         PhenoT1    = BaseGNFemales@pheno[,1],
                         EbvT1      = BaseGNFemales@ebv[, 1],
                         TbvT1      = BaseGNFemales@gv[, 1])
  )
  
  
  
  DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                      Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                      GenMeanT1 = NA, 
                      GenVarT1  = NA)
  DataEvalGN$Program = "GN"
  
  ##############################################################################3
  ##############################################################################3
  # ---- Program PN1  ----
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
        mean(SelCand@gv[SelCand@gender == gender,])
      DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenVarCols] =
        var(SelCand@gv[SelCand@gender == gender,])
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
                           EbvT1      = NA,
                           TbvT1      = SelCand@gv[, 1]))
    
    # Estimate EBVs with blupf90
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
    blupdat <- PedEval[,c("IId", "PhenoT1", "Generation", "Gender")]
    #    blup1dat$Generation[blup1dat$Generation == 0] <- "00"
    write.table(blupdat, "Blupf90.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf90.dat"))
    
    
    
    
    
    # Evaluate Trait 1
    system(command = "./renumf90 < renumParam")
    system(command = "./blupf90 renf90.par")
    system(command = "bash Match_AFTERRenum.sh")
    system(command = paste0("mv renumbered_Solutions renumbered_Solutions_", Generation))
    
    sol <- read.table(paste0("renumbered_Solutions_", Generation))[,2:3]
    sol <- sol[order(match(sol$V2, PedEval$IId)),]
    PedEval$EbvT1 <- sol$V3
    
  
    
    
    # Set EBVs for SelCand
    SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1")])
    
    # Select
    GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                          use = "ebv", trait = function(x) rowMeans(x))
    GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                          use = "ebv", trait = function(x) rowMeans(x))
    # Clean
    rm(SelCand)
    
    
    
  }
  write.table(PedEval, paste0("PedEval_", h2, ".csv"), quote=FALSE, row.names=FALSE)
  write.table(DataEvalGN, paste0("DataEval_", h2, ".csv"), quote=FALSE, row.names=FALSE)
  
 

  
  detach(package:plyr)
  library(dplyr)
  PN11 <- as.data.frame(PedEval %>% 
    dplyr::group_by(Generation, Gender) %>%
    summarize(COR=cor(TbvT1, EbvT1)))
  
  Acc <- rbind(Acc, data.frame(Rep = REP, H2 = h2, Generation = PN11$Generation, Gender = PN11$Gender, rEBV = PN11$COR))

  
  
  #Korelacija PA!
  # Ebv
  PN1Fathers <- unique(PedEval[PedEval$IId %in% PedEval$FId, c("IId", "EbvT1", "TbvT1")])
  colnames(PN1Fathers) <- c("FId", "F_EbvT1",  "F_TbvT1")
  PN1Mothers <- unique(PedEval[PedEval$IId %in% PedEval$MId, c("IId", "EbvT1",  "TbvT1")])
  colnames(PN1Mothers) <- c("MId", "M_EbvT1",  "M_TbvT1")
  PedEval <- merge(PedEval, PN1Fathers, by="FId", all.x=TRUE)
  PedEval$Generation <- as.numeric(PedEval$Generation)
  PedEval <- PedEval[order(PedEval$Generation),]
  PedEval <- merge(PedEval, PN1Mothers, by="MId", all.x=TRUE)
  PedEval$PAEbv1 <- 0.5 * (PedEval$F_EbvT1 + PedEval$M_EbvT1)
  PedEval$PATbv1 <- 0.5 * (PedEval$F_TbvT1 + PedEval$M_EbvT1)
  PedEval <- PedEval[order(PedEval$IId),]

  
  PedEval$MSTe <- PedEval$PAEbv1 - PedEval$EbvT1
  PedEval$MSTt <- PedEval$PATbv1 - PedEval$TbvT1

  MSTCor <- as.data.frame(PedEval %>% 
    dplyr::group_by(Generation, Gender) %>%
    summarize(COR=cor(MSTe, MSTt)))
  
  
  MST <- rbind(MST, data.frame(Rep = REP, H2 = h2, Generation = MSTCor$Generation, Gender = MSTCor$Gender, rMST = MSTCor$COR))
}

write.table(MST, paste0(homedir, "MST_", h2, ".csv"), quote=FALSE, row.names=FALSE)
write.table(Acc, paste0(homedir, "Acc", h2, ".csv"), quote=FALSE, row.names=FALSE)

