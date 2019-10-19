# 
# ---- Environment ----
args = commandArgs(trailingOnly=TRUE)
REP <- args[1]
print(REP)

setwd("/home/jana/Documents/PhD/SimulationAlphaPart/NewModelAllPed/")
library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")
library(pedigree)


#args = commandArgs(trailingOnly=TRUE)
#varE1 = args[1]
#varE2 = args[2]

homeDir = getwd()
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
                     segSites = 10,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")

###################################################################################
###################################################################################
# ---- Simulation/Base population parameters ----
for (rep in c(REP)) {
  system(paste0("mkdir ", homeDir, "/Rep", rep))
  RepDir = paste0(homeDir, "/Rep", rep, "/")
  setwd(RepDir)
for (h2 in c(0.25)) {
  system(paste0('mkdir ', RepDir, '/H_', h2))
  hDir = paste0(RepDir, '/H_', h2, "/")
  system(paste0('cp Essentials/* ', hDir))
  setwd(hDir)
  system('mkdir GN PN1 PN2')
  system(paste0('cp ', homeDir, '/Essentials/* ', hDir, '/PN1/'))
  system(paste0('cp ', homeDir, '/Essentials/* ', hDir, '/PN2/'))
  system(paste0('cp ', homeDir, '/Essentials/* ', hDir, '/GN/'))
  print(getwd())
  SP = SimParam$new(founderPop)
  # VarA = matrix(data = c(1.0, 0.1, 0.1, 1.0), nrow = 2); cov2cor(VarA)
  VarA = matrix(data = c(1.0, 0.0, 0.0, 1.0), nrow = 2); cov2cor(VarA)
  #varE = 1/ h2 -1
  VarE = matrix(data = c(3.0, 0.0, 0.0, 9.0), nrow = 2); cov2cor(VarE)
  
  VarP = VarA + VarE; diag(VarA) / diag(VarP)
  SP$addTraitA(nQtlPerChr = 10, mean = c(0, 0), var = diag(VarA), cor = cov2cor(VarA))
  # SP$addSnpChip(nSnpPerChr = 1000)
  SP$setGender(gender = "yes_rand")
  
  # ---- Base GN population ----
  
  GN = newPop(founderPop)
  BaseGNMales   = GN[GN@gender == "M"]
  BaseGNFemales = GN[GN@gender == "F"]
  rm(GN)
  
  ###################################################################################
  ###################################################################################
  
  # ---- GN burn-in ----
  
  DataBurn = tibble(Generation = rep(1:nGenerationBurn, each=2), Gender = rep(c("F", "M"), nGenerationBurn),
                    GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                    GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
  PedEval <- tibble()
  
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
    
    PedEval = rbind(PedEval,
                    tibble(Generation = Generation,
                           IId        = SelCand@id,
                           FId        = SelCand@father,
                           MId        = SelCand@mother,
                           Gender     = SelCand@gender,
                           Program    = "BurnIn",
                           PhenoT1    = SelCand@pheno[,1],
                           PhenoT2    = SelCand@pheno[,2],
                           EbvT1      = NA,
                           EbvT2      = NA,
                           TbvT1      = SelCand@gv[, 1],
                           TbvT2      = SelCand@gv[, 2]))
    
    # Select
    BaseGNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                              use = "pheno", trait = function(x) rowMeans(x))
    BaseGNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                              use = "pheno", trait = function(x) rowMeans(x))
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
  BasePNFemales = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                            use = "pheno", trait = function(x) rowMeans(x))
  # meanG(GNMales); meanG(GNFemales); meanG(PNFemales)
  
  
  # ---- PN1 ----
  
  PedEvalBurnIn <- PedEval
  
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
  hDir = getwd()
  PedEval <- PedEvalBurnIn  
  # ---- Program PN1  ----
  for (Generation in (1 + nGenerationBurn):(nGenerationEval + nGenerationBurn)) {
    GNEffect1 <- rnorm(1, 0, 0.2 * diag(VarP)[1])
    GNEffect2 <- rnorm(1, 0, 0.2 * diag(VarP)[2])
    PNEffect1 <- rnorm(1, 0, 0.2 * diag(VarP)[1])
    PNEffect2 <- rnorm(1, 0, 0.2 * diag(VarP)[2])
    GenerationEffect1 <- rnorm(1, 0, 0.2 * diag(VarP)[1])
    GenerationEffect2 <- rnorm(1, 0, 0.2 * diag(VarP)[2])

    if (Generation == (1 + nGenerationBurn)) {
      GNFemales = BaseGNFemales
      GNMales = BaseGNMales
      PedEval$Program[PedEval$IId %in% BaseGNFemales@id] <- "GN"
      PedEval$Program[PedEval$IId %in% BaseGNMales@id] <- "GN"
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
    setwd(paste0(hDir, "/GN/"))
    print(getwd())
    if (Generation == (1 + nGenerationBurn)) {
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
    blup1dat <- PedEval[,c("IId", "PhenoT1", "Program", "Generation")]
    blup1dat$Generation[blup1dat$Generation == 0] <- "00"
    write.table(blup1dat, "Blupf901.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf901.dat"))
    
    
    
    # Trait 2
    # Create phenotype file
    blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2", "Program", "Generation")]
    blup2dat$Generation[blup2dat$Generation == 0] <- "00"
    write.table(blup2dat, "Blupf902.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf902.dat"))
    
    ## set residual values in renumf90
    system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
    system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
    
    # Evaluate Trait 1
    system(command = "./renumf90 < renumParam1")
    system(command = "./blupf90 renf90.par")
    system(command = "bash Match_AFTERRenum.sh")
    system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
    
    sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
    sol1 <- sol1[sol1$V2 %in% PedEval$IId,]
    sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
    PedEval$EbvT1 <- sol1$V3
    
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
    print(summary(PedEval$EbvT1))
    print(summary(PedEval$EbvT2))
    print("Selecting GN1 animals")
    GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                          use = "ebv", trait = function(x) rowMeans(x))
    GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                          use = "ebv", trait = function(x) rowMeans(x))
    # Clean
    rm(SelCand)
    
    
    
    ##############################################################################3
    ##############################################################################3
    # ---- Program 1: PN with 100% GNmales, PN does T1 ----
    
    if (Generation == (nGenerationBurn + 1)) {
      PNFemales1 = BasePNFemales
      PedEval$Program[PedEval$IId %in% BasePNFemales@id] <- "PN1"
      PNFemales1ForPure = selectInd(pop = PNFemales1, nInd = PNFemales1@nInd * pPNFemalesPure,
                                    use = "pheno", trait = function(x) rowMeans(x, na.rm = TRUE))
    } else {
      PNFemales1ForPure = selectInd(pop = PNFemales1, nInd = PNFemales1@nInd * pPNFemalesPure,
                                    use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
    }
    

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
    SelCand@pheno[, 2] <- NA
    
    
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
    setwd(paste0(hDir, "/PN1/"))
    print("Changing dir to PN1")
    
    if (Generation == (1 + nGenerationBurn)) {
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
    blup1dat <- PedEval[,c("IId", "PhenoT1", "Program", "Generation")]
    blup1dat$Generation[blup1dat$Generation == 0] <- "00"
    write.table(blup1dat, "Blupf901.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf901.dat"))

    ## set residual values in renumf90
    system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
    system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
    
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
    blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2", "Program", "Generation")]          
    blup2dat$Generation[blup2dat$Generation == 0] <- "00"
    write.table(blup2dat, "Blupf902.dat",
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

    print(nrow(PedEval))
    print(nrow(sol1))
    print(nrow(sol2))
    print(summary(PedEval$EbvT1))
    print(summary(PedEval$EbvT2))    
    
    # Set EBVs for SelCand
    SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
    accuraciesPN1 <- rbind(accuraciesPN1, c("PN1", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
    accuraciesPN1 <- rbind(accuraciesPN1, c("PN1", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
    
    # Select
    PNMales1   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
                         use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
    PNFemales1 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                          use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
    
    # Clean
    rm(SelCand)
  }
  
  PedEval1 = PedEval
  rm(PedEval)
  DataEvalGN1 = DataEvalGN
  rm(DataEvalGN)
  
  write.table(PedEval1, paste0("PedEval1_", h2, ".csv"), quote=FALSE, row.names=FALSE)
  write.table(DataEvalGN1, paste0("DataEval1_", h2, ".csv"), quote=FALSE, row.names=FALSE)
  write.table(accuraciesPN1, paste0("Accuracies1_", h2, ".csv"), quote=FALSE, row.names=FALSE)

  ##############################################################################3
  ##############################################################################3
  # ---- PN2 ----
  DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                      Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                      GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                      GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
  DataEvalGN$Program = "GN"
  
  PedEval <- PedEvalBurnIn
  
  for (Generation in (1 + nGenerationBurn):(nGenerationEval + nGenerationBurn)) {
    if (Generation == (1 + nGenerationBurn)) {
      GNFemales = BaseGNFemales
      GNMales = BaseGNMales
      PedEval$Program[PedEval$IId %in% BaseGNFemales@id] <- "GN"
      PedEval$Program[PedEval$IId %in% BaseGNMales@id] <- "GN"
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
    setwd(paste0(hDir, "/GN/"))
    if (Generation == (1 + nGenerationBurn)) {
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
    blup1dat <-	PedEval[,c("IId", "PhenoT1", "Program", "Generation")]
    blup1dat$Generation[blup1dat$Generation == 0] <- "00"
    write.table(blup1dat, "Blupf901.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf901.dat"))



    # Trait 2
    # Create phenotype file
    blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2", "Program", "Generation")]          
    blup2dat$Generation[blup2dat$Generation == 0] <- "00"
    write.table(blup2dat, "Blupf902.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf902.dat"))    

    ## set residual values in renumf90
    system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
    system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
    
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
    if (Generation == (1 + nGenerationBurn)) {
      PNFemales2 = BasePNFemales
      PNMales2 = BaseGNMales
      PedEval$Program[PedEval$IId %in% BasePNFemales@id] <- "PN2"
      PNFemales2ForPure = selectInd(pop = PNFemales2, nInd = PNFemales2@nInd * pPNFemalesPure,
                                    use = "pheno", trait = function(x) rowMeans(x, na.rm = TRUE))
    } else {
      PNFemales2ForPure = selectInd(pop = PNFemales2, nInd = PNFemales2@nInd * pPNFemalesPure,
                                  use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
    }
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
    setwd(paste0(hDir, "/PN2/"))
    if (Generation == (1 + nGenerationBurn)) {
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
    blup1dat <-	PedEval[,c("IId", "PhenoT1", "Program", "Generation")]
    blup1dat$Generation[blup1dat$Generation == 0] <- "00"
    write.table(blup1dat, "Blupf901.dat",
                quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                na = "0", append = file.exists("Blupf901.dat"))



    ## set residual values in renumf90
    system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
    system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
    
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
    blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2", "Program", "Generation")]          
    blup2dat$Generation[blup2dat$Generation == 0] <- "00"
    write.table(blup2dat, "Blupf902.dat",
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
  write.table(PedEval2, paste0("PedEval2_", h2, ".csv"), quote=FALSE, row.names=FALSE)
  write.table(DataEvalGN2, paste0("DataEval2_", h2, ".csv"), quote=FALSE, row.names=FALSE)
  write.table(accuraciesPN2, paste0("Accuracies2_", h2, ".csv"), quote=FALSE, row.names=FALSE)
 }
}

