library(AlphaSimR)
nMales = 100
nFemales = 10000


nGenerationBurn = 20
nGenerationEval = 20

# ---- Base population genomes ----

founderPop = runMacs(nInd = nMales + nFemales,
                     nChr = 1,
                     segSites = 1000,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")

save.image(file = "~/KISdir/Documents/PhD/Projects/inProgress/AlphaPart/SimulationAccBlup.RData")
#save(founderPop, file = "~/KISdir/Documents/PhD/Projects/inProgress/AlphaPart/FounderPop.RData")

load("~/Documents/PhD/Projects/inProgress/GenomicAlphaPart/FounderPopObject.RData")

SP = SimParam$new(founderPop)
VarA = matrix(data = rep(1.0, 100), nrow = 10); cov2cor(VarA)
VarE = matrix(data = rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 10), nrow = 10)
diag(VarE) = 3.0
cov2cor(VarE)
VarP = VarA + VarE; diag(VarA) / diag(VarP)
SP$addTraitA(nQtlPerChr = 10, mean = rep(0, 10), var = diag(VarA), cor = cov2cor(VarA))

SP$setGender(gender = "yes_sys")

Base = newPop(founderPop)
Base = setPheno(Base, varE = VarE)
# Select males and females
Males = selectInd(Base, gender = "M", nInd = nMales, 
                  use = "pheno", trait = function(x) rowMeans(scale(x)))##Base[Base@gender == "M"]
Females = Base[!(Base@id %in% Males@id)]
Females@gender = "F"

# Create selection candidates
selCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd, nProgeny = 20)
# Set Phenotype for all 10 traits
selCand = setPheno(selCand, varE = VarE)
table(selCand@gender)
table(selCand@father[selCand@gender == "F"])
# Select the daughters
Daughters = selCand[selCand@gender == "F"]
Daughters@nInd
min(as.numeric(Daughters@id))
max(as.numeric(Daughters@id))

# -- Estimate EBVs --
options(width=200)

### --- Required packages ---

## install.packages(pkg=c("pedigreemm", "MatrixModels"))

library(package="pedigreemm")   ## pedigree functions
library(package="MatrixModels") ## sparse matrices
library(package="plyr") 

### --- Data ---
# Obtain Fathers and Daughters from the simulated population
#we have money for 10,000 phenotypes
#10,000 daughters (100 per bull), all phenotyped once
# allData <- data.frame(
#   individual=c(Males@id, daughters@id),
#   father=c(Males@father, daughters@father),
#   phenotype=c(rep(NA, Males@nInd), daughters@pheno[,1]))


###################################################
#1,000 daughters - constant 100 per bull, vary the number of phenotypic records
accuracies1 <- data.frame(nRecords = NA, nDaughtersPerSire = NA, AccSires = NA, AccDaughters = NA)

for (nRecords in c(1, 2, 5, 10)) {
  print(paste0("Number of records is ", nRecords))
  nDaughtersPerBull <- 100
  PED <- data.frame(daughter = Daughters@id, father = Daughters@father)
  SelDaughters <- ddply(PED,.(father),function(x) x[sample(nrow(x),nDaughtersPerBull),])
  daughters <- Daughters[Daughters@id %in% SelDaughters$daughter]
  allData <- data.frame(
    individual=c(Males@id, rep(daughters@id, nRecords)),
    father=c(Males@father, rep(daughters@father, nRecords)),
    phenotype=c(rep(NA, Males@nInd), daughters@pheno[,1:nRecords]))
  ############################################
  
  #order
  allData <- allData[order(as.numeric(as.character(allData$individual))),]
  head(allData)
  tail(head(allData, 100), 10)
  tail(allData)
  table(allData$father)
  
  #allData$father <- as.numeric(allData$father)
  allData$mother = NA
  head(allData)
  table(allData$father)
  allData$father[allData$father == 0] <- NA
  #allData$father <- as.numeric(as.character(allData$father))
  allData$father <- as.factor(allData$father)
  allData$mother <- as.factor(allData$mother)
  allData$individual <- as.factor(allData$individual)
  head(allData)
  tail(head(allData, 100), 10)
  table(allData$father)
  summary(allData$phenotype)
  allData <- allData[,c('individual', 'father', 'mother', 'phenotype')]
  sum(allData$father %in% allData$individual)
  
  
  
  ## Variance components
  sigma2e <- VarE
  sigma2a <- VarA
  (h2 <- sigma2a / (sigma2a + sigma2e))
  
  ### --- Setup data ---
  
  ## Make sure each individual has only one record in pedigree
  ped <- allData[!duplicated(allData$individual), 1:3]
  head(ped)
  sum(ped$father %in% ped$individual)
  
  ## Factors (this eases buliding the design matrix considerably)
  allData$individual <- factor(allData$individual)
  #allData$group      <- factor(allData$group)
  
  
  #dat <- allData[!is.na(allData$phenotype), ]
  dat = allData

  
  ### --- Setup MME ---
  
  ## Phenotype vector
  (y <- dat$phenotype)
  
  ## Design matrix for the "fixed" effects (group)
  #(X <- model.Matrix( ~ group,          data=dat, sparse=TRUE))
  
  ## Design matrix for the "random" effects (individual)
  # library(lme4)
  # mform <- dat$phenotype ~ (1 | dat$individual)
  # (bar.f <- findbars(mform)) # list with 3 terms
  # mf <- model.frame(subbars(mform),data=dat)
  # rt <- mkReTrms(bar.f,mf)
  # names(rt)
  # (Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))
  
  
  ## Variance components
  sigma2e <- VarE[1,1]
  sigma2a <- VarA[1,1]
  (h2 <- sigma2a / (sigma2a + sigma2e))
  
  ### --- Setup data ---
  
  ## Make sure each individual has only one record in pedigree
  ped <- allData[!duplicated(allData$individual), 1:3]
  
  ## Factors (this eases buliding the design matrix considerably)
  allData$individual <- factor(allData$individual)
  #example$group      <- factor(example$group)
  
  ## Phenotype data
  dat <- allData[!is.na(allData$phenotype), ]
  
  ### --- Setup MME ---
  
  ## Phenotype vector
  (y <- dat$phenotype)
  
  ## Design matrix for the "fixed" effects (group)
  #(X <- model.Matrix( ~ group,          data=dat, sparse=TRUE))
  
  ## Design matrix for the "random" effects (individual)
  (Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))
  order <- paste0("individual", ped$individual)
  Z <- Z[,order]
  
  ## Inverse additive relationship matrix
  ped$individual <- as.numeric(as.character(ped$individual))
  ped$father <- as.numeric(as.character(ped$father))
  ped2 <- pedigree(sire=ped$father, dam=ped$mother, label=ped$individual)
  TInv <- as(ped2, "sparseMatrix")
  DInv <- Diagonal(x=1/Dmat(ped2))
  AInv <- crossprod(sqrt(DInv) %*% TInv)
  
  ## Variance ratio
  alpha <- sigma2e / sigma2a
  
  ## Mixed Model Equations (MME)
  
  ## ... Left-Hand Side (LHS) without pedigree prior
  #(LHS0 <- rBind(cBind(crossprod(X),    crossprod(X, Z)),
   #              cBind(crossprod(Z, X), crossprod(Z, Z))))
  
  ## ... Left-Hand Side (LHS) with    pedigree prior
  round(LHS <- crossprod(Z, Z) + AInv * alpha, digits=1)
  
  ## ... Right-Hand Side (RHS)
  (RHS <-       crossprod(Z, y))
  
  ### --- Solutions ---
  
  ## Solve
  LHSInv <- solve(LHS)
  sol <- LHSInv %*% RHS
  
  ## Standard errors
  se <- diag(LHSInv) * sigma2e
  
  ## Reliabilities
  r2 <- 1 - diag(LHSInv) * alpha
  
  ## Accuracies
  r <- sqrt(r2)
  acc <- data.frame(id = 1:length(r), r = r)
  qplot(data = acc, x = id, y = r)
  hist(r)
  
  # (Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))
  # order <- paste0("individual", ped$individual)
  # #Z <- Z[,order]
  # #colnames(Z)[1:100]
  # #colnames(Z)[1000:1100]
  # 
  # ## Inverse additive relationship matrix
  # ped$father <- as.numeric(as.character(ped$father))
  # ped$individual <- as.numeric(as.character(ped$individual))
  # table(ped$father)
  # head(ped)
  # ped2 <- pedigree(sire=ped$father, dam=ped$mother, label=ped$individual)
  # head(ped)
  # TInv <- as(ped2, "sparseMatrix")
  # DInv <- Diagonal(x=1/Dmat(ped2))
  # AInv <- crossprod(sqrt(DInv) %*% TInv)
  # 
  # ## Variance ratio
  # alpha <- diag(sigma2e / sigma2a)[1]
  # 
  # ## Mixed Model Equations (MME)
  # 
  # ## ... Left-Hand Side (LHS) without pedigree prior
  # # (LHS0 <- rBind(cBind(crossprod(X),    crossprod(X, Z)),
  # #                cBind(crossprod(Z, X), crossprod(Z, Z))))
  # 
  # ## ... Left-Hand Side (LHS) with    pedigree prior
  # round(
  #   LHS <- crossprod(Z, Z) + AInv * alpha, digits=1)
  # 
  # ## ... Right-Hand Side (RHS)
  # RHS <- crossprod(Z, y)
  # 
  # ### --- Solutions ---
  # 
  # ## Solve
  # LHSInv <- solve(LHS)
  # sol <- LHSInv %*% RHS
  # 
  # ## Standard errors
  # se <- diag(LHSInv) * sigma2e[1,1]
  # 
  # ## Reliabilities
  # #r2 <- 1 - (se / VarA[1,1])
  # r2a <- 1 - diag(LHSInv) * alpha
  # 
  # ## Accuracies
  # r <- sqrt(r2)
  # hist(r)
  # ## For fathers
  accSire <- mean(r[1:100])
  # r[1:100]
  # ## For daughters
  accDaughters <- mean(r[101:length(r)])
  # r[101:10100]
  # hist(r[101:length(r)])

  accuracies1 <- rbind(accuracies1, c(nRecords, nDaughtersPerBull, accSire, accDaughters))
}

accuracies1$PA <- sqrt(0.25*(accuracies1$AccSires**2 + accuracies1$AccDaughters**2))
accuracies1 <- accuracies1[-1,]
accuracies1
###################################################
# funds for 10,000 phenotyping, vary the number of daughters with phenotype
accuracies <- data.frame(nRecords = NA, nDaughtersPerSire = NA, AccSires = NA, AccDaughters = NA)

for (nRecords in c(1, 2, 5, 10)) {
  print(paste0("Number of records is ", nRecords))
  nDaughtersPerBull <- 100000/100/nRecords
  PED <- data.frame(daughter = Daughters@id, father = Daughters@father)
  SelDaughters <- ddply(PED,.(father),function(x) x[sample(nrow(x),nDaughtersPerBull),])
  daughters <- Daughters[Daughters@id %in% SelDaughters$daughter]
  allData <- data.frame(
    individual=c(Males@id, rep(daughters@id, nRecords)),
    father=c(Males@father, rep(daughters@father, nRecords)),
    phenotype=c(rep(NA, Males@nInd), daughters@pheno[,1:nRecords]))
  ############################################
  
  #order
  allData <- allData[order(as.numeric(as.character(allData$individual))),]
  head(allData)
  tail(head(allData, 100), 10)
  tail(allData)
  table(allData$father)
  
  #allData$father <- as.numeric(allData$father)
  allData$mother = NA
  head(allData)
  table(allData$father)
  allData$father[allData$father == 0] <- NA
  #allData$father <- as.numeric(as.character(allData$father))
  allData$father <- as.factor(allData$father)
  allData$mother <- as.factor(allData$mother)
  allData$individual <- as.factor(allData$individual)
  head(allData)
  tail(head(allData, 100), 10)
  table(allData$father)
  summary(allData$phenotype)
  allData <- allData[,c('individual', 'father', 'mother', 'phenotype')]
  sum(allData$father %in% allData$individual)
  
  
  
  ## Variance components
  sigma2e <- VarE
  sigma2a <- VarA
  (h2 <- sigma2a / (sigma2a + sigma2e))
  
  ### --- Setup data ---
  
  ## Make sure each individual has only one record in pedigree
  ped <- allData[!duplicated(allData$individual), 1:3]
  head(ped)
  sum(ped$father %in% ped$individual)
  
  ## Factors (this eases buliding the design matrix considerably)
  allData$individual <- factor(allData$individual)
  #allData$group      <- factor(allData$group)
  
  
  #dat <- allData[!is.na(allData$phenotype), ]
  dat = allData

  
  ### --- Setup MME ---
  
  ## Phenotype vector
  (y <- dat$phenotype)
  
  ## Design matrix for the "fixed" effects (group)
  #(X <- model.Matrix( ~ group,          data=dat, sparse=TRUE))
  
  ## Design matrix for the "random" effects (individual)
  # library(lme4)
  # mform <- dat$phenotype ~ (1 | dat$individual)
  # (bar.f <- findbars(mform)) # list with 3 terms
  # mf <- model.frame(subbars(mform),data=dat)
  # rt <- mkReTrms(bar.f,mf)
  # names(rt)
  # (Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))
  
  
  ## Variance components
  sigma2e <- VarE[1,1]
  sigma2a <- VarA[1,1]
  (h2 <- sigma2a / (sigma2a + sigma2e))
  
  ### --- Setup data ---
  
  ## Make sure each individual has only one record in pedigree
  ped <- allData[!duplicated(allData$individual), 1:3]
  
  ## Factors (this eases buliding the design matrix considerably)
  allData$individual <- factor(allData$individual)
  #example$group      <- factor(example$group)
  
  ## Phenotype data
  dat <- allData[!is.na(allData$phenotype), ]
  
  ### --- Setup MME ---
  
  ## Phenotype vector
  (y <- dat$phenotype)
  
  ## Design matrix for the "fixed" effects (group)
  #(X <- model.Matrix( ~ group,          data=dat, sparse=TRUE))
  
  ## Design matrix for the "random" effects (individual)
  (Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))
  order <- paste0("individual", ped$individual)
  Z <- Z[,order]
  
  ## Inverse additive relationship matrix
  ped$individual <- as.numeric(as.character(ped$individual))
  ped$father <- as.numeric(as.character(ped$father))
  ped2 <- pedigree(sire=ped$father, dam=ped$mother, label=ped$individual)
  TInv <- as(ped2, "sparseMatrix")
  DInv <- Diagonal(x=1/Dmat(ped2))
  AInv <- crossprod(sqrt(DInv) %*% TInv)
  
  ## Variance ratio
  alpha <- sigma2e / sigma2a
  
  ## Mixed Model Equations (MME)
  
  ## ... Left-Hand Side (LHS) without pedigree prior
  #(LHS0 <- rBind(cBind(crossprod(X),    crossprod(X, Z)),
   #              cBind(crossprod(Z, X), crossprod(Z, Z))))
  
  ## ... Left-Hand Side (LHS) with    pedigree prior
  round(LHS <- crossprod(Z, Z) + AInv * alpha, digits=1)
  
  ## ... Right-Hand Side (RHS)
  (RHS <-       crossprod(Z, y))
  
  ### --- Solutions ---
  
  ## Solve
  LHSInv <- solve(LHS)
  sol <- LHSInv %*% RHS
  
  ## Standard errors
  se <- diag(LHSInv) * sigma2e
  
  ## Reliabilities
  r2 <- 1 - diag(LHSInv) * alpha
  
  ## Accuracies
  r <- sqrt(r2)
  acc <- data.frame(id = 1:length(r), r = r)
  qplot(data = acc, x = id, y = r)
  #hist(r)
  
  # (Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))
  # order <- paste0("individual", ped$individual)
  # #Z <- Z[,order]
  # #colnames(Z)[1:100]
  # #colnames(Z)[1000:1100]
  # 
  # ## Inverse additive relationship matrix
  # ped$father <- as.numeric(as.character(ped$father))
  # ped$individual <- as.numeric(as.character(ped$individual))
  # table(ped$father)
  # head(ped)
  # ped2 <- pedigree(sire=ped$father, dam=ped$mother, label=ped$individual)
  # head(ped)
  # TInv <- as(ped2, "sparseMatrix")
  # DInv <- Diagonal(x=1/Dmat(ped2))
  # AInv <- crossprod(sqrt(DInv) %*% TInv)
  # 
  # ## Variance ratio
  # alpha <- diag(sigma2e / sigma2a)[1]
  # 
  # ## Mixed Model Equations (MME)
  # 
  # ## ... Left-Hand Side (LHS) without pedigree prior
  # # (LHS0 <- rBind(cBind(crossprod(X),    crossprod(X, Z)),
  # #                cBind(crossprod(Z, X), crossprod(Z, Z))))
  # 
  # ## ... Left-Hand Side (LHS) with    pedigree prior
  # round(
  #   LHS <- crossprod(Z, Z) + AInv * alpha, digits=1)
  # 
  # ## ... Right-Hand Side (RHS)
  # RHS <- crossprod(Z, y)
  # 
  # ### --- Solutions ---
  # 
  # ## Solve
  # LHSInv <- solve(LHS)
  # sol <- LHSInv %*% RHS
  # 
  # ## Standard errors
  # se <- diag(LHSInv) * sigma2e[1,1]
  # 
  # ## Reliabilities
  # #r2 <- 1 - (se / VarA[1,1])
  # r2a <- 1 - diag(LHSInv) * alpha
  # 
  # ## Accuracies
  # r <- sqrt(r2)
  # hist(r)
  # ## For fathers
  accSire <- mean(r[1:100])
  # r[1:100]
  # ## For daughters
  accDaughters <- mean(r[101:length(r)])
  # r[101:10100]
  # hist(r[101:length(r)])

  accuracies <- rbind(accuracies, c(nRecords, nDaughtersPerBull, accSire, accDaughters))
}

accuracies
accuracies$Total <- accuracies$nDaughtersPerSire * 100
accuracies$TotalPheno <- accuracies$Total * accuracies$nRecords
accuracies$PA <- sqrt(0.25*(accuracies$AccSires**2 + accuracies$AccDaughters**2))
accuracies <- accuracies[-1,]
accuracies

ped <- data.frame(individual = c(1,1,1,2,2,3,3,3,4,5,5), father = c(0, 0, 0, 1, 1, 1, 1, 1, 2, 4, 4), mother = rep(NA, 11), phenotype = rnorm(11))
ped$individual <- as.factor(ped$individual)
ped$father <- as.factor(ped$father)
ped$mother <- as.factor(ped$mother)
(Z <- model.Matrix(~ individual - 1, data=ped, sparse=TRUE))
