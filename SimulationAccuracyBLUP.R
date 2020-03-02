library(AlphaSimR)
nMales = 10
nFemales = 1000
#nFemales = 100


nGenerationBurn = 20
nGenerationEval = 20

# ---- Base population genomes ----

founderPop = runMacs(nInd = nMales + nFemales,
                     nChr = 1,
                     segSites = 10,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")



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
selCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd, nProgeny = 1)
# Set Phenotype for all 10 traits
selCand = setPheno(selCand, varE = VarE)
table(selCand@gender)
table(selCand@father[selCand@gender == "F"])
# Select the daughters
daughters = selCand[selCand@gender == "F"]


# -- Estimate EBVs --
options(width=200)

### --- Required packages ---

 install.packages(pkg=c("pedigreemm", "MatrixModels"))

library(package="pedigreemm")   ## pedigree functions
library(package="MatrixModels") ## sparse matrices

### --- Data ---
# Obtain Fathers and Daughters from the simulated population
#we have money for 10,000 phenotypes
#10,000 daughters (100 per bull), all phenotyped once
allData <- data.frame(
  individual=c(Males@id, daughters@id),
  father=c(Males@father, daughters@father),
  phenotype=c(rep(NA, Males@nInd), daughters@pheno[,1]))

#1,000 daughters (10 per bull), all phenotyped 10-times
nDaughtersPerBull = 10
nRecords <- 10
SelDaughters <- ddply(allData,.(father),function(x) x[sample(nrow(x),nDaughtersPerBull),])
daughters <- daughters[daughters@id %in% SelDaughters$individual]
allData <- data.frame(
  individual=c(Males@id, rep(daughters@id, nRecords)),
  father=c(Males@father, rep(daughters@father, nRecords)),
  phenotype=c(rep(NA, Males@nInd), daughters@pheno[,1:nRecords]))

#order
allData <- allData[order(as.numeric(as.character(allData$individual))),]
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


dat <- allData[!is.na(allData$phenotype), ]
dat

### --- Setup MME ---

## Phenotype vector
(y <- dat$phenotype)

## Design matrix for the "fixed" effects (group)
#(X <- model.Matrix( ~ group,          data=dat, sparse=TRUE))

## Design matrix for the "random" effects (individual)
(Z <- model.Matrix(~ individual - 1, data=dat, sparse=TRUE))

## Inverse additive relationship matrix
ped$father <- as.numeric(as.character(ped$father))
ped$individual <- as.numeric(as.character(ped$individual))
table(ped$father)
ped2 <- pedigree(sire=ped$father, dam=ped$mother, label=ped$individual)
head(ped)
TInv <- as(ped2, "sparseMatrix")
DInv <- Diagonal(x=1/Dmat(ped2))
AInv <- crossprod(sqrt(DInv) %*% TInv)

## Variance ratio
alpha <- diag(sigma2e / sigma2a)[1]

## Mixed Model Equations (MME)

## ... Left-Hand Side (LHS) without pedigree prior
# (LHS0 <- rBind(cBind(crossprod(X),    crossprod(X, Z)),
#                cBind(crossprod(Z, X), crossprod(Z, Z))))

## ... Left-Hand Side (LHS) with    pedigree prior
round(
  LHS <- crossprod(Z, Z) + AInv * alpha, digits=1)

## ... Right-Hand Side (RHS)
RHS <- crossprod(Z, y)

### --- Solutions ---

## Solve
LHSInv <- solve(LHS)
sol <- LHSInv %*% RHS

## Standard errors
se <- diag(LHSInv) * sigma2e[1,1]

## Reliabilities
r2 <- 1 - diag(LHSInv) * alpha

## Accuracies
r <- sqrt(r2)
r
