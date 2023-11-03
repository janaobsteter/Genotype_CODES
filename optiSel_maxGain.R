
# ---- Code requirements ----

#install.packages("tidyverse")
#install.packages("optiSel")
#install.packages("mipfp")
library(package = "tidyverse") # for tidy data handling
library(package = "optiSel") # for Optimum Contribution Selection via quadratic programming
#library(package = "mipfp") # for Iterative proportional fitting (used in MateAtRandom)
source(file = "Functions.R")

args <- commandArgs(TRUE)

DEGREE <- as.numeric(args[1])
# ---- Data ----

#setwd("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample/")
# Breeding values
Data = read.table(file="CRITERION.txt", sep=",")
colnames(Data) = c("Indiv", "Ebv1")
# ... standardize them
Tmp = scale(Data$Ebv1)
EbvMean = attr(Tmp, "scaled:center")
EbvSd   = attr(Tmp, "scaled:scale")
Data$Ebv1S = Tmp[, 1]

# Add gender
Data$Sex = c("male", "female")[read.table(file="GENDER.txt", sep=" ")[, 2]]

#femaleSample <- sample(Data$Indiv[Data$Sex=="female"], 500, replace=FALSE)
#maleSample <- sample(Data$Indiv[Data$Sex=="male"], 20, replace=FALSE)

#Data <- Data[Data$Indiv %in% c(femaleSample, maleSample),]

# Relationship matrix
Nrm = as.matrix(read.table(file="PedigreeNrm.txt")[, -1])
nInd = nrow(Data)
dimnames(Nrm) = list(Data$Indiv, Data$Indiv)

#Nrm <- Nrm[rownames(Nrm) %in% c(maleSample, femaleSample),]
#Nrm <- Nrm[,colnames(Nrm) %in% c(maleSample, femaleSample)]
#nInd = nrow(Data)
# ... make it coancestry
Kin = Nrm / 2

str(Data); head(Data)
str(Nrm); Nrm[1:10, 1:10]

# Define maximum possible contributions
#MaxUse = rep(4 / (2 * nInd), times=nInd)
MaxUse = ifelse(Data$Sex == "female", NA, 2000/(2*nInd))

names(MaxUse) = Data$Indiv
MaxUseM <- MaxUse[!is.na(MaxUse)]

# Construct candidates object
Candidates = candes(phen=Data, sKin=Kin, N=nInd)

# ---- Optimisation for maximal selection criterion ----

ConstraintsMaxSelCriterion = list(ub=MaxUse, uniform="female")
ContributionsMaxSelCriterion = opticont(method="max.Ebv1", cand=Candidates, con=ConstraintsMaxSelCriterion)
ContributionsMaxSelCriterion$info
ContributionsMaxSelCriterion$mean
ContributionsMaxSelCriterion$parent$nOff = noffspring(ContributionsMaxSelCriterion$parent, N=8640)$nOff 
ContributionsMaxSelCriterion$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))



fathers <- ContributionsMaxSelCriterion$parent[(ContributionsMaxSelCriterion$parent$Sex=="male") & (ContributionsMaxSelCriterion$parent$nOff != 0), c("Indiv", "nOff")]
Ocetje <- sample(rep(fathers$Indiv, fathers$nOff))

write.table(Ocetje, "Ocetje.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
