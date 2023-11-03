home = paste0(getwd(), "/")

library(data.table)
setwd("./SimulatedData/")


print(paste0("I AM RUNNING THIS SCRIPT FROM ", getwd()))
args = commandArgs(trailingOnly=TRUE)
gen = as.numeric(as.character((args[1])))
rep = args[2]
scenarioHome = args[3]
scenarioImport = args[4]
trait1 = args[5]
trait2 = args[6]
strategy = args[7]

print("Reading in chip information")
chip1 <- read.table("Chip1SnpInformation.txt", header=TRUE)
chip2 <- read.table("Chip2SnpInformation.txt", header=TRUE)


chip1$pos <- paste0(chip1$ChromId, "_", chip1$SnpPhysicalPosOnChrom)
chip2$pos <- paste0(chip2$ChromId, "_", chip2$SnpPhysicalPosOnChrom)

chip2U <- !(chip2$pos %in% intersect(chip1$pos, chip2$pos))
chip1U <- (chip2$pos %in% intersect(chip1$pos, chip2$pos))

print("Reading in chip genotypes")
system("tail -n17280 ./AllIndividualsSnpChips/Chip2Genotype.txt | head -n8640 > ./AllIndividualsSnpChips/CurrentGeno1.txt")
system("tail -n8640 ./AllIndividualsSnpChips/Chip2Genotype.txt > ./AllIndividualsSnpChips/CurrentGeno2.txt")
chip2Geno1 <- fread("./AllIndividualsSnpChips/CurrentGeno1.txt", header=FALSE, fill=TRUE)
print(chip2Geno1[1:10, 1:10])
chip2Geno2 <- fread("./AllIndividualsSnpChips/CurrentGeno2.txt", header=FALSE, fill=TRUE)
print(chip2Geno2[1:10, 1:10])
chip2Geno <- rbind(chip2Geno1, chip2Geno2)
print("Merged")
print(chip2Geno[1:10, 1:10])
chip2Geno_ID <- as.data.frame(chip2Geno[,1])
chip2Geno_SNP <- as.data.frame(chip2Geno[,-1])
chip2Geno_Unique <- chip2Geno_SNP[,c(chip2U)]
print("Unique")
print(chip2Geno_Unique[1:10, 1:10])
chip1Geno_Unique <- chip2Geno_SNP[,c(chip1U)]

print("REading in QTNs")
#še za QTN
system("tail -n17280 ./UnrestrictedQtnIndivGenotypes.txt | head -n8640 > ./CurrentQTN1.txt")
system("tail -n8640 ./UnrestrictedQtnIndivGenotypes.txt > ./CurrentQTN2.txt")
qtn1 <- fread(paste0("./CurrentQTN1.txt"), header=FALSE, fill=TRUE)
qtn2 <- fread(paste0("./CurrentQTN2.txt"), header=FALSE, fill=TRUE)
print("QTNs")
print(qtn1[1:10, 1:10])
print(qtn2[1:10, 1:10])

qtn <- rbind(qtn1, qtn2)
print(qtn[1:10, 1:10])

qtn_ID <- as.data.frame(qtn[,1])
qtn_SNP <- as.data.frame(qtn[,-1])

print("Reading in alphas")
#read in average substituion effects
avgSub <- read.table("UnrestrictedQtnAverAlleleSubstEffects.txt", header=TRUE)
alpha1 <- avgSub[[paste0("AlphaNormUnres", trait1)]]
alpha2 <- avgSub[[paste0("AlphaNormUnres", trait2)]]

print("Reading in popSplit")
#SPLIT
pops <- read.csv("../PopulationSplit.txt")
library(readr)
library(reshape2)

print("Reading in empty files")
MeanHet <- read.csv(paste0(home, "MeanHet_Neutral_import.csv"))
MeanHetMarker <- read.csv(paste0(home, "MeanHet_Marker_import.csv"))
MeanHetQTN <- read.csv(paste0(home, "MeanHet_QTN_import.csv"))
genicVariance <- read.csv(paste0(home, "GenicVariance_import.csv"))

print("Calculating heterozygoisty")
for (group in c("home", "import")) {
  print(paste0("Group is", group))
  chip2Group <- chip2Geno_Unique[chip2Geno_ID$V1 %in% pops$ID[pops$Group==group],]
  print(chip2Group[1:10, 1:10])
  chip1Group <- chip1Geno_Unique[chip2Geno_ID$V1 %in% pops$ID[pops$Group==group],]
  print(chip1Group[1:10, 1:10])
  qtnGroup <- qtn_SNP[qtn_ID$V1 %in% pops$ID[pops$Group==group],]
  print("QTN group")
  print(qtnGroup[1:10, 1:10])
  print("Setting the scenario")
  scenario <- as.character(ifelse(group == "home", scenarioHome, scenarioImport))
  print(scenario)
  print(class(scenario))
  print(class(strategy))
  print(paste0("Group is ", group, ", scenario is ", scenario)) 
  print("Binding for neutral")
  MeanHet <- rbind(MeanHet, c(group, strategy, scenario, rep, gen, mean(apply(X = chip2Group, 2,  FUN = function(z) sum(z == 1)) / nrow(chip2Group))))
  print("Binding for marker")
  MeanHetMarker <- rbind(MeanHetMarker, c(group, strategy, scenario, rep, gen, mean(apply(X = chip1Group, 2,  FUN = function(z) sum(z == 1)) / nrow(chip1Group))))
  print("Binding for QTN")
  MeanHetQTN <- rbind(MeanHetQTN, c(group, strategy, scenario, rep, gen, mean(apply(X = qtnGroup, 2,  FUN = function(z) sum(z == 1)) / nrow(qtnGroup))))
  print("Computing genicVar")
  genVarF <- function(x, y) 2 * ((sum(x) / (length(x)*2))  * (( 2 * length(x) - sum(x)) / (length(x)*2))) * y^2
  print("Binding for genicVar")
  print(c(group, strategy, scenario, rep, gen, sum(mapply(genVarF, qtnGroup, alpha1)) ))
  genicVariance <- rbind(genicVariance, c(group, strategy, scenario, rep, gen, sum(mapply(genVarF, qtnGroup, alpha1)), trait1))
  genicVariance <- rbind(genicVariance, c(group, strategy, scenario, rep, gen, sum(mapply(genVarF, qtnGroup, alpha2)), trait2))
}

print("Writing files")
write.csv(MeanHet, paste0(home, "MeanHet_Neutral_import30.csv"), quote=FALSE, row.names=FALSE)
write.csv(MeanHetMarker, paste0(home, "MeanHet_Marker_import30.csv"), quote=FALSE, row.names=FALSE)
write.csv(MeanHetQTN, paste0(home, "MeanHet_QTN_import30.csv"), quote=FALSE, row.names=FALSE)
write.csv(genicVariance, paste0(home, "GenicVariance_import30.csv"), quote=FALSE, row.names=FALSE)


