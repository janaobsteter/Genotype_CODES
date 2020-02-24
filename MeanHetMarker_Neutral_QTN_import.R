#setwd("~/bin/AlphaSim1.05Linux/SimulatedData/")
home = paste0(getwd(), "/")

library(readr)
setwd("./SimulatedData/")


print(paste0("I AM RUNNING THIS SCRIPT FROM ", getwd()))
args = commandArgs(trailingOnly=TRUE)
gen = as.numeric(as.character((args[1])))
rep = args[2]
scenario = args[3]
strategy = args[4]

print("Reading in chip information")
chip1 <- read.table("Chip1SnpInformation.txt", header=TRUE)
chip2 <- read.table("Chip2SnpInformation.txt", header=TRUE)


chip1$pos <- paste0(chip1$ChromId, "_", chip1$SnpPhysicalPosOnChrom)
chip2$pos <- paste0(chip2$ChromId, "_", chip2$SnpPhysicalPosOnChrom)

chip2U <- !(chip2$pos %in% intersect(chip1$pos, chip2$pos))
chip1U <- (chip2$pos %in% intersect(chip1$pos, chip2$pos))

print("Reading in chip genotypes")
chip2Geno <- readr::read_table("./AllIndividualsSnpChips/Chip2Genotype.txt", skip=8640*2*(gen-1),
                               col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
chip2Geno_ID <- as.data.frame(chip2Geno[,1])
chip2Geno_SNP <- as.data.frame(chip2Geno[,-1])
chip2Geno_Unique <- chip2Geno_SNP[,c(chip2U)]
chip1Geno_Unique <- chip2Geno_SNP[,c(chip1U)]

print("REading in QTNs")
#še za QTN
qtn <- readr::read_table(paste0("./UnrestrictedQtnIndivGenotypes.txt"), skip=8640*2*(gen-1),
                       col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
qtn_ID <- as.data.frame(qtn[,1])
qtn_SNP <- as.data.frame(qtn[,-1])

print("Reading in alphas")
#read in average substituion effects
avgSub <- read.table("UnrestrictedQtnAverAlleleSubstEffects.txt", header=TRUE)
alpha <- avgSub$AlphaNormUnres1

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
  chip2Group <- chip2Geno_Unique[chip2Geno_ID$X1 %in% pops$ID[pops$Group==group],]
  chip1Group <- chip1Geno_Unique[chip2Geno_ID$X1 %in% pops$ID[pops$Group==group],]
  qtnGroup <- qtn_SNP[qtn_ID$X1 %in% pops$ID[pops$Group==group],]
  MeanHet <- rbind(MeanHet, c(group, strategy, scenario, rep, gen, mean(apply(X = chip2Group, 2,  FUN = function(z) sum(z == 1)) / nrow(chip2Group))))
  MeanHetMarker <- rbind(MeanHetMarker, c(group, strategy, scenario, rep, gen, mean(apply(X = chip1Group, 2,  FUN = function(z) sum(z == 1)) / nrow(chip1Group))))
  MeanHetQTN <- rbind(MeanHetQTN, c(group, strategy, scenario, rep, gen, mean(apply(X = qtnGroup, 2,  FUN = function(z) sum(z == 1)) / nrow(qtnGroup))))
  genVarF <- function(x, y) 2 * ((sum(x) / (length(x)*2))  * (( 2 * length(x) - sum(x)) / (length(x)*2))) * y^2
  genicVariance <- rbind(genicVariance, c(group, strategy, scenario, rep, gen, sum(mapply(genVarF, qtnGroup, alpha)) ))
}

print("Writing files")
write.csv(MeanHet, paste0(home, "MeanHet_Neutral_import.csv"), quote=FALSE, row.names=FALSE)
write.csv(MeanHetMarker, paste0(home, "MeanHet_Marker_import.csv"), quote=FALSE, row.names=FALSE)
write.csv(MeanHetQTN, paste0(home, "MeanHet_QTN_import.csv"), quote=FALSE, row.names=FALSE)
write.csv(genicVariance, paste0(home, "GenicVariance_import.csv"), quote=FALSE, row.names=FALSE)


