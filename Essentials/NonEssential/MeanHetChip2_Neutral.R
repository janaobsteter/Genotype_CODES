#setwd("~/bin/AlphaSim1.05Linux/SimulatedData/")
home = paste0(getwd(), "/")

library(readr)
setwd("./SimulatedData/")

args = commandArgs(trailingOnly=TRUE)
gen = as.numeric(as.character((args[1])))
rep = args[2]
scenario = args[3]
strategy = args[4]

chip1 <- read.table("Chip1SnpInformation.txt", header=TRUE)
chip2 <- read.table("Chip2SnpInformation.txt", header=TRUE)

chip1$pos <- paste0(chip1$ChromId, "_", chip1$SnpPhysicalPosOnChrom)
chip2$pos <- paste0(chip2$ChromId, "_", chip2$SnpPhysicalPosOnChrom)

chip2U <- !(chip2$pos %in% intersect(chip1$pos, chip2$pos))

chip2Geno <- readr::read_table("./AllIndividualsSnpChips/Chip2Genotype.txt", skip=8640*(gen-1),
                               col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
chip2Geno_ID <- as.data.frame(chip2Geno[,1])
chip2Geno_SNP <- as.data.frame(chip2Geno[,-1])
chip2Geno_Unique <- chip2Geno_SNP[,c(chip2U)]


library(readr)
library(reshape2)


MeanHet <- read.csv(paste0(home, "MeanHet_Neutral.csv"))

MeanHet <- rbind(MeanHet, c(strategy, scenario, rep, gen, mean(apply(X = chip2Geno_Unique, 2,  FUN = function(z) sum(z == 1)) / ncol(chip2Geno_Unique))))

write.csv(MeanHet, paste0(home, "MeanHet_Neutral.csv"), quote=FALSE, row.names=FALSE)

