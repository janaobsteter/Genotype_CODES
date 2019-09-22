par <- read.csv("Documents/Projects/inProgress/Phenotyping/ParameterFile_Simulation.csv")
par$Name <- paste0(par$Reference, par$NoControl, "_", par$Pheno.Geno)
head(par)

par$MaleSelection <- NA
par$MaleSelection[par$Reference == "yes"] <- "yes"

write.table(par, "Documents/Projects/inProgress/Phenotyping/ParameterFile_Simulation.csv", quote=FALSE, row.names=FALSE, sep=",")
