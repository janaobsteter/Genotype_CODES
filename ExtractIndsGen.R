args = commandArgs(trailingOnly=TRUE)
gen <- args[1]

ped <- read.table("PedigreeAndGeneticValues.txt", header=TRUE)
fathers <- unique(ped$Father[ped$Generation == gen])
indivs <- ped$Indiv[ped$Generation == gen]

A <- as.data.frame(c(fathers, indivs))

write.table(A, paste0("Indivs", gen, ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
