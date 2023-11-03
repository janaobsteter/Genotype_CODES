args = commandArgs(trailingOnly=TRUE)
gen <- args[1]

ped <- read.table("PedigreeAndGeneticValues.txt", header=TRUE)
pedGen <- ped[ped$Generation == gen,]
pedGen$OF <- NA

g <- read.table(paste0("Genotypes", gen, ".txt"))


for (row in 1:8640) {
  fG <- as.data.frame(t(g[g$V1 == pedGen$Father[row],]))
  fO <- as.data.frame(t(g[g$V1 == pedGen$Indiv[row],]))
  G <- as.data.frame(cbind(fO[-1,], fG[-1,]))
  colnames(G) <- c("O","F")
  pedGen$OF[row] <- sum((G$O == 0 & G$F == 2) | (G$O == 2 & G$F == 0))
  
}


write.table(pedGen, paste0("Ped", gen, ".txt"), quote=FALSE, row.names=FALSE)
