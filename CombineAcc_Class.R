acc <- data.frame()
for (ref in c( "10K")) {
  for (rep in c(0:9)) {
	for (pheno in c(1, 2, 5, 9, 11)) {
	    accTmp <- read.csv(paste0(ref, "/SU55_permEnv/Class", rep, "_", pheno, "/Accuracies_CatAge.csv"))
	    accTmp <- accTmp[accTmp$AgeCat %in% c("potomciNP0", "telF1", "genTest1", "vhlevljeni1", "cak5", "k3", "k5"),]
	    accTmp$Rep <- rep
	    accTmp$Ref <- ref
	    accTmp$sc <- pheno
	    acc <- rbind(acc, accTmp)
	}
  }
}

write.table(acc, "CombinedAcc_Class_15122019.csv", quote=FALSE, row.names=FALSE)
