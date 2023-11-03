args = commandArgs(trailingOnly=TRUE)
date = args[1]

acc <- data.frame()
for (rep in 0:9) {
  for (control in c(11)) {
    for (scenario in "Class") {
	    accTmp <- read.csv(paste0("10K/SU55_permEnv/", scenario,  rep, "_", control, "/Accuracies_CatAge.csv"))
	    accTmp <- accTmp[accTmp$AgeCat %in% c("potomciNP0", "telF1", "vhlevljeni1", "k3", "k5", "pBM4", "pBM5", "pBM6", "pb6", "pb7", "pb8", "pb9", "cak5"),]
	    accTmp$Rep <- rep
	    accTmp$NoControl <- control
	    acc <- rbind(acc, accTmp)
    }
  }
}

write.table(acc, paste0("CombinedAcc_phenoClass_", date, ".csv"), quote=FALSE, row.names=FALSE)
