args = commandArgs(trailingOnly=TRUE)
date = args[1]

acc <- data.frame()
for (rep in 0:9) {
  for (control in c(1,2,5,8,9,10)) {
    for (pg in c("1_1", "1_2", "2_1")) {
      for (ref in c("False")) {
	    accTmp <- read.csv(paste0("10K/SU55_permEnv/Gen", rep, "_", ref, control, "_", pg, "/Accuracies_CatAge.csv"))
	    accTmp <- accTmp[accTmp$AgeCat %in% c("potomciNP0", "telF1", "genTest1", "k2", "k3", "k4", "k5", "k6", "pBM4", "pBM5", "pBM6", "bm7", "gpb2", "gpb3", "gpb4", "gpb5", "gpb6"),]
	    accTmp$Rep <- rep
	    accTmp$NoControl <- control
	    accTmp$gp <- pg
	    acc <- rbind(acc, accTmp)
	}
    }
  }
}

write.table(acc, paste0("CombinedAcc_phenoFalse_", date, ".csv"), quote=FALSE, row.names=FALSE)
