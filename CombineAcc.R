acc <- data.frame()
for (ref in c("10K")) {
  for (rep in c(0, 2)) {
    for (varE in c(0.1, 0.2, 0.5, 0.6, 1.0, 1.5, 1.8)) {
    accTmp <- read.csv(paste0(ref, "/SU55_permEnv/Gen", rep, "_11_", sprintf("%.1f", varE), "/Accuracies_CatAge.csv"))
    accTmp <- accTmp[accTmp$AgeCat %in% c("potomciNP0", "telF1", "genTest1", "k3", "k5", "pBM4", "pBM5", "pBM6", "gpb2", "gpb3", "gpb4", "gpb5"),]
    accTmp$Rep <- rep
    accTmp$Ref <- ref
    accTmp$varE <- varE
    acc <- rbind(acc, accTmp)


    }
  }
}

write.table(acc, "CombinedAcc.csv", quote=FALSE, row.names=FALSE)
