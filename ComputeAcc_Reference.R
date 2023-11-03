COR = data.frame(Rep = NA, Ref=NA, Pop=NA, Cor=NA, nInd=NA)

for (rep in 0:2) {
  ped = read.table(paste0("FillInBurnIn", rep, "_permEnv/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE)
  ped$Cat <- paste0(ped$cat, ped$age)
  for (ref in c("1K", "5K")) { ##: ##, "3K", "4K", "5K", "6K", "7K", "8K", "9K")) {
    ind <- read.csv(paste0("FillInBurnIn", rep, "_permEnv/IndForGeno_", ref, ".txt"))
    nind <- nrow(ind)
#    sol <- read.table(paste0("FillInBurnIn", rep, "/renumbered_Solutions_", ref))[,2:3]
    sol <- read.table(paste0("FillInBurnIn", rep,"_permEnv/renumbered_Solutions_", ref, "_newModel"))[,2:3]

    colnames(sol) <- c("Indiv", paste0("EBV_", ref))
    ped <- merge(ped, sol, by="Indiv")
    print(nind)
    pedC <- ped[ped$cat == "potomciNP",]
    COR <- rbind(COR, c(rep, ref, "potomciNP", cor(pedC[[paste0("EBV_", ref)]], pedC$gvNormUnres1), nind))
  }
}

write.table(COR, "Correlation_reference_newModel.csv", quote=FALSE, row.names=FALSE)
