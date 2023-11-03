acc <- data.frame(Rep = NA, varE = NA, Ref = NA, Acc = NA)
for (ref in c("1K", "10K")) {
  for (rep in c(0, 2)) {
    for (varE in c(1.5, 1.8)) {
	    genPed <- read.csv(paste0(ref, "/SU55_permEnv/Gen", rep, "_11_", varE, "/GenPed_EBV.txt"))
	    Ped <- read.table(paste0(ref, "/SU55_permEnv/Gen", rep, "_11_", varE, "/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE, sep=" ")
	    Ped <- Ped[Ped$cat %in% c("potomciNP", "genTest", "telF"),]
	    gen <- merge(genPed, Ped[,c("Indiv", "cat")], by="Indiv", all.y=TRUE)
	    fathers <- genPed[genPed$Indiv %in% gen$Father,c("Indiv", "EBV", "gvNormUnres1")]
	    colnames(fathers) <- c("Father", "FatherEBV", "FatherTGV")
            mothers <- genPed[genPed$Indiv %in% gen$Mother,c("Indiv", "EBV", "gvNormUnres1")]
            colnames(mothers) <- c("Mother", "MotherEBV", "MotherTGV")
 	    gen <- merge(gen, fathers, by="Father")
	    gen <- merge(gen, mothers, by="Mother")
	    gen$PAest <- (gen$MotherEBV + gen$FatherEBV) / 2
	    gen$PAtrue <- (gen$MotherTGV + gen$FatherTGV) / 2
	    
 	    acc <- rbind(acc, c(rep, varE, ref, cor(gen$PAest, gen$gvNormUnres1)))
	}
  }
}

write.table(acc, "PAAcc.csv", quote=FALSE, row.names=FALSE)
