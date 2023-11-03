
ped <- read.table("./SU55/Class0/SimulatedData/PedigreeAndGeneticValues.txt", header=TRUE)[,c(1,2)]
INB <- data.frame()
for (str in c("SU55", "SU51", "SU15")) {
	INBALL <- data.frame()
	for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
		for (rep in 0:19) {
			print(c(str, scenario, rep))
			inb_temp <- read.table(paste0("./", str, "/Inbreeding_", scenario, rep, ".txt"))
			colnames(inb_temp) <- c("Indiv", "F")
			inb_temp <- merge(inb_temp, ped, by="Indiv")
			inb_tempA <- aggregate(inb_temp$F ~ inb_temp$Generation, FUN="mean")
			inb_tempA$strategy <- str
			inb_tempA$scenario <- scenario
			inb_tempA$rep <- rep
			INB <- rbind(INB, inb_tempA)
			INBALL <- rbind(INBALL, inb_tempA)
		}
	}
	write.csv(INBALL, paste0("./PedigreeInbreeding/Inbreeding_", str, ".csv"), quote=FALSE, row.names=FALSE)
}
write.csv(INB, "./PedigreeInbreeding/InbreedingALL_14082018.csv", quote=FALSE, row.names=FALSE)



	

