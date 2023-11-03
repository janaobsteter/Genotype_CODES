accALL <- data.frame()
for (str in c("SU55", "SU15", "SU51")) {
	for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
		for (rep in 0:19) {
			tmp <- read.csv (paste0("./", str, "/", scenario, rep, "/Accuracies_Cat.csv"))
			tmp$Strategy <- str
			tmp$Scenario <- scenario
			tmp$Rep <- rep
			accALL <- rbind(accALL, tmp)
		}
	}
}

write.csv(accALL, "./AccuraciesALL_Cat.csv",  row.names=FALSE, quote=FALSE)
