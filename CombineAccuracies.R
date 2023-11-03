accALL <- data.frame()
for (str in c("SU55", "SU15", "SU51")) {
	for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
		tmp <- read.csv (paste0("./", str, "/AccuraciesEBVPA_", scenario, ".csv"))
		tmp$Scenario <- scenario
		tmp$Strategy <- str
		accALL <- rbind(accALL, tmp)
	}
}

write.csv(accALL, "./AccuraciesALL_correlation_14082018.csv",  row.names=FALSE, quote=FALSE)

print("Created: AccuraciesALL_correlation_14082018.csv")
