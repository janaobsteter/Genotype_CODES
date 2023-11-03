GenInts <- data.frame()

for (rep in 0:19) {
	for (scenario in c(15, 30, 45,50, 55, 60, 75))  {
		gi <- read.table(paste0("./Gen", rep,"_", scenario, "OCS/GenInts.txt"), header=TRUE)
		gi$Rep <- rep
		gi$scenario <- scenario	
		gi$strategy <- "OCS"
		GenInts <- rbind(GenInts, gi)
	}
}

write.csv(GenInts, "GENINTS_all_20032019_OCS.csv", quote=FALSE)
