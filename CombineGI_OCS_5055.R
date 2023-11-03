GenInts <- data.frame()

for (rep in 0:9) {
	for (scenario in c(50, 55))  {
		gi <- read.table(paste0("./Gen", rep,"_", scenario, "OCS/GenInts.txt"), header=TRUE)
		gi$Rep <- rep
		gi$scenario <- scenario	
		gi$strategy <- "OCS"
		GenInts <- rbind(GenInts, gi)
	}
}

write.csv(GenInts, "GENINTS_5055_20032019_OCS.csv", quote=FALSE)
