GenInts <- data.frame()

for (strategy in c("SU55", "SU51", "SU15")) {
	for (rep in 0:19) {
		for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))  {
			gi <- read.table(paste0("./", strategy, "/", scenario, rep, "/GenInts.txt"), header=TRUE)
			gi$Rep <- rep
			gi$scenario <- scenario	
			gi$strategy <- strategy
			GenInts <- rbind(GenInts, gi)
		}
	}
}

write.csv(GenInts, "GENINTS_all_14082018.csv", quote=FALSE)
