GenVar <- data.frame()

for (strategy in c("GenGen", "ClassGen")) {
	for (rep in 0) {
		for (scenario in c("0_0", "100_100", "50_50", "10_10", "0_100"))  {
		  for (trait in c(12,13)) {
			g <- read.table(paste0("./10K/SU55_import/", strategy, rep, "_", scenario, trait, "/GenicVariance_import.csv"), header=TRUE)
			g$Rep <- rep
			g$scenario <- scenario	
			g$strategy <- strategy
			g$trait <- trait
			GenVar <- rbind(GenVar, g)
		  }
		}
	}
}

write.csv(GenVar, "GenicVariance_import_28102020.csv", quote=FALSE)
