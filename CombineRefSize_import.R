GenInts <- data.frame()	

for (strategy in c("GenGen", "ClassGen")) {
	for (rep in 0) {
		for (scenario in c("0_0", "100_100", "50_50", "10_10", "0_100"))  {
		  for (trait in c(12,13)) {
			gih <- read.table(paste0("./10K/SU55_import/", strategy, rep, "_", scenario, trait, "/ReferenceSizehome.txt"), header=F)
			gih$Group="home"
			gii <- read.table(paste0("./10K/SU55_import/", strategy, rep, "_", scenario, trait, "/ReferenceSizeimport.txt"), header=F)
			gii$Group="import"
			gi <- rbind(gih, gii)			
			gi$Rep <- rep
			gi$scenario <- scenario	
			gi$strategy <- strategy
			gi$trait <- trait
			GenInts <- rbind(GenInts, gi)
		  }
		}
	}
}

write.csv(GenInts, "ReferenceSize_Import_15122020.csv", quote=FALSE)
