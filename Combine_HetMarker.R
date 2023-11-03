HetAll <- data.frame()

for (str in c("SU15", "SU51")) {
	for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
		for (rep in 0:19) {
			het_temp <- read.csv(paste0("./", str, "/", scenario, rep, "/", "MeanHet_Marker.csv"))
			HetAll <- rbind(HetAll, het_temp[-1,])
		}
	}
}

write.csv(HetAll, "MeanHet_MarkerALL_14082018.csv", quote=FALSE, row.names=FALSE)


	

