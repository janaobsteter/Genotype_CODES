HetAll <- data.frame()

for (scenario in c(15, 30, 45, 60, 75)) {
	for (rep in 0:1) {
		het_temp <- read.csv(paste0("./Gen", rep, "_", scenario, "OCS/", "MeanHet_Neutral.csv"))
		HetAll <- rbind(HetAll, het_temp[-1,])
	}
}

write.csv(HetAll, "MeanHet_NeutralALL_11122018_OCS.csv", quote=FALSE, row.names=FALSE)


	

