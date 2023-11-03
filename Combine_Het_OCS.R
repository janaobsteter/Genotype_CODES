HetAll <- data.frame()

for (scenario in c(15, 30, 45, 60, 75)) {
  for (rep in 0:19) {
    hetN_temp <- read.csv(paste0("./Gen", rep, "_", scenario, "OCS/", "MeanHet_Neutral.csv"))
    hetN_temp$Marker <- "Neutral" 
    hetN_temp$degree <- scenario
    colnames(hetN_temp)[5] <- "Het"
    hetM_temp <- read.csv(paste0("./Gen", rep, "_", scenario, "OCS/", "MeanHet_Marker.csv"))
    hetM_temp$Marker <- "Marker"
    hetM_temp$degree <-	scenario
    colnames(hetM_temp)[5] <- "Het"
    hetQ_temp <- read.csv(paste0("./Gen", rep, "_", scenario, "OCS/", "MeanHet_QTN.csv"))
    hetQ_temp$Marker <- "QTN"
    hetQ_temp$degree <-	scenario
    colnames(hetQ_temp)[5] <- "Het"
    HetAll <- rbind(HetAll, hetN_temp[-1,])
    HetAll <- rbind(HetAll, hetM_temp[-1,])
    HetAll <- rbind(HetAll, hetQ_temp[-1,])
  }
}

write.csv(HetAll, "MeanHet_ALL_11022019_OCS.csv", quote=FALSE, row.names=FALSE)
