######################################################################

OceGen <- data.frame(Generation=NA, NoFathers=NA, Rep=NA, Scenario=NA, Strategy=NA)

for (rep in 0:19) {
   for (strategy in c("SU55", "SU51")) {
     for (scenario in c("Class", "Gen")) {
       WorkingDir = paste0("/home/v1jobste/JanaO/10K_SireUse/", strategy, "/", scenario,rep, "/")
       print(WorkingDir)
       ped <- read.csv(paste0(WorkingDir,'/GenPed_EBV.txt'), header=TRUE)
       for (gen in 40:60) {
           OceGen <- rbind(OceGen, c(gen, length(unique(ped[ped$Generation==gen, "Father"])), rep, scenario, strategy)  )
     	}
  }
}
}
write.table(OceGen, paste0("/home/v1jobste/JanaO/NoFather_Trunc_20032019.txt"), quote=FALSE)
