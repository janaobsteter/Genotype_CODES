######################################################################


for (rep in 0) {
   for (degree in c("15")) {
       WorkingDir = paste0("/home/v1jobste/JanaO/Gen", rep, "_", degree, "OCS/SimulatedData/")
       print(WorkingDir)
       ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=TRUE)
       ped$Generation <- as.numeric(ped$Generation)
       print(head(ped))
       OceGen <- data.frame(Generation=NA, NoFathers=NA, Rep=NA, Degree=NA)
       for (gen in c(40, 41)) {
	   print(c(gen, length(unique(ped[ped$Generation==as.numeric(gen), "Father"])), rep, degree))

           OceGen <- rbind(OceGen, c(gen, length(unique(ped[ped$Generation==as.numeric(gen), "Father"])), rep, degree))
     	}
	print(OceGen)
    	write.table(OceGen, paste0("/home/v1jobste/JanaO/NoFather_OCS/Father", rep, degree, ".txt"), quote=FALSE)
  }
}

