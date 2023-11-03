######################################################################

OceGen <- data.frame(Generation=NA, NoFathers=NA, Rep=NA, Degree=NA)

for (rep in 0:2) {
   for (degree in c("15", "30", "45", "60", "75", "minCoan", "maxGain")) {
       WorkingDir = paste0("/home/v1jobste/JanaO/Gen", rep, "_", degree, "OCS/SimulatedData/")
       print(WorkingDir)
       ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=TRUE)
       for (gen in 40:60) {
           OceGen <- rbind(OceGen, c(gen, length(unique(ped[ped$Generation==gen, "Father"])), rep, degree)  )
     	}
  }
}
write.table(OceGen, paste0("/home/v1jobste/JanaO/NoFather_OCS/FatherNumber.txt"), quote=FALSE)
