######################################################################

OceGen <- data.frame(Generation=NA, NoFathers=NA, Rep=NA, Degree=NA)

for (rep in 0:9) {
   for (degree in c(50,55)) {
       WorkingDir = paste0("/home/v1jobste/JanaO/Gen", rep, "_", degree, "OCS/SimulatedData/")
       print(WorkingDir)
       ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=TRUE)
       for (gen in 40:60) {
           OceGen <- rbind(OceGen, c(gen, length(unique(ped[ped$Generation==gen, "Father"])), rep, degree)  )
     	}
  }
}
write.table(OceGen, paste0("/home/v1jobste/JanaO/NoFather_OCS_20032019_5055.txt"), quote=FALSE)
