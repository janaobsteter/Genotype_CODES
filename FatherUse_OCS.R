homedir <- getwd()

FatherUse <- data.frame(Degree=NA, Rep=NA, Father = NA, Use = NA)
for (scenario in c(50, 55)) {
  for (rep in c(0:9)) {
    setwd(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    print(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    ped <- read.csv("GenPed_EBV.txt")[,1:5]
    fat <- unique(ped$Father[ped$Generation %in% 41:60])
    for (father in fat) {
      use <- max(ped$Generation[ped$Father == father]) - min(ped$Generation[ped$Father == father])
      FatherUse <- rbind(FatherUse, c(scenario, rep, father, use))
        }
      }
}
setwd(homedir)
write.csv(FatherUse, "FatherUse_OCS_20032019_5055.csv", quote=FALSE, row.names=FALSE)
