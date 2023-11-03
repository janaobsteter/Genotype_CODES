homedir <- getwd()

FatherUse <- data.frame(Strategy=NA, Degree=NA, Rep=NA, Father = NA, Use = NA)
for (strategy in c("SU55", "SU51")) {
  for (scenario in c("Class", "Gen")) {
    for (rep in c(0:19)) {
      setwd(paste0(homedir, "/10K_SireUse/", strategy, "/", scenario, rep, "/"))
      print(paste0(homedir, "/10K_SireUse/",strategy, "/", scenario, rep, "/"))
      ped <- read.csv("GenPed_EBV.txt")[,1:5]
      fat <- unique(ped$Father[ped$Generation %in% 41:60])
      for (father in fat) {
        use <- max(ped$Generation[ped$Father == father]) - min(ped$Generation[ped$Father == father])
        FatherUse <- rbind(FatherUse, c(strategy, scenario, rep, father, use))
          }
    }      
  }
}
setwd(homedir)
write.csv(FatherUse, "FatherUse_20032019.csv", quote=FALSE, row.names=FALSE)
