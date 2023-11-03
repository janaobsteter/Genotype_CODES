library(tidyr)

homedir <- "/home/v1jobste/JanaO/10K_SireUse/"

for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      setwd(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      ped <- read.csv("GenPed_EBV.txt")
      print(getwd())
      for (cycle in 40:60) {
        print(paste("Cycle: ", cycle))
        cat <- read.csv(paste0("Categories_gen", cycle, "DF.csv"))
        cats <- gather(cat, category, Indiv)
        colnames(cats) <- c(paste0("Category", cycle), "Indiv")
        ped <- merge(ped, cats, by="Indiv", all.x=TRUE)
      }
      write.csv(ped, paste0(homedir, "PedCat/", strategy, "_", scenario, rep, "_GenPed_EBV_CATS.txt"), quote=FALSE, row.names=FALSE)
      setwd(homedir)

    }
  }
  
}


