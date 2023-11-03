library(dplyr)
library(tidyr)

homedir <- getwd()


for (strategy in c(1, 2, 5, 8, 9, 11)) {
  ACC <- data.frame()
  ACCage <- data.frame()
  for (scenario in c("Class", "Gen")) {
    for (rep in c(0:3)) {
      setwd(paste0(homedir, "/10K/SU55_permEnv/",  scenario, rep, "_", strategy, "/"))
      ped <- read.csv("GenPed_EBV.txt")[,c(1,2,5)]
#      ACC <- data.frame()
#      ACCage <- data.frame()
      print(getwd())
      for (cycle in 40:max(ped$Generation)) {
	print(paste("Cycle: ", cycle))
        pedC <- ped[ped$Generation %in% 1:cycle,]
        sol <- read.table(paste0("renumbered_Solutions_", cycle))[,c(2,3)]
        colnames(sol) <- c("Indiv", "EBV")
        cat <- read.csv(paste0("Categories_gen", cycle, "DF.csv"))
        cats <- gather(cat, category, Indiv)
        acc <- merge(pedC, sol, by="Indiv", all.y=TRUE)
        acc <- merge(acc, cats, all.x=TRUE)
        tmp <- as.data.frame(acc %>% group_by(category) %>% summarize(COR=cor(gvNormUnres1, EBV)))
        tmp$strategy <- strategy
        tmp$scenario <- scenario
        tmp$rep <- rep
        ACC <- rbind(ACC, tmp)
        
        acc$AgeCat <- paste0(acc$category, cycle-acc$Generation)
        tmpAge <- as.data.frame(acc %>% group_by(AgeCat) %>% summarize(COR=cor(gvNormUnres1, EBV)))
        tmpAge$strategy <- strategy
        tmpAge$scenario <- scenario
        tmpAge$rep <- rep
        ACCage <- rbind(ACCage, tmpAge)
      }
      setwd(homedir)
    }
#      write.csv(ACC,paste0( "Accuracy_Cat_permEnv", scenario, ".csv"), quote=FALSE, row.names=FALSE)
#      write.csv(ACCage, paste0("Accuracy_CatAge_permEnv", scenario, ".csv"), quote=FALSE, row.names=FALSE)

  }
  write.csv(ACC,paste0( "Accuracy_Cat", strategy, ".csv"), quote=FALSE, row.names=FALSE)
  write.csv(ACCage, paste0("Accuracy_CatAge", strategy, ".csv"), quote=FALSE, row.names=FALSE)

}


write.csv(ACC, "Accuracy_Cat_permEnv.csv", quote=FALSE, row.names=FALSE)
write.csv(ACCage, "Accuracy_CatAge_permEnv.csv", quote=FALSE, row.names=FALSE)

