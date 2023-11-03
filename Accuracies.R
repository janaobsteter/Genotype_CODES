library(dplyr)
library(tidyr)

homedir <- getwd()


for (strategy in c("SU55", "SU51", "SU15")) {
  ACC <- data.frame()
  ACCage <- data.frame()
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      setwd(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      ped <- read.csv("GenPed_EBV.txt")[,c(1,2,5)]
#      ACC <- data.frame()
#      ACCage <- data.frame()
      print(getwd())
      for (cycle in 40:60) {
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
        
        acc$AgeCat <- paste0(acc$category, 21-acc$Generation)
        tmpAge <- as.data.frame(acc %>% group_by(AgeCat) %>% summarize(COR=cor(gvNormUnres1, EBV)))
        tmpAge$strategy <- strategy
        tmpAge$scenario <- scenario
        tmpAge$rep <- rep
        ACCage <- rbind(ACCage, tmpAge)
      }
      setwd(homedir)
    }
  }
  write.csv(ACC,paste0( "Accuracy_Cat", strategy, ".csv"), quote=FALSE, row.names=FALSE)
  write.csv(ACCage, paste0("Accuracy_CatAge", strategy, ".csv"), quote=FALSE, row.names=FALSE)

}


write.csv(ACC, "Accuracy_Cat.csv", quote=FALSE, row.names=FALSE)
write.csv(ACCage, "Accuracy_CatAge.csv", quote=FALSE, row.names=FALSE)

