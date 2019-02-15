library(tidyr)
homedir <- getwd()


for (strategy in c("SU55", "SU51", "SU15")) {
  MST <- data.frame()
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      setwd(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      print(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      if (file.info("GenPed_EBV.txt")$size > 0) {
        ped <- read.csv("GenPed_EBV.txt")
        print("Reading GenPed")
        #C <- gather(read.csv("Categories_gen21DF_0.csv"))
        #C <- C[!(is.na(C$value)),]
        #colnames(C) <- c("Cat", "Indiv")
        #ped <- merge(ped, C, by="Indiv", all.x=TRUE)
      }  else {
        ped <- read.table("./SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
        print("Reading Pedigree")
      }
      ped$MST <- NA
      
      fathers <-  unique(ped$Father[ped$Generation %in% 40:60])
      mothers <-  unique(ped$Mother[ped$Generation %in% 40:60])
      pedFather <- ped[ped$Indiv %in% fathers,c("Indiv", "gvNormUnres1")]
      pedMother <- ped[ped$Indiv %in% mothers,c("Indiv", "gvNormUnres1")]
      colnames(pedFather) <- c("Father", "gvF")
      colnames(pedMother) <- c("Mother", "gvM")
      
      ped1 <- ped[ped$Generation %in% 40:60,]
      nrow(ped1)
      ped1 <- merge(ped1, pedFather, by="Father", all.x=TRUE)
      nrow(ped1)
      ped1 <- merge(ped1, pedMother, by="Mother", all.x=TRUE)
      nrow(ped1)

      
      ped1$MST <- ped1$gvNormUnres1 -  ((ped1$gvF +  ped1$gvM) / 2)
      
      mstA <- summarySE(data=ped1, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
      mstA$Scenario <- scenario
      mstA$Rep <- rep
      
      #only fathers
      ped2 <- ped[ped$Indiv %in% fathers,]
      nrow(ped2)
      ped2 <- merge(ped2, pedFather, by="Father", all.x=TRUE)
      nrow(ped2)
      ped2 <- merge(ped2, pedMother, by="Mother", all.x=TRUE)
      nrow(ped2)
      
      ped2$MST <- ped2$gvNormUnres1 -  ((ped2$gvF +  ped2$gvM) / 2)
      mstCat <- summarySE(data=ped2, measurevar = "MST", groupvars = c("Generation"))#[,c(1,3,4)]
      mstCat$Scenario <- scenario
      mstCat$Rep <- rep
      
      MST <- rbind(MST, mstA)
      MSTCat <- rbind(MSTCat, mstCat)
      
      setwd(homedir)
    }
  }
    write.csv(MST,paste0( "MST_", scenario, ".csv"), quote=FALSE, row.names=FALSE)
    write.csv(MSTCat ,paste0( "MSTCat_", scenario, ".csv"), quote=FALSE, row.names=FALSE)
}
write.csv(MST, "MST_OCS.csv", quote=FALSE, row.names=FALSE)
write.csv(MSTCat, "MSTCat.csv", quote=FALSE, row.names=FALSE)
  
  
  