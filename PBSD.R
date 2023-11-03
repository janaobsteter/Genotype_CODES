meansd <- data.frame("Scenario" = NA, "MeanSD" = NA, "Rep" = NA)
sd <- data.frame("Scenario" = NA, "SD" = NA, "Rep" = NA)


for (rep in 1:19) {
  pedC <- read.table(paste0("./Class", rep, "/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE)
  pedSLO <-  read.table(paste0("./GenSLO", rep, "/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE)
  pedGSC <- read.table(paste0("./OtherCowsGen", rep, "/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE)
  pedGSBD <-  read.table(paste0("./BmGen", rep, "/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE)
  pedGS <-  read.table(paste0("./Gen", rep, "/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE)
  

  pedc_PB <- pedC[pedC$cat=="pb"  & pedC$Generation %in% 35:60,]
  print(c(rep, "PT"))
  meansd <- rbind(meansd, c("PT", mean(aggregate(pedc_PB$gvNormUnres1 ~ pedc_PB$Generation, FUN="sd")[,2]), rep))
  sd <- rbind(sd, c("PT", sd(pedc_PB$gvNormUnres1), rep))
  
  pedslo_PB <- pedSLO[pedSLO$cat=="pb"  & pedSLO$Generation %in% 35:60,]
  print(c(rep, "GS-PS"))
  meansd <- rbind(meansd, c("GS-PS", mean(aggregate(pedslo_PB$gvNormUnres1 ~ pedslo_PB$Generation, FUN="sd")[,2]), rep))
  sd <- rbind(sd, c("GS-PS", sd(pedslo_PB$gvNormUnres1), rep))
  
  
  pedgsc_PB <- pedGSC[(pedGSC$cat=="pb" | pedGSC$cat=="gpb")  & pedGSC$Generation %in% 35:60,]
  print(c(rep, "GS-C"))
  meansd <- rbind(meansd, c("GS-C", mean(aggregate(pedgsc_PB$gvNormUnres1 ~ pedgsc_PB$Generation, FUN="sd")[,2]), rep))
  sd <- rbind(sd, c("GS-C", sd(pedgsc_PB$gvNormUnres1), rep))
  
  pedgsbd_PB <- pedGSBD[pedGSBD$cat=="pb"  & pedGSBD$Generation %in% 35:60,]
  print(c(rep, "GS-BD"))
  meansd <- rbind(meansd, c("GS-BD", mean(aggregate(pedgsbd_PB$gvNormUnres1 ~ pedgsbd_PB$Generation, FUN="sd")[,2]), rep))
  sd <- rbind(sd, c("GS-BD", sd(pedgsbd_PB$gvNormUnres1), rep))
  
  print(c(rep, "GS"))
  pedgs_PB <- pedGS[pedGS$cat=="gpb"  & pedGS$Generation %in% 35:60,]
  meansd <- rbind(meansd, c("GS", mean(aggregate(pedgs_PB$gvNormUnres1 ~ pedgs_PB$Generation, FUN="sd")[,2]), rep))
  sd <- rbind(sd, c("GS", sd(pedgs_PB$gvNormUnres1), rep))
  

  
}

sd <- sd[-1,]
meansd <-  meansd[-1,]

write.csv(meansd, "MeanSD.csv", row.names=FALSE, quote=FALSE)
write.csv(sd, "SD.csv", row.names=FALSE, quote=FALSE)
