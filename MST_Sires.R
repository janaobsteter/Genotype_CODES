mstAa <- summarySE(data=ped1, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
colnames(mstAa) <- c("Generation", "MSTmean", "MSTsd")
mstAb <- summarySE(data=ped1, measurevar = "PA", groupvars = "Generation")[,c(1,3,4)]
colnames(mstAb) <-        c("Generation", "PAmean", "PAsd")
mstA <- merge(mstAa, mstAb, by="Generation")
mstA$Scenario <- scenario
mstA$Rep <- rep



pedPB <- ped1[ped1$cat %in% c("pb", "gpb"),]
PBmst <- summarySE(data=pedPB, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
colnames(PBmst) <- c("Generation", "MSTmean_PB", "MSTsd_PB")
PBpa <- summarySE(data=pedPB, measurevar = "PA", groupvars = "Generation")[,c(1,3,4)]
colnames(PBpa) <-        c("Generation", "PAmean_PB", "PAsd_PB")
mstA <- merge(mstA, PBmst, by="Generation")
mstA <- merge(mstA, PBpa, by="Generation")


###
#preberi MST.csv

m <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MST.csv")
table(m$Strategy, m$Generation)
m$diffPA <- m$PAmean_PB - m$PAmean
PAa <- summarySE(data=m, measurevar = "diffPA", groupvars = c("Strategy", "Scenario") )[,c()]
colnames(PAa) <- c("Strategy", "Scenario", "PAdiff", "PAdiff_sd")
MSTa <- summarySE(data=m, measurevar = "MSTmean_PB", groupvars = c("Strategy", "Scenario") )[,c()]
colnames(PAa) <- c("Strategy", "Scenario", "MST", "MST_sd")
merge(PAa, MSTa, by=c("Strategy", "Scenario"))
