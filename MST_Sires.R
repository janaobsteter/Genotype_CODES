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
PAa <- summarySE(data=m, measurevar = "diffPA", groupvars = c("Strategy", "Scenario") )[,c(1,2,4,5)]
colnames(PAa) <- c("Strategy", "Scenario", "PAdiff", "PAdiff_sd")
MSTa <- summarySE(data=m, measurevar = "MSTmean_PB", groupvars = c("Strategy", "Scenario") )[,c(1,2,4,5)]
colnames(MSTa) <- c("Strategy", "Scenario", "MST", "MST_sd")
MA <- merge(PAa, MSTa, by=c("Strategy", "Scenario"))
MA$Scenario <- revalue(MA$Scenario, c("Class" = "PT", "GenSLO" = "GS-PS", "OtherCowsGen" = "GS-C", "BmGen" = "GS-BD", "Gen" = "GS"))
MA$Strategy <- revalue(MA$Strategy, c("SU55" = "5 sires/year, use 5 years","SU51" = "5 sires/year, use 1 year","SU15" = "1 sire/year, use 5 years"))
MA$Scenario <- factor(MA$Scenario, levels =c("PT", "GS-PS", "GS-C", "GS-BD", "GS"))
MA$Strategy <- factor(MA$Strategy, levels =c("5 sires/year, use 5 years", "5 sires/year, use 1 year", "1 sire/year, use 5 years"))


library(tidyr)
MA$PlotGroup <- paste0(MA$Strategy, MA$Scenario)
MA %>% gather(MST, ~Scenario, ~Strategy)

MST <- ggplot(MA, aes(x=Scenario, y =MST, group=Strategy, fill=Strategy)) + geom_bar(stat="identity", position="dodge") + theme_bw() +
  theme(axis.text=element_text(size=16), legend.position = "none", 
        axis.title.x=element_text(size=16, vjust=-1), 
        axis.title.y=element_text(size=16, vjust=2), 
        legend.text=element_text(size=16), legend.title=element_text(size=16), legend.box = "horizontal",
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        legend.text.align = 0) + ylab("MST selected sires")
PA <- ggplot(MA, aes(x=Scenario, y =PAdiff, group=Strategy, fill=Strategy)) + geom_bar(stat="identity", position="dodge") + theme_bw() +
  theme(axis.text=element_text(size=16), legend.position = "top", 
        axis.title.x=element_text(size=16, vjust=-1), 
        axis.title.y=element_text(size=16, vjust=2), 
        legend.text=element_text(size=16), legend.title=element_text(size=16), legend.box = "horizontal",
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        legend.text.align = 0) + ylab("PA difference")

multiplot(PA, MST)
