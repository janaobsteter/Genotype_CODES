gi <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//GENINTS_all_16072018.csv")
gi$LINE <- paste0(gi$line, gi$sex)
gi <- gi[gi$Gen %in% 40:60,]
gi$genInt <- as.numeric(as.character(gi$genInt))

giA <- aggregate(gi$genInt ~ gi$LINE + gi$scenario + gi$strategy, FUN="mean")
giAS <- aggregate(giA$`gi$genInt` ~ giA$`gi$scenario` + giA$`gi$strategy`, FUN="sum")
colnames(giAS) <- c("Scenario", "Strategy", "genInt")

giAS$genInt[giAS$Strategy=="SU55"] / giAS$genInt[giAS$Strategy=="SU55"]

#gs-c
cbind(giA[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="OtherCowsGen"),], 
      giA[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="Class"),])

1 - (giA$`gi$genInt`[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="OtherCowsGen")] / 
      giA$`gi$genInt`[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="Class")])

#gs-bd
cbind(giA[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="BmGen"),], 
      giA[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="Class"),])

1 - (giA$`gi$genInt`[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="BmGen")] / 
      giA$`gi$genInt`[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="Class")])

#gs
cbind(giA[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="Gen"),], 
      giA[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="Class"),])

1 - (giA$`gi$genInt`[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="Gen")] / 
      giA$`gi$genInt`[(giA$`gi$LINE`=="sireM") & (giA$`gi$scenario`=="Class")])

cbind(giA[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="Gen"),], 
      giA[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="Class"),])

1 - (giA$`gi$genInt`[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="Gen")] / 
      giA$`gi$genInt`[(giA$`gi$LINE`=="sireF") & (giA$`gi$scenario`=="Class")])

#SU5/5 - SU5/1
SU <- cbind(giA[(giA$`gi$strategy`=="10K_Ref_1Year"),],
  giA[(giA$`gi$strategy`=="10K_Ref_20Rep"),])

SU$diff <- 1 - (giA$`gi$genInt`[(giA$`gi$strategy`=="SU51")] / 
       giA$`gi$genInt`[(giA$`gi$strategy`=="SU55")])
SU[order(-SU$diff),]
