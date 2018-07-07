gi <- read.csv("~/GENINTS_all.csv")
gi$LINE <- paste0(gi$line, gi$sex)
gi <- gi[gi$Gen %in% 40:60,]
gi$genInt <- as.numeric(as.character(gi$genInt))

giA <- aggregate(gi$genInt ~ gi$LINE + gi$scenario + gi$strategy, FUN="mean")
giAS <- aggregate(giA$`gi$genInt` ~ giA$`gi$scenario` + giA$`gi$strategy`, FUN="sum")
colnames(giAS) <- c("Scenario", "Strategy", "genInt")

giAS$genInt[giAS$Strategy=="10K_Ref_1Year"] / giAS$genInt[giAS$Strategy=="10K_Ref_20Rep"]

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

SU$diff <- 1 - (giA$`gi$genInt`[(giA$`gi$strategy`=="10K_Ref_1Year")] / 
       giA$`gi$genInt`[(giA$`gi$strategy`=="10K_Ref_20Rep")])
SU[order(-SU$diff),]
