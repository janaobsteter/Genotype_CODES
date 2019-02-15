gi <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/GENINTS_all_14082018.csv")
gi$LINE <- paste0(gi$line, gi$sex)
gi <- gi[gi$Gen %in% 40:60,]
gi$genInt <- as.numeric(as.character(gi$genInt))

giA <- aggregate(gi$genInt ~ gi$LINE + gi$scenario + gi$strategy + gi$Rep, FUN="mean")
colnames(giA) <- c("Line", "Scenario", "Strategy", "Rep", "genInt")
giA$Strategy <- factor(giA$Strategy, levels =c("SU55", "SU51", "SU15"))
giA$Scenario <- factor(giA$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
giA$genInt <- round(giA$genInt, 1)
giA[order(giA$Strategy, giA$Scenario),][giA$Line=="sireM",]

giPer <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (line in c("sireF", "sireM", "damF", "damM")) {
      for (rep in 0:19) {
        #genetic gain
        base <- giA$genInt[giA$Scenario=="Class" & giA$Strategy=="SU55" & giA$Line==line & giA$Rep==rep]
        giLINE <- giA[giA$Scenario==scenario & giA$Strategy==strategy & giA$Line==line & giA$Rep==rep,]
        giLINE$per_gi <- giLINE$genInt / base
        
        giPer <- rbind(giPer, giLINE)
      }
    }
  }
}
giPer$per_gi <- (giPer$per_gi)*100-100
giPer_a <- summarySE(giPer, measurevar="per_gi", groupvars=c("Strategy", "Scenario", "Line"))[,c(1,2,3,5,6)]
giPer_abs <- summarySE(giPer, measurevar="genInt", groupvars=c("Strategy", "Scenario", "Line"))[,c(1,2,3,5,6)]
colnames(giPer_a) <- c("Strategy", "Scenario", "Line", "per_gi", "per_giSD")
colnames(giPer_abs) <- c("Strategy", "Scenario", "Line", "gi", "giSD")
#giPer_a$per_gi <- round(giPer_a$per_gi)
#giPer_a$per_giSD <- round(giPer_a$per_giSD)

giPer_a$Strategy <- factor(giPer_a$Strategy, levels =c("SU55", "SU51", "SU15"))
giPer_a$Scenario <- factor(giPer_a$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
giPer_a[order(giPer_a$Strategy, giPer_a$Scenario),][giPer_a$Line=="sireM",]
giPer_a[order(giPer_a$Strategy, giPer_a$Scenario),][giPer_a$Line=="sireF",]


giPer_abs$Strategy <- factor(giPer_abs$Strategy, levels =c("SU55", "SU51", "SU15"))
giPer_abs$Scenario <- factor(giPer_abs$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
giPer_abs[,4:5] <- round(giPer_abs[,4:5],2)
giPer_abs[order(giPer_abs$Strategy, giPer_a$Scenario),][giPer_abs$Line=="sireM",]
giPer_abs[order(giPer_abs$Strategy, giPer_a$Scenario),][giPer_abs$Line=="sireF",]

#Sires of dams
model <- lm(per_gi ~ Scenario + Strategy + Scenario:Strategy, data=giPer[giPer$Line %in% c("sireF"),])
model <- lm(genInt ~ Scenario + Strategy + Scenario:Strategy, data=giPer[giPer$Line %in% c("sireF"),])#absolute
marginal = emmeans(model, ~ Scenario:Strategy)
CLD = cld(marginal, sort=FALSE, by="Strategy",
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD
CLD = cld(marginal, sort=FALSE, by="Scenario",
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD

#Sires of sires
model <- lm(per_gi ~ Scenario + Strategy + Scenario:Strategy, data=giPer[giPer$Line %in% c("sireM"),])
model <- lm(genInt ~ Scenario + Strategy + Scenario:Strategy, data=giPer[giPer$Line %in% c("sireM"),])
marginal = emmeans(model, ~ Scenario:Strategy)
CLD = cld(marginal, sort=FALSE, by="Strategy",
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD
CLD = cld(marginal, sort=FALSE, by="Scenario",
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD


giAS <- aggregate(giA$genInt ~ giA$Scenario + giA$Strategy, FUN="sum")
colnames(giAS) <- c("Scenario", "Strategy", "genInt")

giAS$genInt[giAS$Strategy=="SU55"] / giAS$genInt[giAS$Strategy=="SU55"]



#gs-c
cbind(giA[(giA$Line=="sireF") & (giA$Scenario=="OtherCowsGen"),], 
      giA[(giA$Line=="sireF") & (giA$Scenario=="Class"),])

1 - (giA$genInt[(giA$Line=="sireF") & (giA$Scenario=="OtherCowsGen")] / 
       giA$genInt[(giA$Line=="sireF") & (giA$Scenario=="Class")])

#gs-bd
cbind(giA[(giA$Line=="sireM") & (giA$Scenario=="BmGen"),], 
      giA[(giA$Line=="sireM") & (giA$Scenario=="Class"),])

1 - (giA$genInt[(giA$Line=="sireM") & (giA$Scenario=="BmGen")] / 
       giA$genInt[(giA$Line=="sireM") & (giA$Scenario=="Class")])

#gs
cbind(giA[(giA$Line=="sireM") & (giA$Scenario=="Gen"),], 
      giA[(giA$Line=="sireM") & (giA$Scenario=="Class"),])

1 - (giA$genInt[(giA$Line=="sireM") & (giA$Scenario=="Gen")] / 
       giA$genInt[(giA$Line=="sireM") & (giA$Scenario=="Class")])

cbind(giA[(giA$Line=="sireF") & (giA$Scenario=="Gen"),], 
      giA[(giA$Line=="sireF") & (giA$Scenario=="Class"),])

1 - (giA$genInt[(giA$Line=="sireF") & (giA$Scenario=="Gen")] / 
       giA$genInt[(giA$Line=="sireF") & (giA$Scenario=="Class")])

#SU5/5 - SU5/1
SU <- cbind(giA[(giA$Strategy=="SU51"),],
            giA[(giA$Strategy=="SU55"),])

SU$diff <- 1 - (giA$genInt[(giA$Strategy=="SU51")] / 
                  giA$genInt[(giA$Strategy=="SU55")])
SU[order(-SU$diff),]

SUt <- cbind(giAS[giAS$Strategy=="SU55",] , giAS[giAS$Strategy=="SU51",])
SUt$diff <- 1 - (giAS$genInt[(giAS$Strategy=="SU51")] / 
                   giAS$genInt[(giAS$Strategy=="SU55")])
SUt[order(-SUt$diff),]

