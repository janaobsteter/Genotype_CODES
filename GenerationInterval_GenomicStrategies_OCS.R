gi <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/GENINTS_all_14082018.csv")
gi_OCS <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/GENINTS_all_17012019_OCS.csv")
gi_OCS$scenario <- as.character(gi_OCS$scenario)
#zaenkrat izključi OCS30, rep 7
gi_OCS$genInt[gi_OCS$Rep==7 & gi_OCS$strategy=="OCS" & gi_OCS$scenario==30] <- NA


gi <- rbind(gi, gi_OCS)
gi$LINE <- paste0(gi$line, gi$sex)
gi <- gi[gi$Gen %in% 40:60,]
gi$genInt <- as.numeric(as.character(gi$genInt))

gi$scenario <- revalue(gi$scenario, c("Class" = "PT", "GenSLO" = "GS-PS", "OtherCowsGen" = "GS-C", "BmGen" = "GS-BD", "Gen" = "GS",
                                            "15"="15", "30"="30", "45"="45", "60"="60", "75"="75"))

gi$Group <- paste0(gi$strategy, gi$scenario)

#here aggregate the genInt by Generation
giA <- aggregate(gi$genInt ~ gi$LINE + gi$scenario + gi$strategy + gi$Rep, FUN="mean")
colnames(giA) <- c("Line", "Scenario", "Strategy", "Rep", "genInt")
giA$Strategy <- factor(giA$Strategy, levels =c("SU55", "SU51", "SU15",  "OCS"))
#giA$Scenario <- factor(giA$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
giA$genInt <- round(giA$genInt, 1)
giA[order(giA$Strategy, giA$Scenario),][giA$Line=="sireM",]


giA$Group <- paste0(giA$Strategy, giA$Scenario)
giA <- giA[giA$Group %in% c("SU55PT","SU55GS", "SU51GS", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75"),]
abs <- aggregate(giA$genInt ~ giA$Strategy + giA$Scenario + giA$Line, FUN="mean")
abs[abs$`giA$Line` == "sireM",]
abs[abs$`giA$Line` == "sireF",]

giPer <- data.frame()
for (group in c("SU55PT","SU55GS", "SU51GS", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75")) {
  for (line in c("sireF", "sireM", "damF", "damM")) {
    for (rep in 0:19) {
      #genetic gain
      base <- giA$genInt[giA$Scenario=="PT" & giA$Strategy=="SU55" & giA$Line==line & giA$Rep==rep]
      giLINE <- giA[giA$Group==group & giA$Line==line & giA$Rep==rep,]
      giLINE$per_gi <- giLINE$genInt / base
      
      giPer <- rbind(giPer, giLINE)
    }
  }
}

giPer$per_gi <- (giPer$per_gi)*100-100
giPer[giPer$Group=="OCS30",]

giPer_a <- summarySE(giPer, measurevar="per_gi", groupvars=c("Strategy", "Scenario", "Line"))[,c(1,2,3,5,6)]
colnames(giPer_a) <- c("Strategy", "Scenario", "Line", "per_gi", "per_giSD")
#giPer_a$per_gi <- round(giPer_a$per_gi)
#giPer_a$per_giSD <- round(giPer_a$per_giSD)

giPer_a$Strategy <- factor(giPer_a$Strategy, levels =c("SU55", "SU51", "OCS"))
giPer_a$Scenario <- factor(giPer_a$Scenario, levels =c("PT", "GS", 15, 30, 45, 60, 75))
giPer_a[order(giPer_a$Strategy, giPer_a$Scenario),][giPer_a$Line=="sireM",]
giPer_a[order(giPer_a$Strategy, giPer_a$Scenario),][giPer_a$Line=="sireF",]

#končna tabela
giPer_a[giPer_a$Line=="sireM",]
giPer_a[giPer_a$Line=="sireF",]
#giPer$Group <- factor(giPer$Group, level=c("SU55PT","SU55GS","SU51GS","OCS15","OCS30","OCS45","OCS60","OCS75"))

#preveri, kaj se dogaja z OCS30
ocs30 <- gi[gi$strategy=="OCS" & gi$scenario==30,]
summary(ocs30$genInt)
ocs30[ocs30$LINE == "sireF",]
sd(ocs30$genInt[ocs30$LINE == "sireF"])
sd(gi$genInt[gi$strategy=="OCS" & gi$scenario==45 & gi$LINE == "sireF"])

giPer$Group <- factor(giPer$Group, level=c("SU55PT","SU55GS","SU51GS","OCS15","OCS30","OCS45","OCS60","OCS75"))

model <- lm(per_gi ~ Group + Line + Group:Line, data=giPer)
marginal = emmeans(model, ~ Group:Line)
CLD = cld(marginal, sort=FALSE, by="Line",
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD



#nUMBER OF FATHERS
noF <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NoFathers_OCS_11122018.csv")
noF <- noF[noF$Generation %in% 41:60,]
colnames(noF) <- c("Scenario", "Rep", "Generation", "NoFathers")

summarySE(noF, measurevar="NoFathers", groupvars = c("Scenario"))
