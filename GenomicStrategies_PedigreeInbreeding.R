inb <- read.csv("~/InbreedingALL.csv")
colnames(inb)[c(1,2)] <- c("Generation", "F")

ScenarioDF <- data.frame(Strategy=NA, Scenario=NA, Interval = NA, Rep=NA, deltaF=NA, Ne=NA)
for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      for (interval in c(1,21,41)) {
        #print(paste(scenario, rep, int, sep=","))
        genInb <- inb[(inb$strategy==strategy) & (inb$scenario == scenario) & (inb$rep == rep) & (inb$Generation %in% interval:(interval+19)),]
        genInb$y = log(1 - genInb$F)
        genInb$t = genInb$Generation - min(genInb$Generation) + 1
        fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
        # lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
        dF = 1 - exp(coef(fit)[2])
        Ne = as.numeric(1 / (2 * dF))
        ScenarioDF <- rbind(ScenarioDF, c(strategy, scenario, interval, rep, dF, Ne))
      }
    }
  }
}

ScenarioDF <- ScenarioDF[-1,]
ScenarioDF$Interval <- as.numeric(ScenarioDF$Interval)
ScenarioDF$Interval1 <- paste(ScenarioDF$Interval, ScenarioDF$Interval+19, sep=":")
ScenarioDF$Interval1 <- as.factor(ScenarioDF$Interval1)



ScenarioDF$Ne <- as.numeric(ScenarioDF$Ne)
ScenarioDFA <- unique(aggregate(ScenarioDF$Ne ~ ScenarioDF$Interval1 + ScenarioDF$Scenario + ScenarioDF$Strategy, FUN="mean"))
colnames(ScenarioDFA) <- c("Interval[gen]", "Scenario","Strategy", "Ne")
ScenarioDF$deltaF <- as.numeric(ScenarioDF$deltaF)
ScenarioDFA_dF <- unique(aggregate(ScenarioDF$deltaF ~ ScenarioDF$Interval1 + ScenarioDF$Scenario + ScenarioDF$Strategy, FUN="mean"))
colnames(ScenarioDFA_dF) <- c("Interval[gen]", "Scenario","Strategy", "deltaF")
ScenarioA <- merge(ScenarioDFA, ScenarioDFA_dF, by=c("Interval[gen]", "Scenario","Strategy"))
INB <- ScenarioA[ScenarioA$`Interval[gen]` =="41:60",]
INB1 <- ScenarioA[ScenarioA$`Interval[gen]` !="41:60",]
INB1$Strategy <- "FILLIN"
INB1$Scenario <- "FILLIN"
INB1 <- unique(INB1)
INB <- rbind(INB1, INB)
INB$Ne <- round(INB$Ne, 1)
INB$Scenario <- as.factor(INB$Scenario)
levels(INB$Scenario)
levels(INB$Scenario) <- c("GS-BD", "PT", "FILLIN", "GS", "GS-PS", "GS-C")
INB$Strategy <- as.factor(INB$Strategy)
levels(INB$Strategy)
levels(INB$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5", "FILLIN")

write.csv(INB, "~/Documents/PhD/Simulaton/NEs_deltaFs.csv", row.names=FALSE, quote=FALSE)

inbRep <- aggregate(inb$F ~  inb$Generation + inb$scenario + inb$strategy, FUN="mean")
colnames(inbRep) <- c("Generation", "Scenario", "Strategy", "F")
inbRep$Scenario <- as.factor(inbRep$Scenario)

#plot pedigree inbreeding
levels(inbRep$Strategy)
levels(inbRep$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5")
library(ggplot2)
ggplot(data=inbRep, aes(x=Generation, y=F, group = Scenario, colour=Scenario, fill=Scenario)) + 
   geom_path() + 
   xlim (c(40, 60)) + ylim(c(0.16, 0.35)) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) +  
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=16), legend.position = "left",
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20),
        strip.text.y = element_text(size = 14))   + 
  facet_grid(Strategy ~ ., scales = "free_y") + theme(legend.position = "right") 
  
