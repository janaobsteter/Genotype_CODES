inb <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//InbreedingALL_14082018.csv")
colnames(inb)[c(1,2)] <- c("Generation", "F")



'

ScenarioDF <- data.frame()
for (strategy in unique(inb$strategy)) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      for (interval in c(0,20,40)) {
        inb_temp <- inb[(inb$strategy==strategy) & (inb$scenario==scenario) & 
                          (inb$rep==rep) & (inb$Generation %in% interval:(interval+20)),]
        genInb <- aggregate(inb_temp$F ~ inb_temp$Generation, FUN="mean")
        colnames(genInb) <- c("Gen", "inbreeding")
        genInb$y = log(1 - genInb$inbreeding)
        genInb$t = genInb$Gen - min(genInb$Gen) + 1
        fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
        # lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
        dF = as.data.frame(1 - exp(coef(fit)[2]))[1,1]
        Ne = 1 / (2 * dF)
        ScenarioDF$dF <- as.numeric(ScenarioDF$dF)
        ScenarioDF$Ne <- as.numeric(ScenarioDF$Ne)
        ScenarioDF$Rep <- as.numeric(ScenarioDF$Rep)
        ScenarioDF <- rbind(ScenarioDF, data.frame(Strategy=strategy, 
                                                   Scenario = scenario, Interval = interval, 
                                                   Rep = rep, dF = dF, Ne = Ne))
        colnames(ScenarioDF) <- c("Strategy", "Scenario", "Interval", "Rep", "dF", "Ne")
        
'
ScenarioDF <- data.frame(Strategy=NA, Scenario=NA, Interval = NA, Rep=NA, deltaF=NA, Ne=NA)
for (strategy in c("SU55", "SU51", "SU15")) {
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


#povprečni DF
inb_40 <- inb[inb$Generation %in% 40:60,]
inb_avgF <- aggregate(inb_40$F ~inb_40$strategy + inb_40$scenario + inb_40$Generation, FUN="mean")
colnames(inb_avgF) <- c("Strategy", "Scenario", "Generation", "F")
inb_avgF$Strategy <- factor(inb_avgF$Strategy, levels =c("SU55", "SU51", "SU15"))
inb_avgF$Scenario <- factor(inb_avgF$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
inb_avgF$Scenario <- revalue(inb_avgF$Scenario, c("Class" = "PT", "GenSLO" = "GS-PS", "OtherCowsGen" = "GS-C", "BmGen" = "GS-BD", "Gen" = "GS"))
ggplot(data=inb_avgF, aes(x=Generation, y=F, group=Scenario, colour=Scenario)) + geom_path() + 
  facet_grid(Strategy ~ ., scales = "free_y")

inb_40$dF <- NA
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      for (gen in 41:60) {
        #print(paste(scenario, rep, int, sep=","))
        inb_40$dF[(inb_40$strategy==strategy) & (inb_40$scenario == scenario) & (inb_40$rep == rep) & (inb_40$Generation==gen)] <- 
          inb_40$F[(inb_40$strategy==strategy) & (inb_40$scenario == scenario) & (inb_40$rep == rep) & (inb_40$Generation==gen)] - 
          inb_40$F[(inb_40$strategy==strategy) & (inb_40$scenario == scenario) & (inb_40$rep == rep) & (inb_40$Generation==(gen-1))]
      }
    }
  }
}


inbAvg <- aggregate(inb_40$dF ~ inb_40$strategy + inb_40$scenario + inb_40$rep, FUN="mean")
colnames(inbAvg) <- c("Strategy", "Scenario", "Rep", "dF")
inbAvg$Ne <- as.numeric(1 / (2 * inbAvg$dF))
INBavg <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      #genetic gain
      base <- as.numeric(inbAvg$Ne[inbAvg$Scenario=="Class" & inbAvg$Strategy=="SU55"  & inbAvg$Rep==rep])
      inb <- inbAvg[inbAvg$Scenario==scenario & inbAvg$Strategy==strategy & inbAvg$Rep==rep,]
      inb$per_inb <- inb$Ne / base
      
      INBavg <- rbind( INBavg, inb)
    }
  }
}

INBavg$per_inb <- (INBavg$per_inb)*100 - 100
#INBa <- summarySE(INB, measurevar="Ne", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)] to je za ABSOLUTNE
INBavg_a <- summarySE(INBavg, measurevar="per_inb", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
INBavg_abs <- summarySE(INBavg, measurevar="Ne", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
write.csv(INBavg_abs, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Absolute_pedigreeinbreeding_final.csv", quote=FALSE, row.names=FALSE)
colnames(INBavg_a) <- c("Strategy", "Scenario", "per_inb", "per_inbSD")
colnames(INBavg_abs) <- c("Strategy", "Scenario", "inb", "inbSD")
INBavg_a$per_inb <- round(INBavg_a$per_inb)
INBavg_a$per_inbSD <- round(INBavg_a$per_inbSD)

INBavg_a$Strategy <- factor(INBavg_a$Strategy, levels =c("SU55", "SU51", "SU15"))
INBavg_a$Scenario <- factor(INBavg_a$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
INBavg_a[order(INBavg_a$Strategy, INBavg_a$Scenario),]


ScenarioDF <- ScenarioDF[-1,]
ScenarioDF$Interval <- as.numeric(ScenarioDF$Interval)
ScenarioDF$Interval1 <- paste(ScenarioDF$Interval, ScenarioDF$Interval+19, sep=":")
ScenarioDF$Interval1 <- as.factor(ScenarioDF$Interval1)

#to je izračun procenta razlike od PT SU 5/5
Scenario60 <- ScenarioDF[ScenarioDF$Interval==41,]
Scenario60$Ne <- as.numeric(Scenario60$Ne)
INB <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      #genetic gain
      base <- as.numeric(Scenario60$Ne[Scenario60$Scenario=="Class" & Scenario60$Strategy=="SU55"  & Scenario60$Rep==rep])
      inb <- Scenario60[Scenario60$Scenario==scenario & Scenario60$Strategy==strategy & Scenario60$Rep==rep,]
      inb$per_inb <- inb$Ne / base
      
      INB <- rbind( INB, inb)
    }
  }
}

#tukaj najprej izračunaj povprečje in SD
#############################################################
INB$per_inb <- (INB$per_inb)*100 - 100
#INBa <- summarySE(INB, measurevar="Ne", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)] to je za ABSOLUTNE
INBa <- summarySE(INB, measurevar="per_inb", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
colnames(INBa) <- c("Strategy", "Scenario", "per_inb", "per_inbSD")
INBa$per_inb <- round(INBa$per_inb)
INBa$per_inbSD <- round(INBa$per_inbSD)

INBa$Strategy <- factor(INBa$Strategy, levels =c("SU55", "SU51", "SU15"))
INBa$Scenario <- factor(INBa$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
INBa[order(INBa$Strategy, INBa$Scenario),]

INBa[INBa$Strategy=="SU55", "per_inb"] - INBa[INBa$Strategy=="SU15", "per_inb"]

write.csv(INBa, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NEs_relative_21112018.csv", row.names=FALSE, quote=FALSE)
#############################################################


#izračunaj značilnost
############################################################
#naredi isto se za procente
INB$Strategy <- factor(INB$Strategy, levels =c("SU55", "SU51", "SU15"))
INB$Scenario <- factor(INB$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))

INB <- INB[order(INB$Strategy, INB$Scenario),]


model <- lm(per_inb ~ Strategy + Scenario + Strategy : Scenario, data=INB)
marginal = emmeans(model, ~ Strategy:Scenario)
CLD = cld(marginal, by="Strategy",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

CLD = cld(marginal, by="Scenario",
          alpha   = 0.05, sort=FALSE,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
############################################################


#to je za absolutne cifre
Scenario60$deltaF <- as.numeric(Scenario60$deltaF)
ScenarioA <- aggregate(Scenario60$deltaF ~ Scenario60$Strategy + Scenario60$Scenario, FUN="mean")


colnames(ScenarioA) <- c("Strategy", "Scenario", "dF")
ScenarioA$Method <- "Average"
inbAvg$Method <- "Regression"


HETavg_qtn <- HETavg[HETavg$Marker=="QTN",]
HETavg_qtn <- HETavg_qtn[,c(1,2,4)]
HETavg_qtn$Method <- "QTN"
method <- rbind(ScenarioA, inbAvg)
method <- rbind(method, HETavg_qtn)

method$Scenario <- revalue(method$Scenario, c("Class" = "PT", "GenSLO" = "GS-PS", "OtherCowsGen" = "GS-C", "BmGen" = "GS-BD", "Gen" = "GS"))
method$Scenario <- factor(method$Scenario, levels =c("PT", "GS-PS", "GS-C", "GS-BD", "GS"))
ggplot(data=method, aes(y=dF,x=Scenario,  fill=Method)) + geom_bar(position="dodge", stat="identity") +
  facet_grid(Strategy ~ ., scales = "free_y")
'''
ScenarioDF$Ne <- as.numeric(ScenarioDF$Ne)
ScenarioDFA <- unique(aggregate(ScenarioDF$Ne ~ ScenarioDF$Interval1 + ScenarioDF$Scenario + ScenarioDF$Strategy, FUN="mean"))
ScenarioDFA_sd <- unique(aggregate(ScenarioDF$Ne ~ ScenarioDF$Interval1 + ScenarioDF$Scenario + ScenarioDF$Strategy, FUN="sd"))
colnames(ScenarioDFA) <- c("Interval[gen]", "Scenario","Strategy", "Ne")
colnames(ScenarioDFA_sd) <- c("Interval[gen]", "Scenario","Strategy", "Ne_sd")
ScenarioDF$deltaF <- as.numeric(ScenarioDF$deltaF)
ScenarioDFA_dF <- unique(aggregate(ScenarioDF$deltaF ~ ScenarioDF$Interval1 + ScenarioDF$Scenario + ScenarioDF$Strategy, FUN="mean"))
ScenarioDFA_dF_sd <- unique(aggregate(ScenarioDF$deltaF ~ ScenarioDF$Interval1 + ScenarioDF$Scenario + ScenarioDF$Strategy, FUN="sd"))
colnames(ScenarioDFA_dF) <- c("Interval[gen]", "Scenario","Strategy", "deltaF")
colnames(ScenarioDFA_dF_sd) <- c("Interval[gen]", "Scenario","Strategy", "deltaF_sd")
ScenarioA <- merge(ScenarioDFA, ScenarioDFA_sd, by=c("Interval[gen]", "Scenario","Strategy"))
ScenarioA <- merge(ScenarioA, ScenarioDFA_dF, by=c("Interval[gen]", "Scenario","Strategy"))
ScenarioA <- merge(ScenarioA, ScenarioDFA_dF_sd, by=c("Interval[gen]", "Scenario","Strategy"))

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
levels(INB$Strategy) <- c( "FILLIN", "SU 1/5", "SU 5/1", "SU 5/5")
INB[INB$Strategy=="SU 5/1",]

write.csv(INB, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//NEs_deltaFs_14082018.csv", row.names=FALSE, quote=FALSE)


#to je poolano po replikah
inbRep <- aggregate(inb$F ~  inb$Generation + inb$scenario + inb$strategy, FUN="mean")
inbRepSD <- aggregate(inb$F ~  inb$Generation + inb$scenario + inb$strategy, FUN="sd")
colnames(inbRep) <- c("Generation", "Scenario", "Strategy", "F")
colnames(inbRepSD) <- c("Generation", "Scenario", "Strategy", "sdF")
inbRep <- merge(inbRep, inbRepSD, by=c("Generation", "Scenario", "Strategy"))
inbRep$Scenario <- as.factor(inbRep$Scenario)

#plot pedigree inbreeding
levels(inbRep$Strategy)
levels(inbRep$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5")
inbRep$F <- as.numeric(inbRep$F)
library(ggplot2)
ggplot(data=inbRep, aes(x=Generation, y=F, group = Scenario, colour=Scenario, fill=Scenario)) + 
  geom_line() + 
   xlim (c(40, 60)) + ylim(0.15, 0.30) +
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) +  
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=16), legend.position = "left",
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20),
        strip.text.y = element_text(size = 14))   + 
  facet_grid(Strategy ~ ., scales = "free_y") + theme(legend.position = "right") 
 ''' 
