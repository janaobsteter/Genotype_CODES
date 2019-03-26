#Accuracy of sires in OCS
accOCS <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/ACC_OCS.csv")
meanRep <- summarySE(accOCS, measurevar = "Cor", groupvars = c("Degree", "Rep"))
meanDegree <- summarySE(meanRep, measurevar = "Cor", groupvars = "Degree")
meanDegree


#test significance of sires accuracy
acc <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/ACCAge.csv")[,-1]
acc$scenario <- revalue(acc$scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT"))
acc$Group <- paste0(acc$strategy, acc$scenario)
acc <-  acc[acc$Group %in% c("SU55PT", "SU55GT", "SU51GT") & acc$AgeCat %in% c("cak5", "genTest"),]
head(acc)
head(accOCS)

accOCSa <- summarySE(data = accOCS, measurevar = "Cor", groupvars =  c("Degree", "Rep"))
accOCSa$Strategy <- "OCS"

acc <- acc[c("COR", "strategy", "scenario", "rep")]
colnames(acc) <- c("Cor", "Strategy", "Scenario", "Rep")
accOCSa <- accOCSa[c("Degree", "Rep", "Cor", "Strategy")]
colnames(accOCSa)[1] <- "Scenario"

accOCSa$Scenario <- as.factor(accOCSa$Scenario)
acc$Scenario <- as.factor(acc$Scenario)
accCompare <- rbind(acc, accOCSa)
accCompare$Group <- paste0(accCompare$Strategy, accCompare$Scenario)

accCompare$Group <- factor(accCompare$Group, level=c("SU55PT","SU55GT","SU51GT","OCS15","OCS30","OCS45","OCS60","OCS75"))

model <- lm(Cor ~ Group, data=accCompare)
marginal = emmeans(model, ~ Group)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD

######################################################################################################
######################################################################################################
######################################################################################################
#TGVsAll <- read.csv("~/TGVSALL_11062018.csv")
TGVsAll <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/TGVSALL_14082018.csv")
TGVsOCS <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//TGVsAll_OCS_11022019.csv", sep=" ")
#TGVsOCS <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//TGVsAll_OCS_11122018.csv", sep=" ")

TGVsOCS$degree <- as.factor(TGVsOCS$degree)
colnames(TGVsOCS)[colnames(TGVsOCS)=="degree"] <- "scenario"
TGVsOCS$scenario <- as.factor(TGVsOCS$scenario)
TGVsAll <- rbind(TGVsAll, TGVsOCS)
#TGVsAll1 <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_ReferenceSize//Results/TGVSALL_22082018.csv")
TGVsAll$strategy <-TGVsAll$Strategy
TGVsAll$PlotGroup <- paste0(TGVsAll$Strategy, TGVsAll$scenario)


TGVsAll$group <- paste0(TGVsAll$scenario, TGVsAll$Rep)
library(plyr)
TGVsAll$scenario <- revalue(TGVsAll$scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT",
                                                "15"="15", "30"="30", "45"="45", "60"="60", "75"="75"))


#add genetic and genic variance (standardised)
TGVsAll1 <- data.frame()
for (strategy in c("SU55", "SU51", "SU15", "OCS")) {
  for (scenario in c("PT", "GT-PT", "GT-C", "GT-BD", "GT", 15, 30, 45, 60, 75)) {
    for (rep in 0:19) {
      TGVs <- TGVsAll[TGVsAll$Strategy==strategy & TGVsAll$scenario==scenario & TGVsAll$Rep==rep,]
      TGVs$GeneticVarSt <- TGVs$var / TGVs$var[1]
      TGVs$GenicVarSt <- TGVs$AdditGenicVar1 / TGVs$AdditGenicVar1[1]
      
      TGVsAll1 <- rbind(TGVsAll1, TGVs)  
    }
  }
}

TGVsAll <- TGVsAll1

#TGVsAll <- TGVsAll[TGVsAll$Rep %in% 0:2,]
#TGVsAll$strategy <- NA
#TGVsAll$strategy[TGVsAll$Strategy == "10K_Ref_20Rep"] <- "SU55"
#TGVsAll$strategy[TGVsAll$Strategy == "10K_Ref_1Pb"] <- "SU15"
#TGVsAll$strategy[TGVsAll$Strategy == "10K_Ref_1Year"] <- "SU51"
#TGVsAll <- TGVsAll[TGVsAll$Rep %in% 0:19,] #so that you have 20 and not 21 replicates
#TGVsAll <- TGVsAll[TGVsAll$Rep %in% 10:20,]
#naredi povprečja replik
Averages1 <- aggregate(TGVsAll$SDGenicSt ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages2 <- aggregate(TGVsAll$zMeanGenic ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages3 <- aggregate(TGVsAll$SDSt ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages4 <- aggregate(TGVsAll$zMean ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages5 <- aggregate(TGVsAll$SDGenic ~  TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages6 <- aggregate(TGVsAll$sd ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages7 <- aggregate(TGVsAll$zMeanGenic_Genetic ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages8 <- aggregate(TGVsAll$Genic_Genetic_SDSt ~ TGVsAll$Strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
colnames(Averages1) <- c("strategy", "scenario", "Generation", "SDGenicSt") #stand
colnames(Averages2) <- c("strategy","scenario", "Generation", "MeanGenic") #stand
colnames(Averages3) <- c("strategy","scenario", "Generation", "SDGenetic") #stand
colnames(Averages4) <- c("strategy","scenario", "Generation", "MeanGenetic") #standard
colnames(Averages5) <- c("strategy","scenario", "Generation", "SDGenic_noSt")
colnames(Averages6) <- c("strategy","scenario", "Generation", "SDGenetic_noSt")
colnames(Averages7) <- c("strategy","scenario", "Generation", "MeanGenicGenetic") #standardizirana
colnames(Averages8) <- c("strategy","scenario", "Generation", "SDGenic_Genetic") #standardizirana
Averages <- merge(Averages1, Averages2, by=c( "strategy","scenario", "Generation"))
Averages <- merge(Averages, Averages3, by=c( "strategy","scenario", "Generation"))
Averages <- merge(Averages, Averages4, by=c( "strategy","scenario", "Generation"))
Averages <- merge(Averages, Averages5, by=c( "strategy","scenario", "Generation"))
Averages <- merge(Averages, Averages6, by=c( "strategy","scenario", "Generation"))
Averages <- merge(Averages, Averages7, by=c( "strategy","scenario", "Generation"))
Averages <- merge(Averages, Averages8, by=c( "strategy","scenario", "Generation"))


######################################################
######################################################
######################################################
######################################################
#izračunaj procente glede na SU55 PT
OCS60 <- TGVsAll[TGVsAll$Generation==60,]
summarySE(data=OCS60 , measurevar = "zMean", groupvars = c("Strategy", "scenario"))
tgv60 <- data.frame()

for (scenario in c("PT","GT", 15, 30, 45, 60, 75)) {
  for (rep in 0:19) {
    #genetic gain
    base <- OCS60$zMean[OCS60$scenario=="PT" & OCS60$Strategy=="SU55" & OCS60$Rep==rep]
    tgv <- OCS60[OCS60$scenario==scenario & OCS60$Rep==rep,]
    tgv$per_zMean <- tgv$zMean / base

    
    #genic sd
    base <- OCS60$SDGenicSt[OCS60$scenario=="PT" & OCS60$Strategy=="SU55" & OCS60$Rep==rep]
    tgv$per_GenicSD <- (tgv$SDGenicSt / base)
    
    #genetic sd
    base <- OCS60$SDSt[OCS60$scenario=="PT" & OCS60$Strategy=="SU55" & OCS60$Rep==rep]
    tgv$per_GeneticSD <- (tgv$SDSt / base)      
    
    
    tgv60 <- rbind(tgv60, tgv)
    
  }
}

tgv60$PlotGroup <- paste0(tgv60$strategy, tgv60$scenario)
#OCS60 <- tgv60[tgv60$Strategy %in% c("SU55", "OCS"),]
#OCS60 <- tgv60[tgv60$scenario %in% c("PT","GT", 15, 30, 45, 60, 75),]
OCS60 <- tgv60[tgv60$PlotGroup %in% c("SU55PT","SU55GT", "SU51GT", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75"),]
table(OCS60$Strategy, OCS60$scenario)
OCS60$per_zMean <- OCS60$per_zMean*100 - 100 #tukaj ma osnoven scenarij 100: zato odštej 100
OCS60$per_GenicSD <- (OCS60$per_GenicSD)*100 - 100#tukaj ma osnoven scenarij 0
OCS60$per_GeneticSD <- (OCS60$per_GeneticSD)*100 - 100

#genetic mean
MEAN60_OCS <- summarySE(OCS60, measurevar="per_zMean", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
MEAN60_OCS_abs <- summarySE(OCS60, measurevar="zMean", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(MEAN60_OCS) <- c("Strategy", "Scenario", "per_zMean", "per_zMeanSD")
colnames(MEAN60_OCS_abs) <- c("Strategy", "Scenario", "zMean", "zMeanSD")

#genetic sd
SD60_OCS <- summarySE(OCS60, measurevar="per_GeneticSD", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(SD60_OCS) <- c("Strategy", "Scenario", "per_GeneticSD", "per_GeneticSDSD")
SD60_OCS_abs <- summarySE(OCS60, measurevar="SDSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(SD60_OCS_abs) <- c("Strategy", "Scenario", "GeneticSD", "GeneticSDSD")
#SD60a$per_GeneticSD <- round(SD60a$per_GeneticSD, 3)
#SD60a$per_GeneticSDSD <- round(SD60a$per_GeneticSDSD, 3)

MEAN60_OCS <- merge(MEAN60_OCS, SD60_OCS, by=c("Strategy", "Scenario"))
MEAN60_OCS_abs <- merge(MEAN60_OCS_abs, SD60_OCS_abs, by=c("Strategy", "Scenario"))

#genic sd
SD60_OCSa <- summarySE(OCS60, measurevar="per_GenicSD", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
SD60_OCSa_abs <- summarySE(OCS60, measurevar="SDGenicSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
#SD60b <- summarySE(TGV60, measurevar="SDGenicSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(SD60_OCSa) <- c("Strategy", "Scenario", "per_GenicSD", "per_GenicSDSD")
colnames(SD60_OCSa_abs) <- c("Strategy", "Scenario", "GenicSD", "GenicSDSD")
#SD60b$per_GenicSD <- round(SD60b$per_GenicSD, 3)
#SD60b$per_GenicSDSD <- round(SD60b$per_GenicSDSD, 3)

MEAN60_OCS <- merge(MEAN60_OCS, SD60_OCSa, by=c("Strategy", "Scenario"))
MEAN60_OCS_abs <- merge(MEAN60_OCS_abs, SD60_OCSa_abs, by=c("Strategy", "Scenario"))



#EFFICIENCY
#to je povrepčje regresij
library(nlme)
regRep <- data.frame(Rep=NA, Intercept=NA, Slope=NA, Scenario=NA, Strategy=NA)

for (plotgroup in c("SU55PT","SU55GT", "SU51GT", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75")) {
  #df <- TGVsAll[(TGVsAll$scenario==sc & TGVsAll$Strategy %in% c("SU55", "OCS")) & TGVsAll$Rep %in% 0:19,]
  df <- TGVsAll[TGVsAll$PlotGroup==plotgroup & TGVsAll$Rep %in% 0:19,]
  if (nrow(df) > 0) {
    fm1 <- lmList(zMeanGenic ~ SDGenicStNeg | Rep, data=df, pool = TRUE)
    tmp <- data.frame(Rep=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
    colnames(tmp) <- c("Rep", "Intercept", "Slope")
    tmp$Scenario <- unique(df$scenario)
    tmp$Strategy <- unique(df$Strategy)
    regRep <- rbind(regRep, tmp)
  }
}

regRep <- regRep[-1,]

#povprečje regresijskih koefificentov - efficiency
avGTlope <- summarySE(regRep, measurevar="Slope", groupvars=c("Strategy", "Scenario"), na.rm=TRUE)[,c(1,2,4,5)]
colnames(avGTlope) <- c("Strategy", "Scenario", "Slope", "SlopeSD")
avgInt <- summarySE(regRep, measurevar="Intercept", groupvars=c("Strategy", "Scenario"), na.rm=TRUE)[,c(1,2,4,5)]
colnames(avgInt) <- c("Strategy", "Scenario", "Intercept", "InterceptSD")
avgReg <- merge(avGTlope, avgInt, by=c("Strategy", "Scenario"))

avgReg$Eff <- round(avgReg$Slope, 0)
avgReg$EffSd <- round(avgReg$SlopeSD, 0)

#write.csv(avgReg, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//Efficiency_genomicstrategies_OCS_17012018.csv", quote=FALSE)

#efficiency of SU55
EFF <- data.frame()

for (scenario in c("PT", "GT", 15, 30, 45, 60, 75)) {
  for (rep in 0:19) {
    effBase <- regRep$Slope[regRep$Scenario=="PT"  & regRep$Strategy=="SU55" & regRep$Rep==rep]
    eff <- regRep[regRep$Scenario==scenario & regRep$Rep==rep,]
    eff$per_Eff <- (eff$Slope / effBase)
    EFF <- rbind(EFF, eff)
  }
}


EFF <- EFF[!(is.na(EFF$Slope)),]
EFF$per_Eff <-  EFF$per_Eff * 100 -100
EFF$Name <- paste0(EFF$Strategy, EFF$Scenario)
EFF$Name <- factor(EFF$Name, level=c("SU55PT","SU55GT","SU51GT","OCS15","OCS30","OCS45","OCS60","OCS75"))

#efficiency
eff <- summarySE(EFF, measurevar="per_Eff", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
eff_abs <- summarySE(EFF, measurevar="Slope", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
colnames(eff) <- c("Strategy", "Scenario", "per_Eff", "per_EffSD")
colnames(eff_abs) <- c("Strategy", "Scenario", "Eff", "EffSD")
#eff$per_Eff <- round(eff$per_Eff)
#eff$per_EffSD <- round(eff$per_EffSD)

MEAN60_OCS <- merge(MEAN60_OCS, eff, by=c("Strategy", "Scenario"))
MEAN60_OCS_abs <- merge(MEAN60_OCS_abs, eff_abs, by=c("Strategy", "Scenario"))




library(emmeans)
OCS60$scenario <- as.factor(OCS60$scenario)
OCS60$Strategy <- as.factor(OCS60$Strategy)
OCS60$scenario <- factor(OCS60$scenario, levels =c("PT", "GT", 15,30,45,60,75))
MEAN60_OCS$Scenario <- factor(MEAN60_OCS$Scenario, levels =c("PT", "GT", 15,30,45,60,75))
MEAN60_OCS$Strategy <- factor(MEAN60_OCS$Strategy, levels =c("SU55", "SU51",  "OCS"))
MEAN60_OCS_abs$Scenario <- factor(MEAN60_OCS_abs$Scenario, levels =c("PT", "GT", 15,30,45,60,75))
MEAN60_OCS_abs$Strategy <- factor(MEAN60_OCS_abs$Strategy, levels =c("SU55", "SU51",  "OCS"))
OCS60$scenario <- factor(OCS60$scenario, levels =c("PT", "GT", 15,30,45,60,75))
OCS60$Strategy <- factor(OCS60$Strategy, levels =c("SU55", "SU51", "SU15", "OCS"))
EFF$Scenario <- factor(EFF$Scenario, levels =c("PT", "GT", 15,30,45,60,75))
EFF$Strategy <- factor(EFF$Strategy, levels =c("SU55", "SU51", "SU15", "OCS"))

OCS60$Name <- paste0(OCS60$Strategy, OCS60$scenario)
OCS60$Name <- factor(OCS60$Name, level=c("SU55PT","SU55GT","SU51GT","OCS15","OCS30","OCS45","OCS60","OCS75"))
MEAN60_OCS[order(MEAN60_OCS$Strategy, MEAN60_OCS$Scenario),]
MEAN60_OCS_abs[,3:ncol(MEAN60_OCS_abs)] <- round(MEAN60_OCS_abs[,3:ncol(MEAN60_OCS_abs)], 2)
MEAN60_OCS_abs[order(MEAN60_OCS_abs$Strategy, MEAN60_OCS_abs$Scenario),]

#significane of genetic gain
library(emmeans)
model <- lm(per_zMean ~ scenario + Strategy + scenario:Strategy, data=OCS60)
anova(model)
marginal = emmeans(model, ~ scenario:Strategy)
CLD = cld(marginal, sort=FALSE, by="Strategy",
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD

#by name = strategy+scenario
model <- lm(per_zMean ~ Name, data=OCS60)
marginal = emmeans(model, ~ Name)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD

#ABSOLUTE#by name = strategy+scenario
model <- lm(zMean ~ Name, data=OCS60)
marginal = emmeans(model, ~ Name)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD


#significance of genic SD
model <- lm(per_GenicSD ~ scenario, data=OCS60)
marginal = emmeans(model, ~ scenario)

#by name = strategy+scenario
model <- lm(per_GenicSD ~ Name, data=OCS60)
marginal = emmeans(model, ~ Name)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
#by name = strategy+scenario
model <- lm(SDGenicSt ~ Name, data=OCS60)
marginal = emmeans(model, ~ Name)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#significance of genetic SD
model <- lm(per_GeneticSD ~ scenario, data=OCS60)
marginal = emmeans(model, ~ scenario)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD


#significance of efficiency - LETTERS
model <- lm(per_Eff ~ Scenario, data=EFF)
marginal = emmeans(model, ~ Scenario)

#by name = strategy+scenario
model <- lm(per_Eff ~ Name, data=EFF)
marginal = emmeans(model, ~ Name)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#ABSOLUTE#by name = strategy+scenario
model <- lm(Slope ~ Name, data=EFF)
marginal = emmeans(model, ~ Name)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

######################################################
######################################################
######################################################
######################################################


#write.csv(Averages[Averages$Generation==60,], "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//Averages_3strategies_14082018.csv", quote=FALSE)
#AveragesA <- Averages
#to je max min za gensko varianco standardizirano
maxmin <- data.frame(strategy=NA, scenario=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (strategy in unique(Averages$strategy)) {
  for (scenario in unique(Averages$scenario)) {
    maxmin[row,] <- c(strategy, scenario, min(Averages$SDGenicSt[(Averages$strategy==strategy) & (Averages$scenario==scenario)]), 
                      max(Averages$SDGenicSt[(Averages$strategy==strategy) & (Averages$scenario==scenario)]), 
                      min(Averages$MeanGenic[(Averages$strategy==strategy) & (Averages$scenario==scenario)]), 
                      max(Averages$MeanGenic[(Averages$strategy==strategy) & (Averages$scenario==scenario)]))
    row <- row +1
  }
}
maxmin <- maxmin[maxmin$minGenicSD != Inf,]
#to je max min za genetsko varianco standardizirano - ampak ne pada monotono --> problemi pri učnkovitost
#maxmin <- data.frame(strategy=NA, scenario=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
#row = 1
#for (strategy in unique(Averages$strategy)) {
#  for (scenario in unique(Averages$scenario)) {
#    maxmin[row,] <- c(strategy, scenario, min(Averages$SDGenetic[(Averages$strategy==strategy) & (Averages$scenario==scenario)]), 
#                      max(Averages$SDGenetic[(Averages$strategy==strategy) & (Averages$scenario==scenario)]), 
#                      min(Averages$MeanGenic[(Averages$strategy==strategy) & (Averages$scenario==scenario)]), 
#                      max(Averages$MeanGenic[(Averages$strategy==strategy) & (Averages$scenario==scenario)]))
#    row <- row +1
#  }
#}



#genska
maxmin$minGenicSD <- as.numeric(maxmin$minGenicSD)
maxmin$maxGenicSD <- as.numeric(maxmin$maxGenicSD)
maxmin$minTGV <- as.numeric(maxmin$minTGV)
maxmin$maxTGV <- as.numeric(maxmin$maxTGV)


####################################################3
library(plyr)
#maxmin$scenario <- revalue(maxmin$scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT",
 #                                             "15"="15", "30"="30", "45"="45", "60"="60", "75"="75"))




maxminOS <- maxmin[maxmin$scenario %in% c(15, 30, 45, 60, 75, "PT", "GT"),]
#maxminOS$scenario <- factor(maxminOS$scenario, levels =c("PT", 15, 30, 45, 60, 75))
maxminOS$PlotGroup <- paste0(maxminOS$strategy, maxminOS$scenario)
maxminPT <- maxminOS[maxminOS$PlotGroup %in% c("SU55PT", "SU55GT", "SU51GT"),]
maxminOCS <- maxminOS[maxminOS$scenario %in% c(15, 30, 45, 60, 75),]
#maxminOS$strategy <- factor(maxminOS$strategy, levels =c("SU55", "OCS"))

#maxminOS <- maxminOS[order(maxminOS$strategy, maxminOS$scenario),]


TGVsStrategy <- TGVsAll[TGVsAll$PlotGroup %in% c("SU55PT", "SU55GT", "SU51GT", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75"),]
TGVsStrategy$strategy <- factor(TGVsStrategy$strategy, levels =c("SU55", "SU51", "OCS"))
TGVsStrategy$scenario <- factor(TGVsStrategy$scenario, levels =c("PT","GT", 15, 30, 45, 60, 75))
TGVsStrategy$PlotGroup <- factor(TGVsStrategy$PlotGroup, level= c("SU55PT", "SU55GT", "SU51GT", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75"))

TGVsStrategy <- TGVsStrategy[order(TGVsStrategy$PlotGroup),]

# TGVsStrategy$PlotGroup <- factor(TGVsStrategy$PlotGroup, level=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75", "SU55PT", "SU51PT", "SU15PT"))
# TGVsStrategy <- TGVsStrategy[order(TGVsStrategy$PlotGroup),]
# maxminPT$PlotGroup <- factor(maxminPT$PlotGroup, level=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75", "SU55PT", "SU51PT", "SU15PT"))
# maxminPT <- maxminPT[order(maxminPT$PlotGroup),]
# maxminOCS$PlotGroup <- factor(maxminOCS$PlotGroup, level=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75", "SU55PT", "SU51PT", "SU15PT"))
# maxminOCS <- maxminOCS[order(maxminOCS$PlotGroup),]

TGVsStrategy$scenario <- as.factor(as.character(TGVsStrategy$scenario))
maxminOS$scenario <- as.factor(as.character(maxminOS$scenario))
#TGVsAll$strategy <- as.factor(TGVsAll$strategy)
#plot za OCS
#library(ggplot2)

summary(TGVsStrategy$SDGenicSt)
summary(TGVsStrategy$zMeanGenic)
table(TGVsStrategy$PlotGroup)
summary(TGVsStrategy$PlotGroup)
TGVsStrategy$PlotGroup <- as.character(TGVsStrategy$PlotGroup)
ggplot(data = TGVsStrategy, aes(x=SDGenicSt, y=zMeanGenic, group=group, colour=PlotGroup)) + 
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                    name="Converted/Lost genic standard deviation")) +
  geom_line(size=0.2, alpha=0.2) + 
  ylim(0,7) + coord_cartesian(xlim = c(1, 0.75)) + theme_bw() +
  scale_linetype_manual("Breeding program", 
                        breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                        values=c("F1", "longdash",  "dashed", "longdash", "twodash", "solid", "solid", "solid"),
                        labels=c(expression("OCS"["15"]), expression("OCS"["30"]), expression("OCS"["45"]),
                                                                                            expression("OCS"["60"]), expression("OCS"["75"]) ,
                                 "5 sires/year, use 5 years, PT", "5 sires/year, use 5 years, GT", "5 sires/year, use 1 year, GT")) + 
  scale_colour_manual("Breeding program", 
                      breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                      values=c("forestgreen", "orange", "purple", "darkblue", "red3","grey55", "grey35", "black"),
                      labels=c(expression("OCS"["15"]), expression("OCS"["30"]), expression("OCS"["45"]),
                               expression("OCS"["60"]), expression("OCS"["75"]) ,
                               "5 sires/year, use 5 years, PT", "5 sires/year, use 5 years, GT", "5 sires/year, use 1 year, GT")) +
  guides(linetype=guide_legend(nrow=3,  label.position = "right", keyheight = unit(1, "cm"), keywidth = unit(3, "cm"), override.aes = list(alpha = 1, size=2))) +
  xlab("Genic standard deviation") + ylab("Genetic Mean") + 
  theme(axis.text=element_text(size=18), legend.position = "top", 
        axis.title.x=element_text(size=18, vjust=-1), 
        axis.title.y=element_text(size=18, vjust=2), 
        legend.text=element_text(size=18), legend.title=element_text(size=18), legend.box = "horizontal",
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        legend.text.align = 0) +
  geom_segment(data=maxminPT, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                         y=minTGV,  yend=maxTGV,                                    
                                         color=PlotGroup, linetype=PlotGroup, group=PlotGroup), arrow=arrow(type="closed"), show.legend=FALSE, size=1.5) +
  geom_segment(data=maxminOCS, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                         y=minTGV,  yend=maxTGV,                                    
                                         color=PlotGroup, linetype=PlotGroup, group=PlotGroup), arrow=arrow(), show.legend=FALSE, size=1.5) 


# library(gridExtra)
library(grid)
library(gtable)
#Tole je bolj na majavih tleh
maxmin$slope <- (maxmin$maxTGV - maxmin$minTGV) / (maxmin$maxGenicSD - maxmin$minGenicSD)
mA <- aggregate(maxmin$slope ~ maxmin$strategy, FUN="summary")
aggregate(maxmin$minGenicSD ~ maxmin$strategy, FUN="summary")


legend = gtable_filter(ggplotGrob(plotList[[1]]), "guide-box") 

grid.arrange(arrangeGrob(plotList[[1]] + theme(legend.position="none"), 
                         plotList[[2]] + theme(legend.position="none"),
                         plotList[[3]] + theme(legend.position="none"),
                         nrow = 1,
                         top = textGrob("Converted/Lost genic standard deviation", vjust = 1, gp = gpar(fontface = "bold", cex = 1.3)),
                         left = textGrob("Average True Genetic Value", rot = 90, vjust = 1, gp = gpar(fontface = "bold", cex = 1.3)), 
                         bottom = textGrob("Genic standard deviation", rot = 0, vjust = -0.1, gp = gpar(fontface = "bold", cex = 1.3))), 
             legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
             nrow=1)

#this is for the ribbon of SD
MeanAverage <- aggregate(TGVsAll$zMean ~ TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
colnames(MeanAverage) <- c("Strategy", "scenario", "Generation", "MeanTGV")
MeanAverage[MeanAverage$Generation==60,]
MeanAverage$PlotGroup <- paste0(MeanAverage$Strategy, MeanAverage$scenario)
MeanAverage <- MeanAverage[MeanAverage$PlotGroup %in% c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),]
MeanAverage$PlotGroup <- factor(MeanAverage$PlotGroup, level=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"))

gen60 <- TGVsAll[TGVsAll$Generation==60,]
gen60 <- summarySE(data=gen60, measurevar = "zMean", groupvars = c("Strategy", "scenario"))

gen60[gen60$Strategy=="SU51" & gen60$scenario=="GT",]
MeanAverage[MeanAverage$Strategy=="SU51" & MeanAverage$scenario=="GT" & MeanAverage$Generation==60,]

MeanAverage$lower <- NA
MeanAverage$upper <- NA

for (row in 1:nrow(MeanAverage)) {
  MeanAverage$lower[row] <- min(TGVsAll$zMean[(TGVsAll$Generation==MeanAverage$Generation[row]) & (TGVsAll$scenario == MeanAverage$scenario[row]) & 
                                        (TGVsAll$strategy == MeanAverage$Strategy[row])])
  MeanAverage$upper[row] <- max(TGVsAll$zMean[(TGVsAll$Generation==MeanAverage$Generation[row]) & (TGVsAll$scenario == MeanAverage$scenario[row]) & 
                                        (TGVsAll$strategy == MeanAverage$Strategy[row])])
}

MeanAverageSD <- aggregate(TGVsAll$SDSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
MeanAverageSDGenic <- aggregate(TGVsAll$SDGenicSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
colnames(MeanAverageSD) <- c("Strategy", "scenario", "Generation", "SdTGV")
colnames(MeanAverageSDGenic) <- c("Strategy", "scenario", "Generation", "SdGenic")

MeanAverageSD$PlotGroup <- paste0(MeanAverageSD$Strategy, MeanAverageSD$scenario)
MeanAverageSD <- MeanAverageSD[MeanAverageSD$PlotGroup %in% c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),]
MeanAverageSD$PlotGroup <- factor(MeanAverageSD$PlotGroup, level=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"))

MeanAverageSDGenic$PlotGroup <- paste0(MeanAverageSDGenic$Strategy, MeanAverageSDGenic$scenario)
MeanAverageSDGenic <- MeanAverageSDGenic[MeanAverageSDGenic$PlotGroup %in% c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),]
MeanAverageSDGenic$PlotGroup <- factor(MeanAverageSDGenic$PlotGroup, level=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"))


avgReg$Eff <- round(avgReg$Slope, 0)
avgReg[avgReg$Strategy=="SU55",]
avgReg[avgReg$Strategy=="SU51",]
avgReg[avgReg$Strategy=="SU15",]


#genGain plot

ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=PlotGroup, linetype=PlotGroup), size=1.2) + 
  ylim(c(0, 7))  + theme_bw() +
  scale_linetype_manual("Breeding program", breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                        values=c("F1", "longdash",  "dashed", "longdash", "twodash", "solid", "solid", "solid"),
                        labels=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT")) + 
  scale_colour_manual("Breeding program", breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                      values=c("forestgreen", "orange", "purple", "darkblue", "red3","grey40", "grey60", "black"),
                      labels=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT")) + 
  guides(linetype=guide_legend(nrow=2, keyheight = unit(1, "cm"), keywidth = unit(3, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  xlab("Genic standard deviation") + ylab("Genetic Mean") + 
  theme(axis.text=element_text(size=16), legend.position = "top", 
        axis.title.x=element_text(size=16, vjust=-1), 
        axis.title.y=element_text(size=16, vjust=2), 
        legend.text=element_text(size=16), legend.title=element_text(size=16), legend.box = "horizontal",
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))


TGVsAll$group <- paste0(TGVsAll$scenario, TGVsAll$Rep)
TGVsAll <- TGVsAll[!(is.na(TGVsAll$scenario)),]
TGVsAll$Rep <- as.factor(TGVsAll$Rep)
#plotaj replike - preveri, katere replike so zalutale
ggplot(data = TGVsAll[TGVsAll$Strategy %in% c("SU55", "OCS"),], aes(x = Generation, y = zMean, colour=scenario, group=group, linetype=Strategy) ) + geom_line()
ggplot(data = TGVsAll[TGVsAll$scenario==15,], aes(x = Generation, y = zMean, colour=Rep, group=Rep) ) + geom_line()
ggplot(data = TGVsAll[TGVsAll$scenario==30,], aes(x = Generation, y = zMean, colour=Rep, group=Rep) ) + geom_line()
ggplot(data = TGVsAll[TGVsAll$scenario==45,], aes(x = Generation, y = zMean, colour=Rep, group=Rep) ) + geom_line()
ggplot(data = TGVsAll[TGVsAll$scenario==60,], aes(x = Generation, y = zMean, colour=Rep, group=Rep) ) + geom_line()
ggplot(data = TGVsAll[TGVsAll$scenario==75,], aes(x = Generation, y = zMean, colour=Rep, group=Rep) ) + geom_line()

gen60 <- TGVsAll[TGVsAll$Generation==60,]
ggplot(data=gen60[gen60$scenario==15,], aes(x=Rep, y=zMean)) + geom_point()
ggplot(data=gen60[gen60$scenario==30,], aes(x=Rep, y=zMean)) + geom_point()
ggplot(data=gen60[gen60$scenario==45,], aes(x=Rep, y=zMean)) + geom_point()
ggplot(data=gen60[gen60$scenario==60,], aes(x=Rep, y=zMean)) + geom_point()
ggplot(data=gen60[gen60$scenario==75,], aes(x=Rep, y=zMean)) + geom_point()



ggplot(data = TGVsAll[TGVsAll$Strategy %in% c("OCS"),], aes(x = Generation, y = zMean, colour=scenario, group=group, linetype=Strategy) ) + geom_line()

#genetic variance plot
MeanAverageSD$Generation <- as.factor(MeanAverageSD$Generation)
ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  geom_line(data = MeanAverageSD, aes(x=Generation, y=SdTGV, colour=PlotGroup, linetype=PlotGroup), size=1.2) + 
  ylim(c(0.85, 1)) +theme_bw() +
  scale_linetype_manual("Breeding program", breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                        values=c("F1", "longdash",  "dashed", "longdash", "twodash", "solid", "solid", "solid"),
                        labels=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT")) + 
  scale_colour_manual("Breeding program", breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                      values=c("forestgreen", "orange", "purple", "darkblue", "red3","grey40", "grey60", "black"),
                      labels=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT")) + 
  guides(linetype=guide_legend(nrow=2, keyheight = unit(1, "cm"), keywidth = unit(3, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  xlab("Genic standard deviation") + ylab("Genetic Mean") + 
  theme(axis.text=element_text(size=16), legend.position = "top", 
        axis.title.x=element_text(size=16, vjust=-1), 
        axis.title.y=element_text(size=16, vjust=2), 
        legend.text=element_text(size=16), legend.title=element_text(size=16), legend.box = "horizontal",
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

#genic variance plot
ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  geom_line(data = MeanAverageSDGenic, aes(x=Generation, y=SdGenic, colour=PlotGroup, linetype=PlotGroup), size=1.2) + 
  ylim(c(0.75, 1)) +theme_bw() +
  scale_linetype_manual("Breeding program", breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                        values=c("F1", "longdash",  "dashed", "longdash", "twodash", "solid", "solid", "solid"),
                        labels=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT")) + 
                        scale_colour_manual("Breeding program", breaks=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT"),
                      values=c("forestgreen", "orange", "purple", "darkblue", "red3","grey40", "grey60", "black"),
                      labels=c("OCS15", "OCS30", "OCS45", "OCS60", "OCS75","SU55PT", "SU55GT", "SU51GT")) +   
  guides(linetype=guide_legend(nrow=2, keyheight = unit(1, "cm"), keywidth = unit(3, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  xlab("Genic standard deviation") + ylab("Genetic Mean") + 
  theme(axis.text=element_text(size=16), legend.position = "top", 
        axis.title.x=element_text(size=16, vjust=-1), 
        axis.title.y=element_text(size=16, vjust=2), 
        legend.text=element_text(size=16), legend.title=element_text(size=16), legend.box = "horizontal",
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))




####
#MST

m <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MST.csv")
mCat <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MSTCat.csv")
m$Scenario <- revalue(m$Scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT"))
m$Strategy <- revalue(m$Strategy, c("SU55" = "5 sires/year, use 5 years","SU51" = "5 sires/year, use 1 year","SU15" = "1 sire/year, use 5 years"))
m$Scenario <- factor(m$Scenario, levels =c("PT", "GT-PT", "GT-C", "GT-BD", "GT"))
m$Strategy <- factor(m$Strategy, levels =c("5 sires/year, use 5 years", "5 sires/year, use 1 year", "1 sire/year, use 5 years"))


table(m$Strategy)
table(m$Scenario, m$Generation) #klasični scenariji imajo poadtke samo do generaije 54 - ker je to zadnja generacije PB, genomski imajo do 58

m$diff <- m$PAmean_PB - m$PAmean
PA <- summarySE(data=m, measurevar = "diff", groupvars = c("Scenario", "Strategy"))[,c(1,2,4,5)]
colnames(PA) <- c("Scenario", "Strategy", "PAdiff_mean", "PAdiff_sd")
MST <- summarySE(data=m, measurevar = "MSTmean_PB", groupvars = c("Scenario", "Strategy"))[,c(1,2,4,5)]
colnames(MST) <- c("Scenario", "Strategy", "MSTpb_mean", "MSTpb_sd")
MSTt <- summarySE(data=m, measurevar = "MSTmean", groupvars = c("Scenario", "Strategy"))[,c(1,2,4,5)]
colnames(MSTt) <- c("Scenario", "Strategy", "MST_mean", "MST_sd")

VSE <- merge(PA, MST, by=c("Scenario", "Strategy"))
VSE <- merge(VSE, MSTt, by=c("Scenario", "Strategy"))
VSE$Scenario <- revalue(VSE$Scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT"))
VSE$Strategy <- revalue(VSE$Strategy, c("SU55" = "5 sires/year, use 5 years","SU51" = "5 sires/year, use 1 year","SU15" = "1 sire/year, use 5 years"))
VSE$Scenario <- factor(VSE$Scenario, levels =c("PT", "GT-PT", "GT-C", "GT-BD", "GT"))
VSE$Strategy <- factor(VSE$Strategy, levels =c("5 sires/year, use 5 years", "5 sires/year, use 1 year", "1 sire/year, use 5 years"))

VSE <- VSE[order(VSE$Strategy, VSE$Scenario),]
VSE[,3:6] <- round(VSE[,3:6],2)
VSE[,7:8] <- round(VSE[,7:8],4)
VSE[,c(1,2,3,5)]
VSE

library(ggplot2)
VSE$Scenario <- as.factor(VSE$Scenario)
PAplot <- ggplot(data=VSE, aes(y=PAdiff_mean, x=Scenario, group=Strategy, fill=Strategy)) + geom_bar(stat="identity", position="dodge") + 
  ylab("PA selected - PA population") + theme_bw() + theme(legend.position="top", legend.text=element_text(size=14), axis.title = element_text(size=14), axis.text = element_text(size=14)) 
MSTplot <- ggplot(data=VSE, aes(y=MSTpb_mean, x=Scenario, group=Strategy, fill=Strategy)) + geom_bar(stat="identity", position="dodge") + 
  ylab("MST selected") + theme_bw() + theme(legend.position="none", axis.title = element_text(size=14), axis.text = element_text(size=14))
#library(Rmisc)
multiplot(PAplot, MSTplot)


#significance PA diff
library(emmeans)
model <- lm(diff ~ Scenario + Strategy + Scenario:Strategy, data=m)
marginal = emmeans(model, ~ Scenario:Strategy)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05, by = "Strategy",
          Letters = letters, adjust="tukey") 
CLD
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05, by = "Scenario",
          Letters = letters, adjust="tukey") 
CLD

#significance MST
library(emmeans)
model <- lm(MSTmean_PB ~ Scenario + Strategy + Scenario:Strategy, data=m)
marginal = emmeans(model, ~ Scenario:Strategy)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05, by = "Strategy",
          Letters = letters, adjust="tukey") 
CLD
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05, by = "Scenario",
          Letters = letters, adjust="tukey") 
CLD
