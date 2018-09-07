qstatacc <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/AccuraciesALL_correlation.csv")
#acc <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/AccuraciesALL_Cat.csv")
#accAcc <- acc[,paste0("X", 21:40)]
#accNames <- acc[,c("cat", "Strategy", "Scenario", "Rep")]
#accNames$corEBV <- as.numeric(rowMeans(accAcc))
#acc <- accNames
accA <- aggregate(acc$corEBV ~ acc$Scenario + acc$Strategy + acc$Cat + acc$Cycle, FUN="mean")
accASc <- aggregate(acc$corEBV ~ acc$Strategy + acc$Cat + acc$Cycle, FUN="mean")
accSc <- aggregate(acc$corEBV ~ acc$Scenario + acc$Strategy, FUN="mean")
accSc <- aggregate(acc$corEBV ~ acc$Strategy + acc$Cat, FUN="mean")
accS <- aggregate(acc$corEBV ~ acc$Strategy, FUN="mean")
colnames(accA) <- c("Scenario", "Strategy", "Cat",  "Cycle", "corEBV")
colnames(accSc) <- c("Scenario", "Strategy", "Cat", "corEBV")
colnames(accASc) <- c("Strategy", "Cat", "Cycle", "corEBV")
write.csv(accA, "~/Documents/PhD/Simulaton/Accuracies_genomicStrategies.csv", quote=FALSE)

mean(accA$corEBV[accA$Scenario=="Class" & accA$Cat=="potomciNP"])
mean(accA$corEBV[accA$Cat=="genTest"])

mean(accA$corEBV[accA$Scenario %in% c("Class", "GenSLO") & accA$Cat=="pb"])
mean(accA$corEBV[accA$Cat=="pb"])
mean(accA$corEBV[accA$Cat=="gpb"])
accA[accA$Cat=="gpb",]
acc[acc$Strategy=="SU55" & acc$Scenario=="Class" & acc$Cat=="mladi",]


accASc[accASc$Cat=="genTest",]
ggplot(data=accASc[accASc$Cat=="genTest",], aes(x=Cycle, y=corEBV, group=Strategy, colour=Strategy)) + geom_line()


accA[(accA$Cat %in% c("vhlevljeni", "genTest", "mladi", "potomciNP")) & (accA$Strategy=="SU55"),]
comp <- cbind(accSc[(accSc$Cat %in% c("mladi", "genTest")) & (accSc$Strategy=="SU55"),], accSc[(accSc$Cat %in% c("mladi", "genTest")) & (accSc$Strategy=="SU51"),])
comp$diff <- accSc$corEBV[(accSc$Cat %in% c( "genTest")) & (accSc$Strategy=="SU55")] - accSc$corEBV[(accSc$Cat %in% c( "genTest")) & (accSc$Strategy=="SU51")]
diffK <- accSc$corEBV[(accSc$Cat %in% c( "k")) & (accSc$Strategy=="SU55")] - accSc$corEBV[(accSc$Cat %in% c( "k")) & (accSc$Strategy=="SU51")]
mean(comp$diff)
mean(acc$corEBV[acc$Strategy=="SU55"], na.rm=TRUE) - mean(acc$corEBV[acc$Strategy=="SU51"], na.rm=TRUE)
sd(acc$corEBV[acc$Strategy=="SU55"], na.rm=TRUE) 
sd(acc$corEBV[acc$Strategy=="SU51"], na.rm=TRUE)
######################################################################################################
######################################################################################################
######################################################################################################
#TGVsAll <- read.csv("~/TGVSALL_11062018.csv")
TGVsAll <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/TGVSALL_14082018.csv")
#TGVsAll <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_ReferenceSize//Results/TGVSALL_22082018.csv")
TGVsAll$strategy <-TGVsAll$Strategy
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


write.csv(Averages[Averages$Generation==60,], "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//Averages_3strategies_14082018.csv", quote=FALSE)
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

#o je z na genetsko standadrizirano gensko variacno
library(ggplot2)
plotList = list()
number = 1
for (strategy in c("SU55", "SU15", "SU51")) {
  TGVstrategy <- TGVsAll[TGVsAll$Strategy==strategy,]
  TGVstrategy$Group <- paste0(TGVstrategy$scenario, TGVstrategy$Rep)
  maxminS <- maxmin[maxmin$strategy==strategy,]
  maxminSminGenicSD <- as.numeric(maxminS$minGenicSD)
  maxminS$maxGenicSD <- as.numeric(maxminS$maxGenicSD)
  maxminS$minTGV <- as.numeric(maxminS$minTGV)
  maxminS$maxTGV <- as.numeric(maxminS$maxTGV)
  
  plotList[[number]] <- ggplot(data = TGVstrategy, aes(x=SDGenicSt, y=zMeanGenic, group=Group, colour=scenario, linetype=scenario)) + 
    scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                      name="Converted/Lost genic standard deviation")) +
    geom_line(aes(linetype=scenario), size=0.5, alpha=0.2) + ggtitle(unique(TGVstrategy$strategy)) + 
    xlab("Generation") + ylab("True genetic value")  + ylim(0,7) +coord_cartesian(xlim = c(1, 0.85)) +
    scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                          "Scenario", 
                          values=c("solid", "dotted","dashed", "dotdash", "twodash"), 
                          labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
    scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        "Scenario", 
                        values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                        labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
    xlab("Genic standard deviation") + ylab("Average True Genetic Value") + 
    theme(axis.text=element_text(size=16), legend.position = "left", axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16), legend.title=element_text(face="bold", size=16)) +
    geom_segment(data=maxminS, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                          y=minTGV,  yend=maxTGV,                                    
                                          color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=TRUE, size=1.5) 
  number <-  number + 1
}


library(gridExtra)
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

#genetic gain
library(lme4)
Means <- data.frame()
for (str in c("SU55", "SU51", "SU15")) {
  meanStr <- MeanAverage[MeanAverage$Strategy == str,]
  mean <- lmList(MeanTGV ~ Generation | scenario, data=meanStr)
  meanDF <- data.frame(Scenario=rownames(coef(mean)),coef(mean),check.names=FALSE)
  colnames(meanDF) <- c("Scenario", "Mean_Intercept", "Mean_Slope")
  meanDF$Stratey <- str
  Means <- rbind(Means, meanDF)
}

MeansSD <- data.frame()
for (str in c("SU55", "SU51", "SU15")) {
  meanStr <- MeanAverageSD[MeanAverageSD$Strategy == str,]
  mean <- lmList(SdTGV ~ Generation | scenario, data=meanStr)
  meanDF <- data.frame(Scenario=rownames(coef(mean)),coef(mean),check.names=FALSE)
  colnames(meanDF) <- c("Scenario", "Sd_Intercept", "SD_Slope")
  meanDF$Stratey <- str
  MeansSD <- rbind(MeansSD, meanDF)
}


#to je povrepčje regresij
library(nlme)
regRep <- data.frame(Rep=NA, Intercept=NA, Slope=NA, Scenario=NA, Strategy=NA)
for (str in c("SU55", "SU51", "SU15")) {
  for (sc in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    df <- TGVsAll[(TGVsAll$scenario==sc) & (TGVsAll$strategy==str),]
    fm1 <- lmList(zMeanGenic ~ SDGenicSt | Rep, data=df)
    tmp <- data.frame(Rep=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
    colnames(tmp) <- c("Rep", "Intercept", "Slope")
    tmp$Scenario <- sc
    tmp$Strategy <- str
    regRep <- rbind(regRep, tmp)
  }
}

#povprečje regresijskih koefificentov
avgSlope <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgReg <- cbind(avgInt, avgSlope[,2])
colnames(avgReg) <- c("Scenario", "Strategy",  "Intercept", "Slope")
#and Sd
sdSlope <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="sd")
colnames(sdSlope) <- c("Scenario", "Strategy",  "SDSlope")
write.csv(SD, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//StandardDeviation_Efficiency_14082018.csv", quote=FALSE)

#efficiency
regRep <- data.frame(Rep=NA, Intercept=NA, Slope=NA, Scenario=NA, Strategy=NA)
for (str in c("SU55", "SU51", "SU15")) {
  for (sc in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    df <-  TGVsAll[(TGVsAll$scenario==sc) & (TGVsAll$strategy==str),]
    fm1 <- lmList(zMeanGenic ~ SDGenicStNeg | Rep, data=df)
    tmp <- data.frame(Rep=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
    colnames(tmp) <- c("Rep", "Intercept", "Slope")
    tmp$Scenario <- sc
    tmp$Strategy <- str
    regRep <- rbind(regRep, tmp)
  }
}
regRep <- regRep[-1,]


#povprečje regresijskih koefificentov - efficiency
avgSlope <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgSlopeSD <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="sd")
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgIntSD <- aggregate(regRep$Intercept ~ regRep$Scenario + regRep$Strategy, FUN="sd")
avgReg <- merge(avgSlope, avgSlopeSD, by=c("regRep$Scenario", "regRep$Strategy"))
avgReg <- merge(avgReg, avgInt, by=c("regRep$Scenario", "regRep$Strategy"))
avgReg <- merge(avgReg, avgIntSD, by=c("regRep$Scenario", "regRep$Strategy"))
colnames(avgReg) <- c("Scenario", "Strategy",  "Slope", "SlopeSD", "Intercept", "InterceptSD")
avgReg$Eff <- round(avgReg$Slope, 0)
avgReg$EffSd <- round(avgReg$SlopeSD, 0)
write.csv(avgReg, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//Efficiency_genomicstrategies_14082018.csv", quote=FALSE)

#efficiency of SU55
a <- avgReg[avgReg$Strategy=="SU55",]
a[order(a$Slope),]
a <- avgReg[avgReg$Strategy=="SU15",]
a[order(a$Slope),]
a <- avgReg[avgReg$Strategy=="SU51",]
a[order(a$Slope),]
1 - avgReg$Slope[avgReg$Strategy=="SU15"] / avgReg$Slope[avgReg$Strategy=="SU55"]

#significance of efficiency
library(emmeans)
regRep$Scenario <- as.factor(regRep$Scenario)
regRep$Strategy <- as.factor(regRep$Strategy)
regRep1 <- within(regRep, Scenario <- relevel(Scenario, ref = "Class"))
m1 <- lm(Slope~Scenario,data=regRep1[regRep1$Strategy=="SU55",])
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "Scenario")
contrast(m1.grid, method="pairwise")
contrast(m1S, method="eff")
summary(lm(Slope~Scenario,data=regRep1))

#significance of genetic gain
library(emmeans)
Mean60$Scenario <- as.factor(Mean60$Scenario)
Mean60$Strategy <- as.factor(Mean60$Strategy)
Mean601 <- within(Mean60, Scenario <- relevel(Scenario, ref = "Class"))
m1 <- lm(MeanTGV~Scenario,data=Mean601[Mean601$Strategy=="SU51",])
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "Scenario")
contrast(m1.grid, method="pairwise")
contrast(m1S, method="eff")
summary(lm(Slope~Scenario,data=regRep1))

#vrstni red
#preveri, ali je klasična najslabš v vseh, 1 = najslabša, 5 = najboljša
for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
  vec <- c()
  for (rep in 0:19) {
    repDF <- regRep[(regRep$Rep==rep) & (regRep$Strategy=="SU55"),]
    a <- repDF[order(repDF$Slope),]
    #vec <- c(vec, repDF$Scenario[repDF$Slope == min(repDF$Slope)] == "Class")
    vec <- c(vec, which(a$Scenario==scenario))
  }
  print(scenario)
  print(sum(vec == 1))
}

avgReg$Eff <- round(avgReg$Slope, 0)
avgReg[avgReg$Strategy=="SU55",]
avgReg[avgReg$Strategy=="SU51",]
avgReg[avgReg$Strategy=="SU15",]

#Razlike v genic variance loss
maxmin[maxmin$strategy=="SU15","minGenicSD"] -  maxmin[maxmin$strategy=="SU55","minGenicSD"]

#genGain plot
MeanAverage$order <- factor(MeanAverage$Strategy, levels = c("SU55", "SU15", "SU51"))
ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), 
                        labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  geom_line(data = MeanAverage[MeanAverage$Strategy=="SU55",], aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
  #geom_ribbon(data=MeanAverage, aes(x=Generation, ymin=lower, ymax=upper, colour=scenario), alpha=0.1) + 
ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=18), legend.title=element_text(face="bold", size=18), 
        strip.text = element_text(face="bold", size=16))   #+ 
facet_grid(order ~ ., scales = "free_y") + theme(legend.position = "right") 

#genetic variance plot
MeanAverageSD$order <- factor(MeanAverageSD$Strategy, levels = c("SU55", "SU15", "SU51"))

ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), 
                        labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  geom_line(data = MeanAverageSD, aes(x=Generation, y=SdTGV, colour=scenario, linetype=scenario), size=1.2) + 
  ylim(c(0.85, 1)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=16), legend.title=element_text(face="bold", size=18), 
        strip.text = element_text(face="bold", size=16))   + 
  facet_grid(order ~ ., scales = "free_y") + theme(legend.position = "right") 


MeanAverageSDGenic$order <- factor(MeanAverageSDGenic$Strategy, levels = c("SU55", "SU15", "SU51"))
#genic variance plot
ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), 
                        labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  geom_line(data = MeanAverageSDGenic, aes(x=Generation, y=SdGenic, colour=scenario, linetype=scenario), size=1.2) + 
  ylim(c(0.85, 1)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=16), legend.title=element_text(face="bold", size=18), 
        strip.text = element_text(face="bold", size=16))   + 
  facet_grid(order ~ ., scales = "free_y") + theme(legend.position = "right") 


#standard deviations of the measures
#genetic gain in gen 60
TGV60 <- TGVsAll[TGVsAll$Generation==60,]
Mean60 <- aggregate(TGV60$zMean ~ TGV60$strategy + TGV60$scenario + TGV60$Rep, FUN="mean")
colnames(Mean60) <- c("Strategy", "Scenario", "Rep", "MeanTGV")
MEAN60 <- aggregate(Mean60$MeanTGV ~Mean60$Strategy + Mean60$Scenario, FUN="mean")
Sd60 <- aggregate(Mean60$MeanTGV ~Mean60$Strategy + Mean60$Scenario, FUN="sd")
colnames(Sd60) <- c("Strategy", "Scenario", "SdTGV")
colnames(MEAN60) <- c("Strategy", "Scenario", "MeanTGV")

SD <- merge(MEAN60, Sd60, by=c("Strategy", "Scenario")) 
SD$percentage <- round(SD$MeanTGV / SD$MeanTGV[SD$Strategy=="SU55" & SD$Scenario=="Class"]*100,0)
SD[SD$Strategy=="SU55",]
SD[SD$Strategy=="SU15",]
SD[SD$Strategy=="SU51",]
SD$percentage[SD$Strategy=="SU51"] - SD$percentage[SD$Strategy=="SU55"]
write.csv(SD, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//StandardDeviation_GeneticGain_gen60_19072018.csv", quote=FALSE)



'
####
pedQ <- read.table("~/PedQTN.txt", header=TRUE)
varQ <- read.table("~/VarQTN.txt", header=TRUE)


TGVs <- data.frame(Generation=40:60)
ped <- pedQ
#to standardise onto the generation 40 - which is the generation of comparison
ped <- ped[ped$Generation %in% 40:60,]
#obtain mean and sd of genetic values
TGV <- summarySE(ped, measurevar = "gvNormUnres1", groupvars = "Generation")[,c(1,3,4)]
#variance of genetic values
TGV$var <- (TGV$sd)^2
colnames(TGV)[1] <- c("Generation")
#standardise genetic standard devistion
TGV$SDSt <- TGV$sd / TGV$sd[1]
#standardise genetic values with genetic standard deviation
TGV$zMean <- (TGV$gvNormUnres1 - TGV$gvNormUnres1[1]) / TGV$sd[1]
TGVs <- merge(TGVs, TGV, by="Generation")
#read in genic variance
Var <-varQ
#Qtn model 1 is unrestricted
Var <- Var[Var$QtnModel==1,c(1,3)]
TGVs <- merge(TGVs, Var, by="Generation")
#obtain genic standard deviation
TGVs$SDGenic <- (sqrt(TGVs$AdditGenicVar1))
#standarise genic standard devistion
TGVs$SDGenicSt <- TGVs$SDGenic / TGVs$SDGenic[1]
#standardise genetic values with genic standard devistion
TGVs$zMeanGenic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$SDGenic[1]
#reciprocated genic standard deviation
TGVs$SDGenicStNeg <- 1 - (TGVs$SDGenic / TGVs$SDGenic[1])
#genic variance standardised onto genetic variance
koef <- TGVs$var[1] / TGVs$AdditGenicVar1[1]
TGVs$Genic_Genetic_VAR <- TGVs$AdditGenicVar1 * koef
TGVs$Genic_Genetic_SD <- sqrt(TGVs$Genic_Genetic_VAR)
#standardise genic_genetic standard deviation
TGVs$Genic_Genetic_SDSt <- TGVs$Genic_Genetic_SD / TGVs$Genic_Genetic_SD[1]
#standarise genetic values with genic_genetic standard deviation
TGVs$zMeanGenic_Genetic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$Genic_Genetic_SD[1]
#TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
TGVs$scenario <- "Gen_noQTN"
TGVs$Rep <- 21
TGVs$Strategy <- "10K_Ref_20Rep"

#TGVs Gen0_noQTN
TGVQ <- TGVsAll[(TGVsAll$Rep %in% 0:19) & (TGVsAll$scenario=="Gen") & (TGVsAll$Strategy=="10K_Ref_20Rep"),]

TGVc <- rbind(TGVQ, TGVs)
ggplot(data=TGVc, aes(x=Generation, y=zMean, group=Rep, fill=Rep, colour=scenario))  + geom_path() + #y=SDGenicSt, 
  scale_colour_manual("Scenario", breaks = c("Gen_noQTN", "Gen"), values=c("forestgreen", "red"), labels=c("Gen_noQTNsOnChip", "Gen")) + 
  xlab("Generation") + ylab("Genic Standard Deviation") + #ylab("Average True Genetic Value") +
  theme(axis.text=element_text(size=16), legend.position = "left",
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16), legend.title=element_text(size=16))  
'