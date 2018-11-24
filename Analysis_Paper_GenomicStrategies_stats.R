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
######################################################################################################""
setwd("/home/jana/Documents/Projects/inProgress/GenomicStrategies_SireUSe/")

#TGVsAll <- read.csv("~/TGVSALL_11062018.csv")
TGVsAll <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/TGVSALL_14082018.csv")
TGVsAll <- read.csv("TGVSALL_14082018.csv")
#TGVsOCS <- read.csv("~/TGVsAll_OCS_10102018.csv")
#TGVsAll <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_ReferenceSize//Results/TGVSALL_22082018.csv")
TGVsAll$strategy <-TGVsAll$Strategy


#add genetic and genic variance (standardised)
TGVsAll1 <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      base <- TGVsAll$zMean[TGVsAll$scenario=="Class" & TGVsAll$Strategy=="SU55" & TGVsAll$Rep==rep]
      TGVs <- TGVsAll[TGVsAll$Strategy==strategy & TGVsAll$scenario==scenario & TGVsAll$Rep==rep,]
      TGVs$GeneticVarSt <- TGVs$var / TGVs$var[1]
      TGVs$GenicVarSt <- TGVs$AdditGenicVar1 / TGVs$AdditGenicVar1[1]
      
      
      #add percentage
      #     TGVs$per_zMean <- TGVs$zMean / base
      
      TGVsAll1 <- rbind(TGVsAll1, TGVs)  
    }
  }
}

TGVsAll <- TGVsAll1
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
write.csv(Averages[Averages$Generation==60,], "Averages_3strategies_14082018.csv", quote=FALSE)
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
#for (strategy in c("SU55", "SU15", "SU51")) {
strategy=="SU15"
TGVstrategy <- TGVsAll[TGVsAll$Strategy==strategy,]
TGVstrategy$Group <- paste0(TGVstrategy$scenario, TGVstrategy$Rep)
maxminS <- maxmin[maxmin$strategy==strategy,]
maxminSminGenicSD <- as.numeric(maxminS$minGenicSD)
maxminS$maxGenicSD <- as.numeric(maxminS$maxGenicSD)
maxminS$minTGV <- as.numeric(maxminS$minTGV)
maxminS$maxTGV <- as.numeric(maxminS$maxTGV)
  
# plotList[[number]] <- 
ggplot(data = TGVstrategy, aes(x=SDGenicSt, y=zMeanGenic, group=Group, colour=scenario, linetype=scenario)) + 
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                    name="Converted/Lost genic standard deviation")) +
  geom_line(aes(linetype=scenario), size=0.5, alpha=0.2) + #ggtitle(unique(TGVstrategy$strategy)) + 
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
  theme(axis.text=element_text(size=16), legend.position = "left", 
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16), legend.title=element_text(face="bold", size=16)) +
  geom_segment(data=maxminS, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=TRUE, size=1.5)
# number <-  number + 1
# }


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

#standard deviations
MeanAverageSD <- aggregate(TGVsAll$SDSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
MeanAverageSDGenic <- aggregate(TGVsAll$SDGenicSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
MeanAverageSDsd <- aggregate(TGVsAll$SDSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="sd")
MeanAverageSDGenicsd <- aggregate(TGVsAll$SDGenicSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="sd")

colnames(MeanAverageSD) <- c("Strategy", "scenario", "Generation", "GeneticSD")
colnames(MeanAverageSDGenic) <- c("Strategy", "scenario", "Generation", "GenicSD")
colnames(MeanAverageSDsd) <- c("Strategy", "scenario", "Generation", "GeneticSD_SD")
colnames(MeanAverageSDGenicsd) <- c("Strategy", "scenario", "Generation", "GenicSD_SD")

GeneticSD <- merge(MeanAverageSD,MeanAverageSDsd, by=c("Strategy", "scenario", "Generation"))
GenicSD <- merge(MeanAverageSDGenic,MeanAverageSDGenicsd, by=c("Strategy", "scenario", "Generation"))
StandardDeviation <- merge(GeneticSD, GenicSD, by=c("Strategy", "scenario", "Generation"))

TGVsAll <- TGVsAll1
#variances
MeanAverageVARGenetic <- aggregate(TGVsAll$GeneticVarSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
MeanAverageVARGeneticSD <- aggregate(TGVsAll$GeneticVarSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="sd")
MeanAverageVARGenic <- aggregate(TGVsAll$GenicVarSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
MeanAverageVARGenicSD <- aggregate(TGVsAll$GenicVarSt ~  TGVsAll$strategy + TGVsAll$scenario + TGVsAll$Generation, FUN="sd")

colnames(MeanAverageVARGenetic) <- c("Strategy", "scenario", "Generation", "GeneticVAR")
colnames(MeanAverageVARGeneticSD) <- c("Strategy", "scenario", "Generation", "GeneticVAR_SD")
colnames(MeanAverageVARGenic) <- c("Strategy", "scenario", "Generation", "GenicVAR")
colnames(MeanAverageVARGenicSD) <- c("Strategy", "scenario", "Generation", "GenicVAR_SD")


GeneticVar <- merge(MeanAverageVARGenetic,MeanAverageVARGeneticSD, by=c("Strategy", "scenario", "Generation"))
GenicVar <- merge(MeanAverageVARGenic,MeanAverageVARGenicSD, by=c("Strategy", "scenario", "Generation"))
Variance <- merge(GeneticVar, GenicVar, by=c("Strategy", "scenario", "Generation"))


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
    fm1 <- lmList(zMeanGenic ~ SDGenicStNeg | Rep, data=df, pool = TRUE)
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
TGV60$scenario <- as.factor(TGV60$scenario)
TGV60$Strategy <- as.factor(TGV60$Strategy)
TGV60$tr <- paste(TGV60$Strategy, TGV60$scenario, sep="_")
TGV60$tr <- as.factor(TGV60$tr)
TGV601 <- within(TGV60, scenario <- relevel(scenario, ref = "Class"))

dataS = TGV60[TGV60$Strategy=="SU55",]
dataS$Strategy <- as.factor(dataS$Strategy)
dataS$scenario <- as.factor(dataS$scenario)
m1 <- lm(TGV60$zMean ~ TGV60$tr)
m1 <- lm(zMean ~ scenario + Strategy, data=TGV60)
K1 <- glht(m1, mcp(scenario = "Tukey"))$linfct
K2 <- glht(m1, mcp(Strategy = "Tukey"))$linfct
summary(glht(m1, linfct = rbind(K1, K2)))

a1 <- aov(zMean ~ scenario, data=dataS)
posthoc <- TukeyHSD(x=a1, 'scenario', conf.level=0.95)
posthoc
posthoc <- TukeyHSD(x=a1, 'TGV601$Strategy', conf.level=0.95)


model <- lm(zMean ~ scenario + Strategy + scenario : Strategy, data=TGV60)
library(car)
Anova(model,type = "II")
x = residuals(model)
plotNormalHistogram(x)
library(multcompView)
library(lsmeans)
marginal = lsmeans(model, ~ Strategy)


a1 <- aov(m1)
dht <- glht(a1, linfct = mcp(Strategy = "Tukey"))
confint(dht)
summary(dht, test = adjusted("Shaffer"))

### same as TukeyHSD
TukeyHSD(amod, "scenario:Strategy")
### set up linear hypotheses for all-pairs of both factors


tmp <- expand.grid(Strategy = unique(TGV60$Strategy),
                   scenario = unique(TGV60$scenario))
X <- model.matrix(~ scenario * Strategy, data = tmp)
glht(m1, linfct = X)

Tukey <- contrMat(table(TGV60$Strategy), "Tukey")
K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
rownames(K1) <- paste(levels(TGV60$Strategy)[1], rownames(K1), sep = ":")
K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
rownames(K2) <- paste(levels(TGV60$scenario)[2], rownames(K2), sep = ":")
K <- rbind(K1, K2)
colnames(K) <- c(colnames(Tukey), colnames(Tukey))
summary(glht(m1, linfct = K %*% X))


library(agricolae)
print(HSD.test(m1, 'TGV601$Strategy'))
#m1 <- lm(zMean~Strategy,data=TGV601[TGV601$scenario=="Class",])
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "scenario")
m1S <- lsmeans(m1.grid, "Strategy")
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
tgv60 <- data.frame()
EFF <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      #genetic gain
      base <- TGV60$zMean[TGV60$scenario=="Class" & TGV60$Strategy=="SU55" & TGV60$Rep==rep]
      tgv <- TGV60[TGV60$scenario==scenario & TGV60$Strategy==strategy & TGV60$Rep==rep,]
      tgv$per_zMean <- tgv$zMean / base
      
      #genic variance
      base <- TGV60$GenicVarSt[TGV60$scenario=="Class" & TGV60$Strategy=="SU55" & TGV60$Rep==rep]
      tgv$per_GenicVar <- (tgv$GenicVarSt / base)
      
      #genetic variance
      base <- TGV60$GeneticVarSt[TGV60$scenario=="Class" & TGV60$Strategy=="SU55" & TGV60$Rep==rep]
      tgv$per_GeneticVar <- (tgv$GeneticVarSt / base)      
      
      effBase <- regRep$Slope[regRep$Scenario=="Class"  & regRep$Strategy=="SU55" & regRep$Rep==rep]
      eff <- regRep[regRep$Scenario==scenario & regRep$Strategy==strategy & regRep$Rep==rep,]
      eff$per_Eff <- (eff$Slope / effBase)
      
      tgv60 <- rbind(tgv60, tgv)
      EFF <- rbind(EFF, eff)
    }
  }
}
TGV60 <- tgv60

TGV60$per_zMean <- TGV60$per_zMean*100 - 100 #tukaj ma osnoven scenarij 100: zato odštej 100
TGV60$per_GenicVar <- (1 - TGV60$per_GenicVar)*100 #tukaj ma osnoven scenarij 0
TGV60$per_GeneticVar <- (1-TGV60$per_GeneticVar)*100
EFF$per_Eff <-  EFF$per_Eff * 100 -100

#genetic gain
MEAN60_abs <- summarySE(TGV60, measurevar="zMean", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5,6)]
MEAN60 <- summarySE(TGV60, measurevar="per_zMean", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5,6)]
colnames(MEAN60) <- c("Strategy", "Scenario", "per_zMean", "per_zMean_qL", "per_zMean_qH")
colnames(MEAN60_abs) <- c("Strategy", "Scenario", "zMean", "zMean_qL", "zMean_q")
MEAN60$per_zMean <- round(MEAN60$per_zMean)
MEAN60$per_zMean_qL <- round(MEAN60$per_zMean_qL)
MEAN60$per_zMean_qH <- round(MEAN60$per_zMean_qH)
head(MEAN60)


#genetic variance
VAR60a <- summarySE(TGV60, measurevar="per_GeneticVar", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(VAR60a) <- c("Strategy", "Scenario", "per_GeneticVar", "per_GeneticVarSD")
VAR60a$per_GeneticVar <- round(VAR60a$per_GeneticVar, 3)
VAR60a$per_GeneticVarSD <- round(VAR60a$per_GeneticVarSD, 3)

MEAN60 <- merge(MEAN60, VAR60a, by=c("Strategy", "Scenario"))

#genic variance
VAR60b <- summarySE(TGV60, measurevar="per_GenicVar", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(VAR60b) <- c("Strategy", "Scenario", "per_GenicVar", "per_GenicVarSD")
VAR60b$per_GenicVar <- round(VAR60b$per_GenicVar, 3)
VAR60b$per_GenicVarSD <- round(VAR60b$per_GenicVarSD, 3)

MEAN60 <- merge(MEAN60, VAR60b, by=c("Strategy", "Scenario"))


#efficiency
eff <- summarySE(EFF, measurevar="per_Eff", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
colnames(eff) <- c("Strategy", "Scenario", "per_Eff", "per_EffSD")
eff$per_Eff <- round(eff$per_Eff)
eff$per_EffSD <- round(eff$per_EffSD)

MEAN60 <- merge(MEAN60, eff, by=c("Strategy", "Scenario"))

#for displaying
MEAN60$Strategy <- factor(MEAN60$Strategy, levels =c("SU55", "SU51", "SU15"))
MEAN60$Scenario <- factor(MEAN60$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
MEAN60[order(MEAN60$Strategy, MEAN60$Scenario),]

MEAN60[MEAN60$Strategy=="SU55",]
MEAN60[MEAN60$Strategy=="SU15",]
MEAN60[MEAN60$Strategy=="SU51",]
MEAN60$per_zMean[MEAN60$Strategy=="SU51"] - MEAN60$per_zMean[MEAN60$Strategy=="SU55"]
MEAN60$per_zMean[MEAN60$Strategy=="SU15"] - MEAN60$per_zMean[MEAN60$Strategy=="SU55"]
write.csv(MEAN60, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//StandardDeviation_GeneticGain_gen60_16112018.csv", quote=FALSE)



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
