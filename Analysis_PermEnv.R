
################
#Funkcija summarySE
#####################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[colnames(datac) == "mean"] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################

acc <- read.csv("~/Documents/PhD/Projects/inProgress/AmountOfPhenotyping//Results/ACCAge.csv")
library(plyr)
acc$scenario <-  revalue(acc$scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT"))
acc$scenario <- factor(acc$scenario, levels =c("PT", "GT-PT", "GT-C", "GT-BD", "GT"))
acc <- acc[order(acc$scenario),]
#acc <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/AccuraciesALL_Cat.csv")
#accAcc <- acc[,paste0("X", 21:40)]
#accNames <- acc[,c("cat", "Strategy", "Scenario", "Rep")]
#accNames$corEBV <- as.numeric(rowMeans(accAcc))
#acc <- accNames
accA <- summarySE(data=acc, measurevar = "COR" , groupvars = c("scenario", "strategy", "AgeCat"))[,c(1,2,3,5,6)]
colnames(accA) <- c("Scenario", "Strategy", "CatAge", "meanAcc", "sdAcc")
table(accA$CatAge)

accA[accA$Strategy=="SU55" & accA$CatAge %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]
accA[accA$Strategy=="SU51" & accA$CatAge %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]
accA[accA$Strategy=="SU15" & accA$CatAge %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]

accSig <- acc[acc$AgeCat %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]


#significance of significances
library(emmeans)
model <- lm(COR ~ strategy + scenario + AgeCat + strategy : scenario : AgeCat, data=accSig)
marginal = emmeans(model, ~ strategy:scenario:AgeCat)
CLD = cld(marginal, by=c("AgeCat", "strategy"),
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05, sort=FALSE,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#mean za gentest in cak5
mean(acc$COR[acc$AgeCat == "cak5"])
mean(acc$COR[acc$AgeCat == "genTest1"])

bias <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Bias_CatAge.csv")[,-1]
colnames(bias) <- c("Intercept", "Bias", "strategy", "scenario", "rep", "AgeCat")

bias <- bias[bias$AgeCat %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]

biasA <- summarySE(data=bias, measurevar = "Bias" , groupvars = c("scenario", "strategy", "AgeCat"))[,c(1,2,3,5,6)]
colnames(biasA) <- c("Scenario", "Strategy", "CatAge", "meanBias", "sdBias")
table(biasA$CatAge)

biasA[biasA$Strategy=="SU55" & biasA$CatAge %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]
biasA[biasA$Strategy=="SU51" & biasA$CatAge %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]
biasA[biasA$Strategy=="SU15" & biasA$CatAge %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]

biasSig <- bias[bias$AgeCat %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]
######################################################################################################
######################################################################################################
######################################################################################################""
setwd("/home/jana/Documents/PhD/Projects/inProgress/AmountOfPhenotyping//")

#TGVsAll <- read.csv("~/TGVSALL_11062018.csv")
TGVsAll <- read.table("~/Documents/PhD/Projects/inProgress/AmountOfPhenotyping//Results/TGVsAll_permEnv_SU55_17072019.csv", header=TRUE)
TGVsAll <- read.table("~/TGVsAll_permEnv_SU55_28072019.csv", header=TRUE)
#TGVsAll <- read.csv("~/Documents/Projects/inProgress/GenomicStrategies_SireUSe/TGVSALL_14082018.csv")
#TGVsOCS <- read.csv("~/TGVsAll_OCS_10102018.csv")
#TGVsAll <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_ReferenceSize//Results/TGVSALL_22082018.csv")
TGVsAll$strategy <-TGVsAll$Strategy

agg <- summarySE(data = TGVsAll, measurevar = "zMean", groupvars = c("scenario", "Repeats", "Generation"))
agg$GROUP <- paste0(agg$scenario, agg$Repeats)
agg$Repeats <- as.factor(agg$Repeats)
library(ggplot2)
ggplot(data = agg, aes(x=Generation, y=zMean, group=GROUP, colour=Repeats, linetype=scenario)) + geom_line()

#add genetic and genic variance (standardised)
TGVsAll1 <- data.frame()
for (strategy in c(1,2,5,8, 9, 11)) {
  for (scenario in c("Class", "Gen")) {
    for (rep in 0:0) {
      TGVs <- TGVsAll[TGVsAll$Repeats==strategy & TGVsAll$scenario==scenario & TGVsAll$Rep==rep,]
      TGVs$GeneticVarSt <- TGVs$var / TGVs$var[1]
      TGVs$GenicVarSt <- TGVs$AdditGenicVar1 / TGVs$AdditGenicVar1[1]
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


#write.csv(Averages[Averages$Generation==60,], "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//Averages_3strategies_14082018.csv", quote=FALSE)
#write.csv(Averages[Averages$Generation==60,], "Averages_3strategies_14082018.csv", quote=FALSE)
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


TGVsAll$Repeats <- as.factor(TGVsAll$Repeats)
ggplot(data=TGVsAll, aes(x=Generation, y=zMean, group=Repeats, colour=Repeats)) + geom_line()

#o je z na genetsko standadrizirano gensko variacno
library(ggplot2)
plotList = list()
number = 1
for (strategy in c("SU55", "SU51", "SU15")) {
  strategy==strategy
  TGVstrategy <- TGVsAll[TGVsAll$Strategy==strategy,]
  
  TGVstrategy$Group <- paste0(TGVstrategy$scenario, TGVstrategy$Rep)
  TGVstrategy$Strategy <- factor(TGVstrategy$Strategy, levels =c("SU55", "SU51", "SU15"))
  TGVstrategy$scenario <- factor(TGVstrategy$scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
  TGVstrategy <- TGVstrategy[order(TGVstrategy$Strategy, TGVstrategy$scenario),]
  maxminS <- maxmin[maxmin$strategy==strategy,]
  maxminSminGenicSD <- as.numeric(maxminS$minGenicSD)
  maxminS$maxGenicSD <- as.numeric(maxminS$maxGenicSD)
  maxminS$minTGV <- as.numeric(maxminS$minTGV)
  maxminS$maxTGV <- as.numeric(maxminS$maxTGV)
  
  STRATEGY <- ifelse(strategy == "SU55",'5 sires/year, use 5 years', ifelse(strategy =="SU51", "5 sires/year, use 1 year", "1 sire/year, use 5 years" ))
  
  plotList[[number]] <- 
  ggplot(data = TGVstrategy, aes(x=SDGenicSt, y=zMeanGenic, group=Group, colour=scenario, linetype=scenario)) + 
    scale_x_reverse(sec.axis=sec_axis(trans=~ . -1.)) +                                   
                                      #name="Converted/Lost genic standard deviation")) +
    geom_line(aes(linetype=scenario), size=0.2, alpha=0.4) + ggtitle(STRATEGY) + 
        ylim(0,7) + coord_cartesian(xlim = c(1, 0.85)) + theme_bw() +
    scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                          "Breeding program", 
                          values=c("solid", "longdash","dashed", "twodash", "F1"), 
                          labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
    scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        "Breeding program", 
                        values=c("black", "forestgreen", "dodgerblue2", "purple", "red3"), 
                        labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
    guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(3, "cm"), override.aes = list(alpha = 1, size=1.2))) +
   # xlab("Genic standard deviation") + ylab("Average True Genetic Value") + 
    theme(axis.text=element_text(size=16), legend.position = "top", 
          axis.title=element_blank(), legend.text=element_text(size=16), legend.title=element_text(size=16),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
          plot.margin = margin(t = 0, r = 10, b = 10, l = 10)) +
    geom_segment(data=maxminS, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                          y=minTGV,  yend=maxTGV,                                    
                                          color=scenario, linetype=scenario, group=scenario),  arrow=arrow(length=unit(0.7, 'cm')), 
                                          show.legend=FALSE, size=1.4, alpha=1)
  number <-  number + 1
}


library(gridExtra)
library(grid)
library(gtable)
#Tole je bolj na majavih tleh
#maxmin$slope <- (maxmin$maxTGV - maxmin$minTGV) / (maxmin$maxGenicSD - maxmin$minGenicSD)
#mA <- aggregate(maxmin$slope ~ maxmin$strategy, FUN="summary")
#aggregate(maxmin$minGenicSD ~ maxmin$strategy, FUN="summary")
# 
#install.packages("ggpubr")
library(ggpubr)

#legend = gtable_filter(ggplotGrob(plotList[[1]]), "axis-b") 

figure <- ggarrange(plotList[[1]] + theme(legend.position="none"), 
          plotList[[2]] + theme(legend.position="none"),
          plotList[[3]] + theme(legend.position="none"),
          common.legend = TRUE, legend="top", nrow=1, ncol=3, vjust=-2)

annotate_figure(figure, 
                top = textGrob("Converted/Lost genic standard deviation", vjust = 7.1, gp = gpar(cex = 1.5)),
                left = textGrob("Genetic mean", rot = 90, vjust = 0.7, hjust=0.8, gp = gpar(cex = 1.5)), 
                bottom = textGrob("Genic standard deviation", rot = 0, vjust = -0, gp = gpar( cex = 1.5)) 
                )
# grid.arrange(arrangeGrob(plotList[[1]] + theme(legend.position="none"), 
#                          plotList[[2]] + theme(legend.position="none"),
#                          plotList[[3]] + theme(legend.position="none"),
#                          nrow = 1,
#                          top = textGrob("Converted/Lost genic standard deviation", vjust = 1.5, gp = gpar(cex = 1.5)),
#                          left = textGrob("Genetic mean", rot = 90, vjust = 0.8, gp = gpar(cex = 1.5)), 
#                          bottom = textGrob("Genic standard deviation", rot = 0, vjust = -0, gp = gpar( cex = 1.5))), 
#              legend, nrow=2,heights=c(10, 1),
#              widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
#              nrow=1)

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
avGTlope <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgReg <- cbind(avgInt, avGTlope[,2])
colnames(avgReg) <- c("Scenario", "Strategy",  "Intercept", "Slope")
#and Sd
sdSlope <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="sd")
colnames(sdSlope) <- c("Scenario", "Strategy",  "SDSlope")
#write.csv(SD, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//StandardDeviation_Efficiency_14082018.csv", quote=FALSE)

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
avGTlope <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avGTlopeSD <- aggregate(regRep$Slope ~ regRep$Scenario + regRep$Strategy, FUN="sd")
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario + regRep$Strategy, FUN="mean")
avgIntSD <- aggregate(regRep$Intercept ~ regRep$Scenario + regRep$Strategy, FUN="sd")
avgReg <- merge(avGTlope, avGTlopeSD, by=c("regRep$Scenario", "regRep$Strategy"))
avgReg <- merge(avgReg, avgInt, by=c("regRep$Scenario", "regRep$Strategy"))
avgReg <- merge(avgReg, avgIntSD, by=c("regRep$Scenario", "regRep$Strategy"))
colnames(avgReg) <- c("Scenario", "Strategy",  "Slope", "SlopeSD", "Intercept", "InterceptSD")
avgReg$Eff <- round(avgReg$Slope, 0)
avgReg$EffSd <- round(avgReg$SlopeSD, 0)
#write.csv(avgReg, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//Efficiency_genomicstrategies_14082018.csv", quote=FALSE)




#efficiency of SU55
a <- avgReg[avgReg$Strategy=="SU55",]
a[order(a$Slope),]
a <- avgReg[avgReg$Strategy=="SU15",]
a[order(a$Slope),]
a <- avgReg[avgReg$Strategy=="SU51",]
a[order(a$Slope),]
1 - avgReg$Slope[avgReg$Strategy=="SU15"] / avgReg$Slope[avgReg$Strategy=="SU55"]


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
      
      #genic standard
      base <- TGV60$SDGenicSt[TGV60$scenario=="Class" & TGV60$Strategy=="SU55" & TGV60$Rep==rep]
      tgv$per_GenicSD <- (tgv$SDGenicSt / base)
      
      #genetic standard deviation
      base <- TGV60$SDSt[TGV60$scenario=="Class" & TGV60$Strategy=="SU55" & TGV60$Rep==rep]
      tgv$per_GeneticSD <- (tgv$SDSt / base)      

      
      tgv60 <- rbind(tgv60, tgv)

    }
  }
}
TGV60 <- tgv60


EFF <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      effBase <- regRep$Slope[regRep$Scenario=="Class"  & regRep$Strategy=="SU55" & regRep$Rep==rep]
      eff <- regRep[regRep$Scenario==scenario & regRep$Strategy==strategy & regRep$Rep==rep,]
      eff$per_Eff <- (eff$Slope / effBase)
      EFF <- rbind(EFF, eff)
    }
  }
}

EFF <- EFF[!(is.na(EFF$Slope)),]

TGV60$per_zMean <- TGV60$per_zMean*100 - 100 #tukaj ma osnoven scenarij 100: zato odštej 100
TGV60$per_GenicVar <- (TGV60$per_GenicVar)*100 - 100#tukaj ma osnoven scenarij 0
TGV60$per_GeneticVar <- (TGV60$per_GeneticVar)*100 - 100
TGV60$per_GenicSD <- (TGV60$per_GenicSD)*100 - 100#tukaj ma osnoven scenarij 0
TGV60$per_GeneticSD <- (TGV60$per_GeneticSD)*100 - 100
EFF$per_Eff <-  EFF$per_Eff * 100 -100

#genetic gain
MEAN60_abs <- summarySE(TGV60, measurevar="zMean", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
MEAN60 <- summarySE(TGV60, measurevar="per_zMean", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(MEAN60) <- c("Strategy", "Scenario", "per_zMean", "per_zMean_SD")
colnames(MEAN60_abs) <- c("Strategy", "Scenario", "zMean", "zMean_SD")
MEAN60$per_zMean <- round(MEAN60$per_zMean)
MEAN60$per_zMean_qL <- round(MEAN60$per_zMean_qL)
MEAN60$per_zMean_qH <- round(MEAN60$per_zMean_qH)
head(MEAN60)


#genetic variance
VAR60a <- summarySE(TGV60, measurevar="per_GeneticVar", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
VAR60a_abs <- summarySE(TGV60, measurevar="GeneticVarSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(VAR60a) <- c("Strategy", "Scenario", "per_GeneticVar", "per_GeneticVarSD")
colnames(VAR60a_abs) <- c("Strategy", "Scenario", "GeneticVar", "GeneticVarSD")
VAR60a$per_GeneticVar <- round(VAR60a$per_GeneticVar, 3)
VAR60a$per_GeneticVarSD <- round(VAR60a$per_GeneticVarSD, 3)

MEAN60 <- merge(MEAN60, VAR60a, by=c("Strategy", "Scenario"))
MEAN60_abs <- merge(MEAN60_abs, VAR60a_abs, by=c("Strategy", "Scenario"))

#genic variance
VAR60b <- summarySE(TGV60, measurevar="per_GenicVar", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
VAR60b_abs <- summarySE(TGV60, measurevar="GenicVarSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(VAR60b) <- c("Strategy", "Scenario", "per_GenicVar", "per_GenicVarSD")
colnames(VAR60b_abs) <- c("Strategy", "Scenario", "GenicVar", "GenicVarSD")
VAR60b$per_GenicVar <- round(VAR60b$per_GenicVar, 3)
VAR60b$per_GenicVarSD <- round(VAR60b$per_GenicVarSD, 3)
MEAN60 <- merge(MEAN60, VAR60b, by=c("Strategy", "Scenario"))
MEAN60_abs <- merge(MEAN60_abs, VAR60b_abs, by=c("Strategy", "Scenario"))

#genetic sd
SD60a <- summarySE(TGV60, measurevar="per_GeneticSD", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
SD60a_abs <- summarySE(TGV60, measurevar="SDSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(SD60a) <- c("Strategy", "Scenario", "per_GeneticSD", "per_GeneticSDSD")
colnames(SD60a_abs) <- c("Strategy", "Scenario", "GeneticSD", "GeneticSDSD")
SD60a$per_GeneticSD <- round(SD60a$per_GeneticSD, 3)
SD60a$per_GeneticSDSD <- round(SD60a$per_GeneticSDSD, 3)

MEAN60 <- merge(MEAN60, SD60a, by=c("Strategy", "Scenario"))
MEAN60_abs <- merge(MEAN60_abs, SD60a_abs, by=c("Strategy", "Scenario"))

#genic sd
SD60b <- summarySE(TGV60, measurevar="per_GenicSD", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
SD60b_abs <- summarySE(TGV60, measurevar="SDGenicSt", groupvars=c("Strategy", "scenario"))[,c(1,2,4,5)]
colnames(SD60b) <- c("Strategy", "Scenario", "per_GenicSD", "per_GenicSDSD")
colnames(SD60b_abs) <- c("Strategy", "Scenario", "GenicSD", "GenicSDSD")
SD60b$per_GenicSD <- round(SD60b$per_GenicSD, 3)
SD60b$per_GenicSDSD <- round(SD60b$per_GenicSDSD, 3)

MEAN60 <- merge(MEAN60, SD60b, by=c("Strategy", "Scenario"))
MEAN60_abs <- merge(MEAN60_abs, SD60b_abs, by=c("Strategy", "Scenario"))


#efficiency
eff <- summarySE(EFF, measurevar="per_Eff", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
eff_abs <- summarySE(EFF, measurevar="Slope", groupvars=c("Strategy", "Scenario"))[,c(1,2,4,5)]
colnames(eff) <- c("Strategy", "Scenario", "per_Eff", "per_EffSD")
colnames(eff_abs) <- c("Strategy", "Scenario", "Eff", "EffSD")
eff$per_Eff <- round(eff$per_Eff)
eff$per_EffSD <- round(eff$per_EffSD)

MEAN60 <- merge(MEAN60, eff, by=c("Strategy", "Scenario"))
MEAN60_abs <- merge(MEAN60_abs, eff_abs, by=c("Strategy", "Scenario"))


#for displaying
MEAN60$Strategy <- factor(MEAN60$Strategy, levels =c("SU55", "SU51", "SU15"))
MEAN60$Scenario <- factor(MEAN60$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
MEAN60 <- MEAN60[order(MEAN60$Strategy, MEAN60$Scenario),]

EFF$Strategy <- factor(EFF$Strategy, levels =c("SU55", "SU51", "SU15"))
EFF$Scenario <- factor(EFF$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
EFF <- EFF[order(EFF$Strategy, EFF$Scenario),]

MEAN60_abs$Strategy <- factor(MEAN60_abs$Strategy, levels =c("SU55", "SU51", "SU15"))
MEAN60_abs$Scenario <- factor(MEAN60_abs$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
MEAN60_abs[,3:ncol(MEAN60_abs)] <- round(MEAN60_abs[,3:ncol(MEAN60_abs)], 2)
MEAN60_abs[order(MEAN60_abs$Strategy, MEAN60_abs$Scenario),]


#write.csv(MEAN60, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//StandardDeviation_GeneticGain_gen60_16112018.csv", quote=FALSE)


#plot for genic / genetic variance
colnames(VAR60a) <- c("Strategy", "Scenario", "Value", "SD")
VAR60a$Type <- "Genetic"
VAR60a$measure <- "Var"
colnames(VAR60a_abs) <- c("Strategy", "Scenario", "Value", "SD")
VAR60a_abs$Type <- "Genetic"
VAR60a_abs$measure <- "Var"

colnames(VAR60b) <- c("Strategy", "Scenario", "Value", "SD")
VAR60b$Type <- "Genic"
VAR60b$measure <- "Var"

colnames(VAR60b_abs) <- c("Strategy", "Scenario", "Value", "SD")
VAR60b_abs$Type <- "Genic"
VAR60b_abs$measure <- "Var"

colnames(SD60a) <- c("Strategy", "Scenario", "Value", "SD")
SD60a$Type <- "Genetic"
SD60a$measure <- "SD"

colnames(SD60a_abs) <- c("Strategy", "Scenario", "Value", "SD")
SD60a_abs$Type <- "Genetic"
SD60a_abs$measure <- "SD"


colnames(SD60b) <- c("Strategy", "Scenario", "Value", "SD")
SD60b$Type <- "Genic"
SD60b$measure <- "SD"

colnames(SD60b_abs) <- c("Strategy", "Scenario", "Value", "SD")
SD60b_abs$Type <- "Genic"
SD60b_abs$measure <- "SD"


VARSD <- rbind(VAR60a, VAR60b)
VARSD <- rbind(VARSD, SD60a)
VARSD <- rbind(VARSD, SD60b)

VARSD_abs <- rbind(VAR60a_abs, VAR60b_abs)
VARSD_abs <- rbind(VARSD_abs, SD60a_abs)
VARSD_abs <- rbind(VARSD_abs, SD60b_abs)


VARSD$Scenario <- revalue(VARSD$Scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT"))
VARSD$Strategy1 <- revalue(VARSD$Strategy, c("SU55" = "5 sires/year, use 5 years","SU51" = "5 sires/year, use 1 year","SU15" = "1 sire/year, use 5 years"))
VARSD$Scenario <- factor(VARSD$Scenario, levels =c("PT", "GT-PT", "GT-C", "GT-BD", "GT"))
VARSD$Strategy1 <- factor(VARSD$Strategy1, levels =c("5 sires/year, use 5 years", "5 sires/year, use 1 year", "1 sire/year, use 5 years"))
VARSD$Type <- factor(VARSD$Type, levels =c("Genic", "Genetic"))
VARSD <- VARSD[order(VARSD$Strategy1, VARSD$Scenario),]

VARSD_abs$Scenario <- revalue(VARSD_abs$Scenario, c("Class" = "PT", "GenSLO" = "GT-PT", "OtherCowsGen" = "GT-C", "BmGen" = "GT-BD", "Gen" = "GT"))
VARSD_abs$Strategy1 <- revalue(VARSD_abs$Strategy, c("SU55" = "5 sires/year, use 5 years","SU51" = "5 sires/year, use 1 year","SU15" = "1 sire/year, use 5 years"))
VARSD_abs$Scenario <- factor(VARSD_abs$Scenario, levels =c("PT", "GT-PT", "GT-C", "GT-BD", "GT"))
VARSD_abs$Type <- factor(VARSD_abs$Type, levels =c("Genic", "Genetic"))


wrapit <- function(text) {
  wtext <- paste(strwrap(text,width=40),collapse=" \n ")
  return(wtext)
}


#barplot for genicSd
VARSD$wrapped_text <- llply(VARSD$Strategy1, wrapit)
VARSD$wrapped_text <- unlist(VARSD$wrapped_text)
VARSD_abs$wrapped_text <- llply(VARSD_abs$Strategy1, wrapit)
VARSD_abs$wrapped_text <- unlist(VARSD_abs$wrapped_text)
library(grid)
ggplot(data=VARSD[VARSD$measure=="SD",], aes(y=Value,x=Scenario,  fill=Type)) + geom_bar(position="dodge", stat="identity") +
  scale_fill_manual("Standard deviation", breaks = c("Genic", "Genetic"), 
                      values=c("darkgreen", "steelblue1")) + ylim(-15, 6) +
  xlab("Breeding program") + ylab("Percentage change to the baseline scenario") +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=3)) + theme_bw() + 
  theme(axis.text=element_text(size=14), legend.position = "top",
        axis.title=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14), 
        strip.text = element_text( size=16))   +
        theme(panel.spacing = unit(2, "lines")) + 
        geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), width=.2,
                position=position_dodge(.9)) +
       guides(fill=guide_legend(nrow=1)) +
       facet_grid(Strategy1 ~ ., scales = "free_y", labeller = label_wrap_gen(width=15)) 
#absolute
ggplot(data=VARSD_abs[VARSD_abs$measure=="SD",], aes(y=Value,x=Scenario,  fill=Type)) + geom_bar(position="dodge", stat="identity") +
  scale_fill_manual("Standard deviation", breaks = c("Genic", "Genetic"), 
                      values=c("darkgreen", "steelblue1")) + ylim(0, 1) +
  xlab("Breeding program") + ylab("Percentage change to the baseline scenario") +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=3)) + theme_bw() + 
  theme(axis.text=element_text(size=14), legend.position = "top",
        axis.title=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14), 
        strip.text = element_text( size=16))   +
        theme(panel.spacing = unit(2, "lines")) + 
        geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), width=.2,
                position=position_dodge(.9)) +
       guides(fill=guide_legend(nrow=1)) +
       facet_grid(Strategy1 ~ ., scales = "free_y", labeller = label_wrap_gen(width=15)) 
############################################################################################
############################################################################################
############################################################################################
#SIGNIFICANCE-s

#significance of efficiency
# library(emmeans)
# regRep$Scenario <- as.factor(regRep$Scenario)
# regRep$Strategy <- as.factor(regRep$Strategy)
# regRep1 <- within(regRep, Scenario <- relevel(Scenario, ref = "Class"))
# m1 <- lm(Slope~Scenario,data=regRep1[regRep1$Strategy=="SU55",])
# m1.grid <- ref_grid(m1)
# anova(m1)
# m1S <- lsmeans(m1.grid, "Scenario")
# contrast(m1.grid, method="pairwise")
# contrast(m1S, method="eff")
# summary(lm(Slope~Scenario,data=regRep1))

#significance of genetic gain
library(emmeans)
TGV60$scenario <- as.factor(TGV60$scenario)
TGV60$Strategy <- as.factor(TGV60$Strategy)

TGV60$Strategy <- factor(TGV60$Strategy, levels =c("SU55", "SU51", "SU15"))
TGV60$scenario <- factor(TGV60$scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))

#model <- lm(zMean ~ Strategy + scenario + Strategy : scenario, data=TGV60) to je model za ABSOLUTNE vrednsoti
model <- lm(per_zMean ~ Strategy + scenario + Strategy : scenario, data=TGV60)
model <- lm(zMean ~ Strategy + scenario + Strategy : scenario, data=TGV60) #absolute
library(car)
Anova(model,type = "II")
x = residuals(model)
library(rcompanion)
#plotNormalHistogram(x)
library(multcompView)
library(emmeans)
marginal = emmeans(model, ~ Strategy:scenario)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, by="Strategy",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#ALI tole - to so samo p vrednsoti
### same as TukeyHSD
TukeyHSD(aov(model), "Strategy:scenario")
TukeyHSD(aov(model), "Strategy")
TukeyHSD(aov(model), "scenario")
### set up linear hypotheses for all-pairs of both factors

#significance of genic variance
model <- lm(per_GenicVar ~ Strategy + scenario + Strategy : scenario, data=TGV60)
marginal = emmeans(model, ~ Strategy:scenario)
CLD = cld(marginal, by="Strategy",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05, sort=FALSE,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#significance of genetic variance
model <- lm(per_GeneticVar ~ Strategy + scenario + Strategy : scenario, data=TGV60)
marginal = emmeans(model, ~ Strategy:scenario)
CLD = cld(marginal, by="Strategy",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05, sort=FALSE,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#significance of genic SD
model <- lm(per_GenicSD ~ Strategy + scenario + Strategy : scenario, data=TGV60)
marginal = emmeans(model, ~ Strategy:scenario)
CLD = cld(marginal, by="Strategy",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05, sort=FALSE,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#significance of genetic SD
model <- lm(per_GeneticSD ~ Strategy + scenario + Strategy : scenario, data=TGV60)
marginal = emmeans(model, ~ Strategy:scenario)
CLD = cld(marginal, by="Strategy",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05, sort=FALSE,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#significance of efficiency - LETTERS
#regRep$Strategy <- factor(regRep$Strategy, levels =c("SU55", "SU51", "SU15"))
#regRep$Scenario <- factor(regRep$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
#regRep <- within(regRep, Scenario <- relevel(Scenario, ref = "Class"))
model <- lm(per_Eff ~ Strategy + Scenario + Strategy : Scenario, data=EFF)
model <- lm(Slope ~ Strategy + Scenario + Strategy : Scenario, data=EFF)
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



#m1 <- lm(zMean~Strategy,data=TGV601[TGV601$scenario=="Class",])
# m1.grid <- ref_grid(m1)
# anova(m1)
# m1S <- lsmeans(m1.grid, "scenario")
# m1S <- lsmeans(m1.grid, "Strategy")
# contrast(m1.grid, method="pairwise")
# contrast(m1S, method="eff")
# summary(lm(Slope~Scenario,data=regRep1))

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
                        labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  geom_line(data = MeanAverage[MeanAverage$Strategy=="SU55",], aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
  #geom_ribbon(data=MeanAverage, aes(x=Generation, ymin=lower, ymax=upper, colour=scenario), alpha=0.1) + 
ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=18), legend.title=element_text(face="bold", size=18), 
        strip.text = element_text(face="bold", size=16))   + 
facet_grid(order ~ ., scales = "free_y") + theme(legend.position = "right") 

#genetic variance plot
MeanAverageSD$order <- factor(MeanAverageSD$Strategy, levels = c("SU55", "SU15", "SU51"))

ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), 
                        labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
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
                        labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GT-PT", "GT-C", "GT-BD", "GT")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  geom_line(data = MeanAverageSDGenic, aes(x=Generation, y=SdGenic, colour=scenario, linetype=scenario), size=1.2) + 
  ylim(c(0.85, 1)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=16), legend.title=element_text(face="bold", size=18), 
        strip.text = element_text(face="bold", size=16))   + 
  facet_grid(order ~ ., scales = "free_y") + theme(legend.position = "right") 




#PLOT FOR MEASURES  - comparsion
library(reshape)
MEAN60$Group <- paste0(MEAN60$Strategy, MEAN60$Scenario)
MEAN60m <- melt(MEAN60[,3:ncol(MEAN60)], id = "Group")
MEAN60_abs$Group <- paste0(MEAN60_abs$Strategy, MEAN60_abs$Scenario)
MEAN60abs_m <- melt(MEAN60_abs[,3:ncol(MEAN60_abs)], id = c("Group"))

MEAN60abs_m$Group <- as.factor(MEAN60abs_m$Group)
#genic Sd and genetic gain
ggplot(data = MEAN60abs_m[MEAN60abs_m$variable %in% c("zMean", "GenicSD"),], aes(x = Group, y = value, group=variable, colour=variable )) + geom_line()
ggplot(data = MEAN60m[MEAN60m$variable %in% c("per_zMean", "per_GenicSD"),], aes(x = Group, y = value, group=variable, colour=variable )) + geom_line()

#add generation interval
#giPer
giPer_a$Group <- paste0(giPer_a$Strategy, giPer_a$Scenario)
giPerm <- melt(giPer_a[,3:6], id = c("Group", "Line"))
giPerM <- giPerm[giPerm$Line=="sireM",]
giPerF <- giPerm[giPerm$Line=="sireF",]

plotM <- rbind(MEAN60m, giPerM[,c(1,3,4)])
plotF <- rbind(MEAN60m, giPerF[,c(1,3,4)])

#generation interval and genetic gain
ggplot(data = plotM[plotM$variable %in% c("per_zMean", "per_gi"),], aes(x = Group, y = value, group=variable, colour=variable )) + geom_line()
ggplot(data = plotM[plotM$variable %in% c("per_GenicSD", "per_gi", "per_zMean"),], aes(x = Group, y = value, group=variable, colour=variable )) + geom_line()
ggplot(data = plotM[plotF$variable %in% c("per_zMean", "per_gi"),], aes(x = Group, y = value, group=variable, colour=variable )) + geom_line()


MEAN60_abs <- merge(MEAN60_abs, giPer_abs[giPer_abs$Line == "sireM",], by=c("Strategy", "Scenario"))
colnames(MEAN60_abs)[colnames(MEAN60_abs) == "Line"] <- "SireM"
colnames(MEAN60_abs)[colnames(MEAN60_abs) == "gi"] <- "gi_sireM"
colnames(MEAN60_abs)[colnames(MEAN60_abs) == "giSD"] <- "gi_sireMSD"
MEAN60_abs <- merge(MEAN60_abs, giPer_abs[giPer_abs$Line == "sireF",], by=c("Strategy", "Scenario"))
colnames(MEAN60_abs)[colnames(MEAN60_abs) == "Line"] <- "SireF"
colnames(MEAN60_abs)[colnames(MEAN60_abs) == "gi"] <- "gi_sireF"
colnames(MEAN60_abs)[colnames(MEAN60_abs) == "giSD"] <- "gi_sireFSD"


model <- glm(MEAN60_abs$GenicSD ~ MEAN60_abs$gi_sireM + MEAN60_abs$Strategy)
summary(model)
