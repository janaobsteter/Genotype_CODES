library(ggplot2)
p = ggplot(data=TGVsAll, mapping=aes(x=zSdGenic, y=zMean, color=scenario, linetype=scenario)) +  geom_path(size=1/5, alpha=.5) + 
#  # Conv, PYT, Head, TwoPartTS, TwoPartTS+, TwoPartOCS  
scale_color_manual(values=c("black", "red", "green", "blue", "yellow"), name="") +  
scale_linetype_manual(values=c(      1,    1,    1,    1,    2,    1), name="") +  
xlab("Genic standard deviation") + ylab("Genetic mean") +  
scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                    name="Converted/Lost genic standard deviation")) +  
theme_bw() +  theme(legend.position="top") +  
guides(colour=guide_legend(nrow=1)) +  
geom_segment(data=tmp5, mapping=aes(x=zSdGenicMax, xend=zSdGenicMin,                                      
                                     y=zMeanStartSdGenic, yend=zMeanEndSdGenic,                                      
                                     color=scenario3, linetype=scenario3)) +  
geom_segment(data=tmp5, mapping=aes(x=zSdGenicMax, xend=zSdGenicMin,                                      
                                      y=zMeanStartSdGenic, yend=zMeanEndSdGenic,                                      
                                      color=scenario3, linetype=scenario3),               
               arrow=arrow(), show.legend=FALSE)



efficiency = function(x, f) {
  ret = NA
  if (nrow(x) > 5) {
    x = arrange_(x, "Generation")
    # fit = lm(formula=f, data=x)
    fit = MASS:::rlm(formula=f, data=x, maxit=1000)
    # fit = MASS:::lqs(formula=f, data=x)
    ret = -coef(fit)[2]
  }
  ret
}
datByRunStage = dat %>%
  group_by(run, rep, scenario, sel, self, size, scaled, crit, nCycles, target, stage) %>%
  nest()
data <- TGVsAll
datByRunStage = datByRunStage %>%
  mutate(efficiencyG=map(data, efficiency, f=zMean      ~ zSdGenic),
         efficiencyGenic=map(data, efficiency, f=zMeanGenic ~ zSdGenic))
# zMean je standardiziran genetski napredek (TBV - mean(TBV_gen_start)) / sd(TBV_gen_start)
# zMeanGenic je standardiziran genetski napredek (TBV - mean(TBV_gen_start)) / sqrt(2*sum(p*q*alpha)_gen_start


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
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################

TGVsAll <- data.frame()

WorkingDir = '/home/jana/Simulation_Rep0'
for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
  TGVs <- data.frame(Generation=40:60)
  ped <- read.table(paste0(WorkingDir,'/Pedigree', scenario, '.txt'), header=T)
  ped <- ped[ped$Generation %in% 40:60,]
  TGV <- summarySE(ped, measurevar = "gvNormUnres1", groupvars = "Generation")[,c(1,3,4)]
  colnames(TGV)[1] <- c("Generation")
  TGV$zMean <- (TGV$gvNormUnres1 - TGV$gvNormUnres1[1]) / TGV$sd[1]
  TGVs <- merge(TGVs, TGV, by="Generation")
  Var <- read.table(paste0(WorkingDir,'/Variance', scenario, '.txt'), header=T)
  Var <- Var[Var$QtnModel==1,c(1,3)]
  TGVs <- merge(TGVs, Var, by="Generation")
  TGVs$zSdGenic <- (sqrt(TGVs$AdditGenicVar1)) 
  TGVs$SDGenicSt <- TGVs$AdditGenicVar1 / TGVs$AdditGenicVar1[1] 
  TGVs$SDSt <- TGVs$sd / TGVs$sd[1] 
  TGVs$zMeanGenic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$AdditGenicVar1[1]
  #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
  TGVs$scenario <- scenario
  #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
  TGVsAll <- rbind(TGVsAll, TGVs)
}

#40 - 60 generacije
TGVsAll <- read.table("~/Documents/PhD/Simulaton/RefPopSize/TGVsAll_10KRef_20Gen.csv", header=TRUE) # to je stara OtherCowsGen - napačen GI
TGVsAll <- read.table("~/Documents/PhD/Simulaton/RefPopSize/TGVsAll_10KRef_new.csv", header=TRUE) # to je stara OtherCowsGen - napačen GI
TGVsAll <- read.table("TGVsAll.csv", header=TRUE)
TGVsAll <- read.table("~/Documents/PhD/Simulaton/RefPopSize/TGVsAll_1KRef_20Gen.csv", header=TRUE)
TGVsAll <- read.table("~/Documents/PhD/Simulaton/RefPopSize/TGVsAll_10KRef_30Rep.csv", header=TRUE)
TGVsAll <- read.table("~/Documents/PhD/Simulaton//TGVsAll_10KRef_2501.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/Genotipi/Genotipi_CODES//TGVsAll_10KRef_2501.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/TGVsAll_10KRef_2801_1Pb.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/OCS/TGVsAll_10KRef_OCS.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/Documents/PhD/Simulaton//TGVsAll_10KRef_2801_1Pb.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/Documents/PhD/Simulaton//TGVsAll_10KRef_2801_1Year.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/Genotipi/Genotipi_CODES//TGVsAll_10KRef_2801_1Year.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
#1 - 60 generacije
TGVsAll <- read.table("~/Documents/PhD/Simulaton/RefPopSize/TGVsAll_10KRef_60Gen.csv", header=TRUE)
TGVsAll$Group <- paste(TGVsAll$scenario, TGVsAll$Rep)
#TGVsAll <- TGVsAll[TGVsAll$Rep %in% c(0,1,2,3,5,6,7,8,9,10),]
TGVsAll <- TGVsAll[TGVsAll$Generation %in% 40:60,]


#Dodaj še nove spremenljivke
#in spravi gensko in genetsko varianco na isto točko
groupDF <- data.frame()
for (group in unique(TGVsAll$Group)) {
  groupT <- subset(TGVsAll[TGVsAll$Group==group,])
  groupT$SDGenicStNeg <- 1- groupT$AdditGenicVar1 / groupT$AdditGenicVar1[1]
  groupDF <- rbind(groupDF, groupT)
}

#standadise onto the 40 generation
st40 <- data.frame()
for (rep in 0:20) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    repScenario <- TGVsAll[(TGVsAll$Rep==rep) & (TGVsAll$scenario==scenario) & (TGVsAll$Generation %in% 40:60),]
    repScenario$zMean <- (repScenario$gvNormUnres1 - repScenario$gvNormUnres1[1]) / repScenario$sd[1]
    repScenario$SDGenicSt <- repScenario$AdditGenicVar1 / repScenario$AdditGenicVar1[1] 
    repScenario$SDSt <- repScenario$sd / repScenario$sd[1] 
    repScenario$zMeanGenic <- (repScenario$gvNormUnres1 - repScenario$gvNormUnres1[1]) / repScenario$AdditGenicVar1[1]
    repScenario$SDGenicStNeg <- 1- repScenario$AdditGenicVar1 / repScenario$AdditGenicVar1[1]
    #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
    #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
    st40 <- rbind(st40, repScenario)
  }
}

TGVsAll <- st40
varDF <- data.frame()
for (group in unique(TGVsAll$Group)) {
  groupT <- subset(TGVsAll[TGVsAll$Group==group,])
  groupT$var <- (groupT$sd)^2
  koef <- groupT$var[1] / groupT$AdditGenicVar1[1]
  groupT$GenicVAR_genetic <- groupT$AdditGenicVar1 * koef
  groupT$GenicSD_genetic <- sqrt(groupT$GenicVAR_genetic)
  groupT$zMeanGenic <- (groupT$gvNormUnres1 - groupT$gvNormUnres1[1]) / groupT$zSdGenic[1]
  groupT$SDGenicSt <- groupT$zSdGenic / groupT$zSdGenic[1] 
  groupT$MeanGENIC_genetic <- (groupT$gvNormUnres1 - groupT$gvNormUnres1[1]) / groupT$GenicSD_genetic[1]
  groupT$GENICSd_St <- groupT$GenicSD_genetic / groupT$GenicSD_genetic[1] 
  
  varDF <- rbind(varDF, groupT)
}
TGVsAll <- varDF
'''
library(MASS)
library(ggplot2)
lm <- ggplot(data = TGVsAll, aes(x=TGVsAll$zSdGenic, y=TGVsAll$zMeanGenic, colour=TGVsAll$scenario)) + geom_path() + scale_x_reverse() + 
  geom_smooth(method='lm', se=FALSE) + 
  scale_color_hue("Shema", labels=c("Conventional", "GenomicSLO", "GenBulls on Other Cows", "GenBulls on Bull Dams", "GenBulls on All Cows")) + 
  xlab("Genic sd") + ylab("Mean genetic gain") + ggtitle("Genic variance standardisation")
lm2 <- ggplot(data = TGVsAll, aes(x=TGVsAll$zSdGenic, y=TGVsAll$zMean, colour=TGVsAll$scenario)) + geom_path() + scale_x_reverse() + 
  geom_smooth(method='lm', se=FALSE) + 
  scale_color_hue("Shema", labels=c("Conventional", "GenomicSLO", "GenBulls on Other Cows", "GenBulls on Bull Dams", "GenBulls on All Cows")) + 
  xlab("Genic sd") + ylab("Mean genetic gain") + ggtitle("Genetic variance standardisation")
rlm <- ggplot(data = TGVsAll, aes(x=TGVsAll$zSdGenic, y=TGVsAll$zMeanGenic, colour=TGVsAll$scenario)) + geom_path() + scale_x_reverse() + geom_smooth(method='rlm', se=FALSE)
library(Rmisc)
multiplot(lm, rlm, cols=2)
multiplot(lm, lm2, cols=2)
'''
#TGVsAll <- TGVsAll[TGVsAll$Rep %in% 10:20,]


#TGVsAll <- TGVsAll[TGVsAll$Rep %in% 10:20,]
#naredi povprečja replik
Averages1 <- aggregate(TGVsAll$GENICSd_St ~TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages2 <- aggregate(TGVsAll$zMean ~TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages3 <- aggregate(TGVsAll$SDSt ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages4 <- aggregate(TGVsAll$zMean ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages5 <- aggregate(TGVsAll$zSdGenic ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages6 <- aggregate(TGVsAll$sd ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages7 <- aggregate(TGVsAll$MeanGENIC_genetic ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
Averages8 <- aggregate(TGVsAll$GENICSd_St ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
colnames(Averages1) <- c("scenario", "Generation", "SDGenic")
colnames(Averages2) <- c("scenario", "Generation", "MeanGenic")
colnames(Averages3) <- c("scenario", "Generation", "SDGenetic")
colnames(Averages4) <- c("scenario", "Generation", "MeanGenetic")
colnames(Averages5) <- c("scenario", "Generation", "SDGenic_noSt")
colnames(Averages6) <- c("scenario", "Generation", "SDGenetic_noSt")
colnames(Averages7) <- c("scenario", "Generation", "MeanGenicGenetic")
colnames(Averages8) <- c("scenario", "Generation", "SDGenic_Genetic")
Averages <- merge(Averages1, Averages2, by=c( "scenario", "Generation"))
Averages <- merge(Averages, Averages3, by=c( "scenario", "Generation"))
Averages <- merge(Averages, Averages4, by=c( "scenario", "Generation"))
Averages <- merge(Averages, Averages5, by=c( "scenario", "Generation"))
Averages <- merge(Averages, Averages6, by=c( "scenario", "Generation"))
Averages <- merge(Averages, Averages7, by=c( "scenario", "Generation"))
Averages <- merge(Averages, Averages8, by=c( "scenario", "Generation"))

#AveragesA <- Averages
#to je max min za gensko varianco
maxmin <- data.frame(scenario=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (scenario in unique(Averages$scenario)) {
    maxmin[row,] <- c(scenario, min(Averages$SDGenic[Averages$scenario==scenario]), max(Averages$SDGenic[Averages$scenario==scenario]), min(Averages$MeanGenic[Averages$scenario==scenario]), max(Averages$MeanGenic[Averages$scenario==scenario]))
    row <- row +1
}
#to je max min za gensko varianco standardizirano na genetsko varianco
maxmin <- data.frame(scenario=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (scenario in unique(Averages$scenario)) {
    maxmin[row,] <- c(scenario, min(Averages$SDGenic_Genetic[Averages$scenario==scenario]), max(Averages$SDGenic_Genetic[Averages$scenario==scenario]), 
                      min(Averages$MeanGenicGenetic[Averages$scenario==scenario]), max(Averages$MeanGenicGenetic[Averages$scenario==scenario]))
    row <- row +1
}


#to je max-min po ponovitvah, zgoraj je povprečje
maxmin <- data.frame(scenario=NA, rep = NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (rep in unique(TGVsAll$Rep)) {
  for (scenario in unique(Averages$scenario)) {
    maxmin[row,] <- c(scenario, rep, min(Averages$SDGenic_Genetic[Averages$scenario==scenario]), max(Averages$SDGenic_Genetic[Averages$scenario==scenario]), 
                      min(Averages$MeanGenicGenetic[Averages$scenario==scenario]), max(Averages$MeanGenicGenetic[Averages$scenario==scenario]))
    row <- row +1
  }
}


#to je max min za genetsko varianco
maxmin <- data.frame(scenario=NA, minGeneticSD=NA, maxGeneticSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (scenario in unique(Averages$scenario)) {
    maxmin[row,] <- c(scenario, min(Averages$SDGenetic[Averages$scenario==scenario]), max(Averages$SDGenetic[Averages$scenario==scenario]), min(Averages$MeanGenetic[Averages$scenario==scenario]), max(Averages$MeanGenetic[Averages$scenario==scenario]))
    row <- row +1
}




#genska
maxmin$minGenicSD <- as.numeric(maxmin$minGenicSD)
maxmin$maxGenicSD <- as.numeric(maxmin$maxGenicSD)
maxmin$minTGV <- as.numeric(maxmin$minTGV)
maxmin$maxTGV <- as.numeric(maxmin$maxTGV)
#genetksa
maxmin$minGeneticSD <- as.numeric(maxmin$minGeneticSD)
maxmin$maxGeneticSD <- as.numeric(maxmin$maxGeneticSD)
maxmin$minTGV <- as.numeric(maxmin$minTGV)
maxmin$maxTGV <- as.numeric(maxmin$maxTGV)


TGVsAllOrig <- TGVsAll
TGVsAllOrig$design <- "bull55"
AveragesOrig <- Averages
AveragesOrig$design <- "bull55"
TGVsAll1Pb <- TGVsAll
TGVsAll1Pb$design <- "bull15"
Averages1Pb <- Averages
Averages1Pb$design <- "bull15"
TGVsAll1Year <- TGVsAll
TGVsAll1Year$design <- "bull51"
Averages1Year <- Averages
Averages1Year$design <- "bull51"

TGVSALL <- rbind(TGVsAllOrig, TGVsAll1Pb, TGVsAll1Year)
#tukaj spravi gensko varianco in genetsko na isto točko

library(ggplot2)
#write.table(TGVsAll, "~/Simulation_Rep0/TGVsAll.csv", quote=FALSE, row.names=FALSE)
#TGVsAll <- read.table("~/Documents/WCGALP/TGVsAll.csv", header=TRUE)
#To je plot zMean (standardizirana na gensko variacno) na genetsko varianco

TGVsAll$Sim <- ifelse(TGVsAll$Rep %in% 0:10, "1", "2")
TGVsAll$Sim <- as.factor(TGVsAll$Sim)


TGVsAll1 <- TGVsAll[TGVsAll$Sim == 1,]
TGVsAll2 <- TGVsAll[TGVsAll$Sim == 2,]

ggplot(data = TGVsAll1, aes(x=SDGenicSt, y=zMeanGenic, group=Group, colour=Sim, linetype=scenario)) + 
  geom_line(aes(linetype=TGVsAll1$scenario), size=0.5, alpha=0.3) + scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                                                                                     name="Converted/Lost genic standard deviation")) +
  #geom_smooth( se=FALSE, formula=y~x+1, method="lm") + 
  xlab("Generation") + ylab("True genetic value")  + 
  coord_cartesian(xlim = c(1, 0.8)) +

  scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        "Scenario", 
                        values=c("solid", "dotted","dashed", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Genic standard deviation") + ylab("Average True Genetic Value") + guides(linetype = guide_legend(override.aes = list(size=10))) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=FALSE, size=1.5)

#o je z na genetsko standadrizirano gensko variacno
Year1Eff <- ggplot(data = TGVsAll, aes(x=GENICSd_St, y=MeanGENIC_genetic, group=Degree, colour=Degree, linetype=scenario)) + 
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                                                                                     name="Converted/Lost genic standard deviation")) +
  geom_line(aes(linetype=degree), size=0.5, alpha=0.4) +
  xlab("Generation") + ylab("True genetic value")  + ylim(0,7) +coord_cartesian(xlim = c(1, 0.85)) +
  scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        "Scenario", 
                        values=c("solid", "dotted","dashed", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Genic standard deviation") + ylab("Average True Genetic Value") + 
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=18)) +
  geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=FALSE, size=1.5)




ggplot(data = TGVsAll[TGVsAll$Group=="Class 0",], aes(x=GENICSd_St, y=MeanGENIC_genetic)) + 
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                    name="Pretvorjena/izgubljena genska standardna deviacija")) +
  #geom_smooth( se=FALSE, formula=y~x+1, method="lm") + 
  xlab("Generation") + ylab("True genetic value") +  coord_cartesian(xlim = c(1, 0.95)) +
  scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        "Scenario", 
                        values=c("solid", "dotted","dashed", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("black", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Genska standardna deviacija") + ylab("Plemenska vrednost") + 
  theme(axis.text=element_text(size=16), legend.position = "left",
        axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title=element_text(size=18)) +
  geom_segment(data=maxmin[maxmin$scenario=="Class",], mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=FALSE, size=1.5)

e#To je plot zMean (standardizirana na GENETSKO variacno) 
ggplot(data = TGVsAll, aes(x=SDSt, y=zMean, group=Group, colour=scenario, linetype=scenario)) + 
  geom_line(aes(linetype=TGVsAll$scenario), size=0.5, alpha=0.3) + scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                                                                                     name="Converted/Lost genic standard deviation")) +
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        labels= c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"), "Scenario", 
                        values=c("solid", "dotted","dashed", "dotdash", "twodash")) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      labels= c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"), "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1")) + 
  xlab("Genic standard deviation") + ylab("Average True Genetic Value") +  #+
  geom_segment(data=maxmin, mapping=aes(x=maxGeneticSD, xend=minGeneticSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=FALSE, size=1.5)
#To je plot GENETSKE variance po generaijch po scenariih
genetic <- ggplot(data = TGVsAll, aes(x=Generation, y=AdditGenicVar1, group=Group, colour=scenario, linetype=scenario)) +  geom_line(aes(linetype=scenario), size=1) + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Additive Genetic Variance")  #+
  geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minGenicSD,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario),      arrow=arrow(), show.legend=FALSE)
#To je plot GENSKE variance po generaijch po scenariih
genic <- ggplot(data = TGVsAll, aes(x=Generation, y=zSdGenic, group=Group, colour=scenario, linetype=scenario)) +  geom_line(aes(linetype=scenario), size=1) + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Additive Genic Variance")  #+
  geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minGenicSD,  yend=maxTGV,                                    
                                        color=scenario, linetype=scenario, group=scenario),      arrow=arrow(), show.legend=FALSE)

MeanAverage <- aggregate(TGVsAll$zMean ~ TGVsAll$scenario + TGVsAll$Generation, FUN="mean")
MeanAverageSD <- aggregate(TGVsAll$zMean ~ TGVsAll$scenario + TGVsAll$Generation, FUN="sd")
colnames(MeanAverage) <- c("scenario", "Generation", "MeanTGV")
#genetic gain
mean <- lmList(MeanTGV ~ Generation | scenario, data=MeanAverage)
meanDF <- data.frame(Scenario=rownames(coef(mean)),coef(mean),check.names=FALSE)
colnames(meanDF) <- c("Scenario", "Mean_Intercept", "Mean_Slope")
meanDF$test <- "Reference10000"
#efficiency
#to je regresija na povprečja
avg <- lmList(MeanGenic ~ SDGenic | scenario, data=Averages)
Avg <- data.frame(Rep=rownames(coef(avg)),coef(avg),check.names=FALSE)
colnames(Avg) <- c("ScTGVsAllenario", "Intercept", "Slope")
Avg$method <- "Regression On Averages"

#to je regresija na vse
tot <- lmList(zMeanGenic ~ SDGenicSt | scenario, data=TGVsAll)
Tot <- data.frame(Rep=rownames(coef(tot)),coef(tot),check.names=FALSE)
colnames(Tot) <- c("Scenario", "Intercept", "Slope")
Tot$method <- "Regression on All Data"

#to je povrepčje regresij
library(nlme)
regRep <- data.frame(Rep=NA, Intercept=NA, Slope=NA, Scenario=NA)
for (sc in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
  df <- TGVsAll[TGVsAll$scenario==sc,]
  fm1 <- lmList(zMeanGenic ~ SDGenicSt | Rep, data=df)
  tmp <- data.frame(Rep=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
  colnames(tmp) <- c("Rep", "Intercept", "Slope")
  tmp$Scenario <- sc
  regRep <- rbind(regRep, tmp)
}

#povprečje regresijskih koefificentov
avgSlope <- aggregate(regRep$Slope ~ regRep$Scenario, FUN=mean)
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario, FUN=mean)
avgReg <- cbind(avgInt, avgSlope[,2])
colnames(avgReg) <- c("Scenario", "Intercept", "Slope")
avgReg$method <- "Average regression"
##################################################3
#TO JE ŠE z 1 - sd --> da dobiš intercept 0
#to je povrepčje regresij
library(nlme)
regRep <- data.frame(Rep=NA, Intercept=NA, Slope=NA, Scenario=NA)
for (sc in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
  df <- varDF[varDF$scenario==sc,]
  fm1 <- lmList(zMeanGenic ~ SDGenicStNeg | Rep, data=df)
  tmp <- data.frame(Rep=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
  colnames(tmp) <- c("Rep", "Intercept", "Slope")
  tmp$Scenario <- sc
  regRep <- rbind(regRep, tmp)
}
regRep <- regRep[-1,]

#preveri, ali je klasična najslabš v vseh
for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
  vec <- c()
  for (rep in 0:20) {
    repDF <- regRep[regRep$Rep==rep,]
    a <- repDF[order(repDF$Slope),]
    #vec <- c(vec, repDF$Scenario[repDF$Slope == min(repDF$Slope)] == "Class")
    vec <- c(vec, which(a$Scenario==scenario))
  }
  print(scenario)
  print(sum(vec == 5))
}
#povprečje regresijskih koefificentov
avgSlope <- aggregate(regRep$Slope ~ regRep$Scenario, FUN="mean")
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario, FUN="mean")
avgReg <- cbind(avgInt, avgSlope[,2])
colnames(avgReg) <- c("Scenario", "Eff_Intercept", "Eff_Slope")
avgReg$method <- "Average regression"
avgReg$test <- "Reference10000"
#########################################################
#združi tabelo z povporečji in z učinkovitostjo
testDF <- merge(meanDF, avgReg, by=c("Scenario", "test"))
write.csv(testDF, "~/Documents/PhD/Simulaton/RefPopSize/MeanAndEff_10000Ref.csv", quote=FALSE, row.names=FALSE)
#####################
#združi za različne teste
ref1K <- read.csv("~/Documents/PhD/Simulaton/RefPopSize/MeanAndEff_1000Ref.csv")
ref1K$meanDiff <- round((ref1K$Mean_Slope / min(ref1K$Mean_Slope) - 1) * 100, 1)
ref1K$effDiff <- round((ref1K$Eff_Slope / min(ref1K$Eff_Slope) - 1) * 100, 1)
target <- c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")
ref1K <- ref1K[match(target, ref1K$Scenario),]
ref1K$test <- "1K"

ref10K <- read.csv("~/Documents/PhD/Simulaton/RefPopSize/MeanAndEff_10000Ref.csv")
ref10K$meanDiff <- round((ref10K$Mean_Slope / min(ref10K$Mean_Slope) - 1) * 100, 1)
ref10K$effDiff <- round((ref10K$Eff_Slope / min(ref10K$Eff_Slope) - 1) * 100, 1)
ref10K <- ref10K[match(target, ref10K$Scenario),]
ref10K$test <- "10K"

ref <- rbind(ref1K, ref10K)
write.table(ref, "~/Documents/PhD/Simulaton/RefPopSize/Compare_1k_10K_RefSize.csv", quote=FALSE, row.names=FALSE)
####################################################
#zračunaj značilnost razlik naklonov med scenariji
library(lsmeans)
regRep$Scenario <- as.factor(regRep$Scenario)
regRep1 <- within(regRep, Scenario <- relevel(Scenario, ref = "Class"))
m1 <- lm(Slope~Scenario,data=regRep)
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "Scenario")
contrast(m1.grid, method="pairwise")
contrast(m1S, method="eff")
summary(lm(Slope~Scenario,data=regRep1))

AllMethods <- rbind(Avg, Tot, avgReg)
AllMethodsO <- AllMethods[order(AllMethods$Scenario),]
write.csv(AllMethodsO, "~/Documents/WCGALP/Efficiencies_DifferentMethods.csv", quote=FALSE, row.names=FALSE)

AllM <- melt(AllMethods, measure.vars = "Slope")


#plot napredka po replikaj - GENETIC GAIN
TGVsAll <- TGVsAll[TGVsAll$scenario %in% c("Class", "Gen"),]
MeanAverage <- MeanAverage[MeanAverage$scenario %in% c("Class", "Gen"),]
genGainPlot <- ggplot() + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
  ggtitle("a") + ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=20), legend.position = "left",
                                         axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20)) #legend.position="bottom", 

MeanAverage$Generation1 <- MeanAverage$Generation - 40
ggplot() + 
  xlab("Generacija") + ylab("Plemenska vrednost")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Klasična", "Genomic A", "Genomic B", "Genomic C", "Genomska")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("red3", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("Klasična", "Genomic A", "Genomic B", "Genomic C", "Genomska")) + 
  geom_line(data = MeanAverage, aes(x=Generation1, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=20), legend.position = "left",
                                         axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_blank()) #legend.position="bottom", 


Sim1 <- ggplot() + geom_line(data = TGVsAll[TGVsAll$Sim==1,], aes(x=Generation, y=zMean, group=Group, colour=scenario, linetype=scenario), size=0.5, alpha=0.3) + # geom_line(aes(linetype=scenario), size=0.5, alpha=0.4) + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  #geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
  ggtitle("a") + ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=14), legend.position = "left",
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14)) #legend.position="bottom", 
Sim2 <- ggplot() + geom_line(data = TGVsAll[TGVsAll$Sim==2,], aes(x=Generation, y=zMean, group=Group, colour=scenario, linetype=scenario), size=0.5, alpha=0.3) + # geom_line(aes(linetype=scenario), size=0.5, alpha=0.4) + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  #geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
  ggtitle("a") + ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=14), legend.position = "left",
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14)) #legend.position="bottom", 

ggplot() + geom_line(data = TGVsAll, aes(x=Generation, y=zMean, group=Group, colour=scenario, linetype=scenario), size=0.5, alpha=0.3) + # geom_line(aes(linetype=scenario), size=0.5, alpha=0.4) + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Average True Genetic Value") +
  #geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
  ggtitle("a") + ylim(c(0, 7)) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=14), legend.position = "left",
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14)) #legend.position="bottom", 

library(Rmisc)
multiplot(Sim1, Sim2)
##########################
#Tukaj oceni efektivno velikost populacije - z gensko varianco

TGVs1 <- TGVsAll
'''
TGVsAll <- TGVsAll[TGVsAll$Generation %in% 30:40,]
Variances <- read.table("~/GenicVARIANCE.txt", header=TRUE) #to je samo zato, dasi preverila, če dela ok
Variances <- Variances[Variances$QtnModel==1,]
colnames(Variances)[9:10] <- c("Rep", "Scenario")
Variances <- Variances[Variances$Generation != "Generation",]
Variances$Group <- paste0(Variances$Scenario, Variances$Rep)
Variances$AdditGenicVar1 <- as.numeric(as.character(Variances$AdditGenicVar1))
Variances$sdGenic <- sqrt(Variances$AdditGenicVar1)
Variances$Generation <- as.numeric(as.character(Variances$Generation))

VAR <- Variances
Variances <- Variances[Variances$Generation %in% 20:35,]
Nes <- data.frame()

for (group in unique(TGVsAll$Group)) {
  df <- subset(TGVsAll[TGVsAll$Group==group,])
  fit = glm(df$zSdGenic ~ df$Generation, family = Gamma(link = "log"))
  tmp <- data.frame(Scenario=unique(df$scenario), rep = unique(df$Rep), Intercept=coef(fit)[1], Slope=coef(fit)[2], check.names=FALSE)
  Nes <- rbind(Nes, tmp)
}
'''

#po intervalih
Nes <- data.frame()
for (group in unique(TGVsAll$Group)) {
  for (int in c(0,20,40)) {
    TGVsAll1 <- TGVsAll
    df <- TGVsAll[(TGVsAll$Group == group) & (TGVsAll$Generation %in% int:(int+20)),]
    fit = glm(df$AdditGenicVar1 ~ df$Generation, family = Gamma(link = "log"))
    tmp <- data.frame(Scenario=unique(df$scenario), Interval = paste0(int, ":", int+20), rep = unique(df$Rep), method="GenicVar", Intercept=coef(fit)[1], Slope=coef(fit)[2], check.names=FALSE)
    Nes <- rbind(Nes, tmp)
    }
}

Nes$dC = 1 - exp(Nes$Slope)
Nes$Ne = 1 / (2 * Nes$dC)

NEs <- aggregate(Nes$Ne ~Nes$Scenario + Nes$Interval + Nes$method, FUN="mean")
colnames(NEs) <- c("Scenario", "Interval[gen]", "Method", "Ne")

###############
#for Variances
for (group in unique(Variances$Group)) {
  df <- subset(Variances[Variances$Group==group,])
  fit = glm(df$sdGenic ~ df$Generation, family = Gamma(link = "log"))
  tmp <- data.frame(Scenario=unique(df$Scenario), rep = unique(df$Rep), Intercept=coef(fit)[1], Slope=coef(fit)[2], check.names=FALSE)
  Nes <- rbind(Nes, tmp)
}

Nes$dC = 1 - exp(Nes$Slope)
Nes$Ne = 1 / (2 * Nes$dC)

NEs <- aggregate(Nes$Ne ~Nes$Scenario, FUN="mean")
TGVsAll <- TGVs1
Variances <- VAR

#############################################3
#Tukaj oceni efektivno velikost populacije - z GENETSKO varianco
Nes1 <- data.frame()
for (group in unique(TGVsAll$Group)) {
  for (int in c(0,20,40)) {
    TGVsAll1 <- TGVsAll
    df <- TGVsAll[(TGVsAll$Group == group) & (TGVsAll$Generation %in% int:(int+20)),]
    fit = glm((df$sd)^2 ~ df$Generation, family = Gamma(link = "log"))
    tmp <- data.frame(Scenario=unique(df$scenario), Interval = paste0(int, ":", int+20), rep = unique(df$Rep), method="GeneticVar", Intercept=coef(fit)[1], Slope=coef(fit)[2], check.names=FALSE)
    Nes1 <- rbind(Nes1, tmp)
  }
}

Nes1$dC = 1 - exp(Nes1$Slope)
Nes1$Ne = 1 / (2 * Nes1$dC)

NEs1 <- aggregate(Nes1$Ne ~Nes1$Scenario + Nes1$Interval + Nes1$method, FUN="mean")
colnames(NEs1) <- c("Scenario", "Interval[gen]", "Method", "Ne")

NES <- rbind(NEs, NEs1)

#what is happening with CLass genetic variance??????
Class0 <- TGVsAll[TGVsAll$Group=="Class 0",]
ggplot(data=Class0, aes(x=Generation, y=sd)) + geom_path()
####################################################3
#tukaj preveri značilnost razlik med TGV v zadnji generaciji
library(emmeans)
gen60 <- TGVSALL[TGVSALL$Generation==60,]
m2 <- lm(zMean~scenario + design + scenario*design,data=gen60)
m2.grid <- ref_grid(m2)
anova(m2)
m2S <- lsmeans(m2.grid, c("scenario","design"))
contrast(m2.grid, method="eff")
contrast(m2S, method="pairwise")
summary(lm(zMean~scenario,data=gen60))

write.csv(TGVSALL, "TGVSALL.csv", quote=FALSE, row.names=FALSE)
#########################################
#Generacijski interval
GI <- read.table("~/GI.txt", header=TRUE) #to je sam do 59, ker za 60 ne moreš še imet GI
GI <- GI[GI$Gen != "Gen",]
GI$path <- paste0(GI$line, GI$sex)
GI$genInt <- as.numeric(as.character(GI$genInt))
colnames(GI) <- c("Gen", "genInt", "line", "sex", "Group", "path")
colnames(GI) <- c("Gen", "genInt", "line", "sex", "scenario", "path")
GI$Group <- as.character(GI$Group)
GI$scenario <- strsplit(GI$Group, "_")[[1]][2]
GI$Rep <- strsplit(GI$Group, "_")[[1]][1]

GIa <- aggregate(GI$genInt ~ GI$Gen + GI$path + GI$scenario, FUN="mean")
names <- data.frame(Orig = c("Class", "GenSLO",  "OtherCowsGenGen","BmGen", "Gen"), ID = c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"))
colnames(GIa) <- c("Gen", "Path", "Scenario", "genInt")
ScOrder <- c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")
require(gdata)
GIa$Scenario <- reorder.factor(GIa$Scenario, new.order=ScOrder)
require(dplyr)
GIa <- GIa %>%  arrange(Scenario)


'''
plotList <- list()
number = 1
for (scenario in unique(GIa$Class)) {
  GIplot <- GIa[GIa$Class == scenario,]
  plotList [[number]] <- ggplot(GIplot, aes(x=Gen, y=genInt, colour=Path, linetype =Path, group=Path)) + geom_line(aes(linetype=Path), size=1) + 
    scale_linetype_manual("", values=c("solid", "dashed", "dotted", "dotdash"), labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + 
    scale_colour_manual("", values=c("palevioletred2", "purple2", "royalblue3", "darkolivegreen4"), labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire"))+ ggtitle("b") +
    xlab("Generation") + ylab("Generation interval [years]") + xlim(c(40,60)) +
    facet_grid(Scenario ~ ., scales = "free_y") + theme(legend.position = "bottom", legend.text=element_text(size=12)) + 
    theme(axis.title.x=element_text(margin=margin(15,0,0,0),size=14), axis.title.y=element_text(size=14), axis.text=element_text(size=11), legend.title=element_text(size=12))
  number = number + 1
}
multiplot(plotList[[2]], plotList[[4]], plotList[[5]], plotList[[1]], plotList[[3]])
'''

names1 <- names
colnames(names1) <- c("Scenario",  "SCENARIO")
names1$Scenario <- gsub("OtherCowsGenGen", "OtherCowsGen", names1$Scenario)
GIa <- merge(GIa, names1, by="Scenario")

GIa$Scenario <- as.factor(GIa$Scenario)
GIa$SCENARIO <- as.factor(GIa$SCENARIO)
GIa$Gen <- as.numeric(as.character(GIa$Gen))
GIa$genInt <- as.numeric(GIa$genInt)
GIa$Scenario <- as.factor(GIa$Scenario)
GIPlot <- ggplot(GIa, aes(x=Gen, y=genInt, colour=Path, linetype =Path, group=Path)) + geom_line(aes(linetype=Path), size=1) + 
  scale_linetype_manual("", values=c("solid", "dashed", "dotted", "dotdash"), labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + 
  scale_colour_manual("", values=c("palevioletred2", "purple2", "royalblue3", "darkolivegreen4"), labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire"))+ ggtitle("b") +
  xlab("Generation") + ylab("Generation interval [years]") + xlim(c(40,60)) +
  facet_grid(SCENARIO ~ .) + theme(legend.position = "bottom", legend.text=element_text(size=12)) + 
  theme(axis.title.x=element_text(margin=margin(15,0,0,0),size=14), axis.title.y=element_text(size=14), 
        axis.text=element_text(size=11), legend.title=element_text(size=11))+
    guides(group=guide_legend(nrow=2), fill=guide_legend(nrow=2), colour=guide_legend(nrow=2), linetype=guide_legend(nrow=2)) +
  theme(legend.position="bottom")
#TO JE ZDJ COMBINED GAIN in GENINT
library(gridExtra)
library(grid)
gl <- lapply(1:1, function(ii) grobTree(rectGrob(), textGrob(ii)))
a <- ggplot() + theme(panel.background = element_blank())
grid.arrange(genGainPlot, GIPlot, layout_matrix = rbind(c(2,2,2,3,3)))

#tukaj pridobi številke za generation interval
GIAverages <- GI[GI$Gen %in% 40:60,]
GIAverages <- aggregate(GI$genInt ~ GI$path + GI$scenario, FUN=mean)
GIAverages <- GIAverages[order(GIAverages$`GI$path`),]
colnames(GIAverages) <- c("path", "scenario", "genInt")

GISum <- aggregate(GIAverages$genInt ~ GIAverages$scenario, FUN=sum)


  #############
  geom_segment(data=tmp5, mapping=aes(x=zSdGenicMax, xend=zSdGenicMin,                                      
                                      y=zMeanStartSdGenic, yend=zMeanEndSdGenic,                                      
                                      color=scenario3, linetype=scenario3)) +  
    geom_segment(data=tmp5, mapping=aes(x=zSdGenicMax, xend=zSdGenicMin,                                      
                                        y=zMeanStartSdGenic, yend=zMeanEndSdGenic,                                      
                                        color=scenario3, linetype=scenario3),               
                 arrow=arrow(), show.legend=FALSE)
  ###############

library(lme4)
lmList(zMean ~ SDGenicSt | scenario, data=TGVsAll)


#################
#plot genic and genetic variance through 60 generations
#################
TGVsAll <- TGVsAll[,c("Generation", "sd", "zSdGenic")]
a <- melt(TGVsAll, id.vars = 'Generation')
ggplot(data = a, aes(x=Generation, y=value, group = variable, colour=variable)) + geom_point() + geom_smooth() + 
  scale_colour_discrete("Variance", labels=c("Genetic", "Genic"))
Genetic <- ggplot(data = TGVsAll, aes(x=Generation, y=sd, group = scenario, colour=scenario)) + geom_point() + geom_smooth(se=FALSE) +
  scale_colour_discrete("Variance", labels=c("Class1", "Gen for all cows", "Gen for selection for progeny testing","Gen for bull dams",  "Gen for other cows")) + 
  ylab("Genetic variance")
Genic <- ggplot(data = TGVsAll, aes(x=Generation, y=zSdGenic, group = scenario, colour=scenario)) + geom_point() + geom_smooth(se=FALSE) +
  scale_colour_discrete("Variance", labels=c("Class1", "Gen for all cows", "Gen for selection for progeny testing","Gen for bull dams",  "Gen for other cows")) + 
  ylab("Genic variance")
