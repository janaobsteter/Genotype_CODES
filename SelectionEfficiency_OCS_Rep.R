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
TGVsAll3 <- read.table("~/Genotipi/Genotipi_CODES//TGVsAll_10KRef_2501.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/TGVsAll_10KRef_2801_1Pb.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/OCS/TGVsAll_10KRef_OCS.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/TGVsAll_10KRef_OCS.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll1 <- read.table("~/TGVsAll_optiSel.csv", header=TRUE, sep=" ") #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
colnames(TGVsAll1)[11] <- "Degree"
TGVsAll2 <- read.table("~/TGVsAll_AlphaMate.csv", header=TRUE, sep=" ") #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
colnames(TGVsAll2)[11] <- "Degree"
TGVsAll1 <- read.table("~/TGVsAll_optiSel_Rep.csv", header=TRUE, sep=" ") #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll2 <- read.table("~/TGVsAll_AlphaMate_Rep.csv", header=TRUE, sep=" ") #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll1$Program <- "optiSel"
TGVsAll2$Program <- "AlphaMate"
TGVsAll <- read.table("~/Documents/PhD/Simulaton//TGVsAll_10KRef_2801_1Pb.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/Documents/PhD/Simulaton//TGVsAll_10KRef_2801_1Year.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
TGVsAll <- read.table("~/Genotipi/Genotipi_CODES//TGVsAll_10KRef_2801_1Year.csv", header=TRUE) #če imaš 60 generacij, je standardizirano na 1. generacijo!!!
#1 - 60 generacije
TGVsAll3 <- TGVsAll3[TGVsAll3$Rep %in% 0:2,]
TGVsAll3 <- TGVsAll3[TGVsAll3$scenario=="Gen",]
colnames(TGVsAll3)[10] <- "Degree"
TGVsAll3$Program <- "SU55"
TGVsAll <- read.table("~/Documents/PhD/Simulaton/RefPopSize/TGVsAll_10KRef_60Gen.csv", header=TRUE)
TGVsAll <- rbind(TGVsAll1, TGVsAll2)
TGVsAll <- rbind(TGVsAll, TGVsAll3)
TGVsAll$Group <- paste0(TGVsAll$Rep, TGVsAll$Program, TGVsAll$Degree)
TGVsAll$GroupP <- paste0(TGVsAll$Program, TGVsAll$Degree)
#TGVsAll <- TGVsAll[TGVsAll$Rep %in% c(0,1,2,3,5,6,7,8,9,10),]
TGVsAll <- TGVsAll[TGVsAll$Generation %in% 40:60,]




#Dodaj še nove spremenljivke


#standadise onto the 40 generation
st40 <- data.frame()
for (group in unique(TGVsAll$Group)) {
      repScenario <- TGVsAll[(TGVsAll$Group==group) & (TGVsAll$Generation %in% 40:60),]
      repScenario$zMean <- (repScenario$gvNormUnres1 - repScenario$gvNormUnres1[1]) / repScenario$sd[1]
      repScenario$SDGenicSt <- repScenario$AdditGenicVar1 / repScenario$AdditGenicVar1[1] 
      repScenario$SDSt <- repScenario$sd / repScenario$sd[1] 
      repScenario$zMeanGenic <- (repScenario$gvNormUnres1 - repScenario$gvNormUnres1[1]) / repScenario$AdditGenicVar1[1]
      repScenario$SDGenicStNeg <- 1- repScenario$AdditGenicVar1 / repScenario$AdditGenicVar1[1]
      #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
      #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
      st40 <- rbind(st40, repScenario)
}


#standadise onto the 40 generation
st40 <- data.frame()
for (rep in 0:2) {
  for (degree in c(15, 30, 45, 60, 75)) { #, "SU15", "SU51"
    for (program in c("optiSel", "AlphaMate", "SU55")) {
      repScenario <- TGVsAll[(TGVsAll$Rep==rep) & (TGVsAll$Program==program) & (TGVsAll$Degree==degree) & (TGVsAll$Generation %in% 40:60),]
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
Averages1 <- aggregate(TGVsAll$GENICSd_St ~TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages2 <- aggregate(TGVsAll$zMean ~TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages3 <- aggregate(TGVsAll$SDSt ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages4 <- aggregate(TGVsAll$zMean ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages5 <- aggregate(TGVsAll$zSdGenic ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages6 <- aggregate(TGVsAll$sd ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages7 <- aggregate(TGVsAll$MeanGENIC_genetic ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages8 <- aggregate(TGVsAll$GENICSd_St ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
Averages9 <- aggregate(TGVsAll$SDGenicStNeg ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
colnames(Averages1) <- c("Degree", "Generation", "SDGenic")
colnames(Averages2) <- c("Degree", "Generation", "MeanGenic")
colnames(Averages3) <- c("Degree", "Generation", "SDGenetic")
colnames(Averages4) <- c("Degree", "Generation", "MeanGenetic")
colnames(Averages5) <- c("Degree", "Generation", "SDGenic_noSt")
colnames(Averages6) <- c("Degree", "Generation", "SDGenetic_noSt")
colnames(Averages7) <- c("Degree", "Generation", "MeanGenicGenetic")
colnames(Averages8) <- c("Degree", "Generation", "SDGenic_Genetic")
colnames(Averages9) <- c("Degree", "Generation", "NegGenic")
Averages <- merge(Averages1, Averages2, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages3, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages4, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages5, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages6, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages7, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages8, by=c( "Degree", "Generation"))
Averages <- merge(Averages, Averages9, by=c( "Degree", "Generation"))

#AveragesA <- Averages
#to je max min za gensko varianco
maxmin <- data.frame(scenario=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (scenario in unique(Averages$scenario)) {
    maxmin[row,] <- c(scenario, min(Averages$SDGenic[Averages$scenario==scenario]), max(Averages$SDGenic[Averages$scenario==scenario]), min(Averages$MeanGenic[Averages$scenario==scenario]), max(Averages$MeanGenic[Averages$scenario==scenario]))
    row <- row +1
}
#to je max min za gensko varianco standardizirano na genetsko varianco
maxmin <- data.frame(degree=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (degree in unique(Averages$Degree)) {
    maxmin[row,] <- c(degree, (Averages$SDGenic_Genetic[(Averages$Degree==degree) & (Averages$Generation==60)]), (Averages$SDGenic_Genetic[(Averages$Degree==degree) & (Averages$Generation==40)]), 
                      (Averages$MeanGenicGenetic[(Averages$Degree==degree) & (Averages$Generation==40)]), (Averages$MeanGenicGenetic[(Averages$Degree==degree) & (Averages$Generation==60)]))
    row <- row +1
}
#to je max min za gensko varianco standardizirano na genetsko varianco
maxmin <- data.frame(degree=NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (degree in unique(Averages$Degree)) {
    maxmin[row,] <- c(degree, min(Averages$SDGenic_Genetic[Averages$Degree==degree]), max(Averages$SDGenic_Genetic[Averages$Degree==degree]), 
                      min(Averages$MeanGenicGenetic[Averages$Degree==degree]), max(Averages$MeanGenicGenetic[Averages$Degree==degree]))
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


#tukaj spravi gensko varianco in genetsko na isto točko

library(ggplot2)
#write.table(TGVsAll, "~/Simulation_Rep0/TGVsAll.csv", quote=FALSE, row.names=FALSE)
#TGVsAll <- read.table("~/Documents/WCGALP/TGVsAll.csv", header=TRUE)
#To je plot zMean (standardizirana na gensko variacno) na genetsko varianco



TGVsAll$Degree <- as.factor(TGVsAll$Degree)
maxmin$degree <- as.factor(maxmin$degree)

library(stringr)
maxmin$program <- str_sub(maxmin$degree, -1, -2)

TGVsAll$GroupP <- as.factor(TGVsAll$GroupP)
TGVsAll$Program <- as.factor(TGVsAll$Program)
maxmin$Degree <- c(rep(c(15, 30, 45, 60, 75), 2), "Gen")
maxmin$Program <- c(rep("AlphaMate", 5), rep("optiSel", 5), "SU55")
#o je z na genetsko standadrizirano gensko variacno
Year1Eff <- ggplot(data = TGVsAll, aes(x=GENICSd_St, y=MeanGENIC_genetic, group=GroupP, colour=Degree, linetype=Program)) + 
  geom_line(size=0.5, alpha=0.3) +
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,name="Pretvorjen/Izgubljen genski standardni odklon")) +
  ylim(-1,9) +coord_cartesian(xlim = c(1, 0.80)) +

  scale_colour_manual(breaks = c(15, 30, 45, 60, 75, "Gen"),  #optiSel15", "optiSel30", "optiSel45", "optiSel60", "optiSel75",  "AlphaMate15", "AlphaMate30", "AlphaMate45", "AlphaMate60", "AlphaMate75", "SU55Gen
                      "Tarčne stopinje", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1","black"), 
                      labels=c("15", "30", "45", "60", "75", "SLO")) + 
  xlab("Genski standardni odklon") + ylab("Povprečna prava genetska vrednost") + 
  scale_linetype_manual(breaks = c("optiSel", "AlphaMate", "SU55"), 
                        labels= c("optiSel", "AlphaMate", "SU55"), "Program", 
                        values=c("solid", "dashed","solid")) +
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=20), legend.position = "left",
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20)) +
   geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=Degree, group=degree, linetype=Program), arrow=arrow(), show.legend=TRUE, size=1.5)
#to je samo optiSel 
TGVsAll$degree <- substring(TGVsAll$Degree, 1,2)
TGVsAll$program <- substring(TGVsAll$Degree, 3)

ggplot(data = TGVsAll[TGVsAll$program=="optiSel",], aes(x=GENICSd_St, y=MeanGENIC_genetic, group=degree, colour=degree)) + 
  geom_line(size=1, alpha=0.5) +
  scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                                                                                     name="Converted/Lost genic standard deviation")) +
  #geom_smooth( se=FALSE, formula=y~x+1, method="lm") + 
  xlab("Generation") + ylab("True genetic value")  + ylim(-1,9) +coord_cartesian(xlim = c(1, 0.80)) +

  scale_colour_manual(breaks = c(15, 30, 45, 60, 75), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1", "black"), 
                      labels=c(15, 30, 45, 60, 75)) + 
  xlab("Genic standard deviation") + ylab("Average True Genetic Value") + 
  scale_linetype_manual(breaks = c(15, 30, 45, 60, 75), 
                        labels= c(15, 30, 45, 60, 75), "Scenario", 
                        values=c("solid",  "solid", "solid", "solid", "solid", "solid")) +
  theme(axis.text=element_text(size=18), legend.position = "left",
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=18)) +
  geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=degree, group=degree, linetype=degree), arrow=arrow(), show.legend=TRUE, size=1.5)


#To je plot zMean (standardizirana na GENETSKO variacno) 
ggplot(data = TGVsAll, aes(x=SDSt, y=zMean, group=Degree, colour=Degree, linetype=Degree)) + 
  geom_line(aes(linetype=TGVsAll$Degree), size=0.5, alpha=0.3) + scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                                                                                     name="Converted/Lost genic standard deviation")) +
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual(breaks = c(15, 30, 45, 60, 75), 
                        labels= c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"), "Scenario", 
                        values=c("solid", "dotted","dashed", "dotdash", "twodash")) + 
  scale_colour_manual(breaks = c(15, 30, 45, 60, 75), 
                      labels= c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"), "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1")) + 
  xlab("Genic standard deviation") + ylab("Average True Genetic Value") +  #+
  geom_segment(data=maxmin, mapping=aes(x=maxGeneticSD, xend=minGeneticSD,
                                        y=minTGV,  yend=maxTGV,                                    
                                        color=degree, linetype=degree, group=degree), arrow=arrow(), show.legend=FALSE, size=1.5)
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
genic <- ggplot(data = TGVsAll, aes(x=Generation, y=zSdGenic, group=GroupP, colour=GroupP, linetype=GroupP)) +  
  geom_line(aes(linetype=GroupP), size=1) + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual("Scenario", 
                        breaks = c("15optiSel", "30optiSel", "45optiSel", "60optiSel", "75optiSel", 
                                   "15AlphaMate", "30AlphaMate", "45AlphaMate", "60AlphaMate", "75AlphaMate", "SU55Gen"), 
                        values=c("dashed", "dashed", "dashed", "dashed", "dashed","solid", "solid","solid", "solid", "solid",  "solid"), 
                        labels=c("15optiSel", "30optiSel", "45optiSel", "60optiSel", "75optiSel", 
                                 "15AlphaMate", "30AlphaMate", "45AlphaMate", "60AlphaMate", "75AlphaMate", "SU55Gen")) + 
  scale_colour_manual("Scenario", 
                      breaks = c("15optiSel", "30optiSel", "45optiSel", "60optiSel", "75optiSel", 
                                 "15AlphaMate", "30AlphaMate", "45AlphaMate", "60AlphaMate", "75AlphaMate", "SU55Gen"), 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1","forestgreen", "dodgerblue2", "purple", "red3", "orange1", "black"), 
                      labels=c("15optiSel", "30optiSel", "45optiSel", "60optiSel", "75optiSel", 
                               "15AlphaMate", "30AlphaMate", "45AlphaMate", "60AlphaMate", "75AlphaMate", "SU55Gen")) + 
  xlab("Generation") + ylab("Additive Genic Variance")  +
  geom_segment(data=maxmin, mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                        y=minGenicSD,  yend=maxTGV,                                    
                                        color=degree, linetype=degree, group=degree),      arrow=arrow(), show.legend=FALSE)

MeanAverage <- aggregate(TGVsAll$zMean ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
MeanAverageSD <- aggregate(TGVsAll$zMean ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="sd")
MeanAverageGenic <- aggregate(TGVsAll$SDGenicSt ~ TGVsAll$GroupP + TGVsAll$Generation, FUN="mean")
colnames(MeanAverage) <- c("degree", "Generation", "MeanTGV")
colnames(MeanAverageSD) <- c("degree", "Generation", "SdTGV")
colnames(MeanAverageGenic) <- c("degree", "Generation", "Genic")
#genetic gain
library(nlme)
mean <- lmList(MeanTGV ~ Generation | degree, data=MeanAverage)
meanSD <- lmList(SdTGV ~ Generation | degree, data=MeanAverageSD)
meanGenic <- lmList(Genic ~ Generation | degree, data=MeanAverageGenic)
# meanDF <- datma.frame(degree=rownames(coef(mean)),coef(mean),check.names=FALSE)
colnames(meanDF) <- c("Degree", "Mean_Intercept", "Mean_Slope")
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
for (group in (unique(TGVsAll$GroupP))) {
  df <- TGVsAll[TGVsAll$GroupP==group,]
  fm1 <- lmList(zMeanGenic ~ SDGenicStNeg | Rep, data=df)
  tmp <- data.frame(Rep=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
  tmp$Scenario <- unique(df$GroupP)
  colnames(tmp) <- c("Rep", "Intercept", "Slope", "Scenario")
  regRep <- rbind(regRep, tmp)
}

#povprečje regresijskih koefificentov
avgSlope <- aggregate(regRep$Slope ~ regRep$Scenario, FUN="mean")
avgInt <- aggregate(regRep$Intercept ~ regRep$Scenario, FUN="mean")
avgReg <- cbind(avgInt, avgSlope[,2])
colnames(avgReg) <- c("Scenario", "Intercept", "Slope")
avgReg$method <- "Average regression"
##################################################3
#TO JE ŠE z 1 - sd --> da dobiš intercept 0
#to je povrepčje regresij
library(nlme)
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
regRep <- data.frame(Degree=NA, Intercept=NA, Slope=NA)


fm1 <- lmList(zMeanGenic ~ SDGenicStNeg | Degree, data=TGVsAll)
tmp <- data.frame(Degree=rownames(coef(fm1)),coef(fm1),check.names=FALSE)
colnames(tmp) <- c("Degree", "Intercept", "Slope")
regRep <- tmp




#zračunaj značilnost razlik naklonov med scenariji
library(lsmeans)
regRep$Degree <- as.factor(regRep$Degree)
regRep1 <- within(regRep, Degree <- relevel(Degree, ref = 15))
m1 <- lm(Slope~Degree,data=regRep)
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "Scenario")
contrast(m1.grid, method="pairwise")
contrast(m1S, method="eff")
summary(lm(Slope~Scenario,data=regRep1))

AllMethods <- rbind(Avg, Tot, avgReg)
AllMethodsO <- AllMethods[order(AllMethods$Scenario),]



#plot napredka po replikaj - GENETIC GAIN
regexp <- "[[:digit:]]+"
MeanAverage$Degree <- str_extract(MeanAverage$degree, regexp)
MeanAverage$Degree[MeanAverage$Degree==55] <- "Gen"
regexp <- "[[:alpha:]]+"
MeanAverage$Program <- str_extract(MeanAverage$degree, regexp)
MeanAverage$Program[MeanAverage$Program=="SU"] <- "SU55"


genGainPlot <- ggplot() + 
  scale_colour_manual(breaks = c(15, 30, 45, 60, 75, "Gen"),  #optiSel15", "optiSel30", "optiSel45", "optiSel60", "optiSel75",  "AlphaMate15", "AlphaMate30", "AlphaMate45", "AlphaMate60", "AlphaMate75", "SU55Gen
                      "Tarčne Stopinje", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1","black"), 
                      labels=c("15", "30", "45", "60", "75", "SLO_selekcija")) + 
  xlab("Genski standardni odklon") + ylab("Povprečna prava genetska vrednost") + 
  scale_linetype_manual(breaks = c("optiSel", "AlphaMate", "SU55"), 
                        labels= c("optiSel", "AlphaMate", "SLO_selekcija"), "Program", 
                        values=c("solid", "dashed","solid")) +
  xlab("Generacija") + ylab("Povprečna prava genetska vrednost") +
  geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=Degree, linetype=Program), size=1.2) + 
  ggtitle("a") + 
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=20), legend.position = "left",
                                         axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20)) #legend.position="bottom", 



regexp <- "[[:digit:]]+"
MeanAverageGenic$Degree <- str_extract(MeanAverageGenic$degree, regexp)
MeanAverageGenic$Degree[MeanAverageGenic$Degree==55] <- "Gen"
regexp <- "[[:alpha:]]+"
MeanAverageGenic$Program <- str_extract(MeanAverageGenic$degree, regexp)
MeanAverageGenic$Program[MeanAverageGenic$Program=="SU"] <- "SU55"
genicPlot <- ggplot() + 
  scale_linetype_manual(breaks = c("optiSel", "AlphaMate", "SU55"), 
                        labels= c("optiSel", "AlphaMate", "SLO_selekcija"), "Program", 
                        values=c("solid", "dashed","solid")) +
  scale_colour_manual(breaks = c(15, 30, 45, 60, 75, "Gen"),  #optiSel15", "optiSel30", "optiSel45", "optiSel60", "optiSel75",  "AlphaMate15", "AlphaMate30", "AlphaMate45", "AlphaMate60", "AlphaMate75", "SU55Gen
                      "Tarčne Stopinje", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1","black"), 
                      labels=c("15", "30", "45", "60", "75", "SLO_selekcija")) + 
  xlab("Genski standardni odklon") + ylab("Povprečna prava genetska vrednost") + 

  xlab("Generacija") + ylab("Genska varianca") +
  geom_line(data = MeanAverageGenic, aes(x=Generation, y=Genic, colour=Degree, linetype=Program), size=1.2) + 
  ggtitle("a") + 
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=20), legend.position = "left",
                                         axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20)) #legend.position="bottom", 


