library(readr)
library(reshape2)

MeanHet <- matrix(ncol=5)


for (rep in 0:1) {
  scenarioHet <- matrix(ncol=90)
  M <- readr::read_table("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen1/SimulatedData/AllIndividualsSnpChips/Last20Gen_small.txt",
                         col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
  for (gen in 1:20) {
    genM <- M[((gen-1)*8640 +1):(gen*8640),]
    
    MeanHet <- rbind(MeanHet, c("_strategy_", "_scenario_", rep, gen, mean(apply(X = genM, 2,  FUN = function(z) sum(z == 1)) / ncol(genM))))
    scenarioHet <- rbind(scenarioHet, c("_strategy_", "_scenario_", rep, gen, apply(X = genM, 2,  FUN = function(z) sum(z == 1)) / nrow(genM) ))
  }
  #write.csv(scenarioHet, paste0("/home/v1jobste/JanaO/Heterozygosity/ScenarioHET__strategy___scenario__Rep", rep, ".csv"), quote=FALSE)
}



write.csv(MeeanHet, paste0("/home/v1jobste/JanaO/Heterozygosity/ScenarioMeanHET__strategy___scenario_.csv"), quote=FALSE)


#povprečna hetero/homozigotnost --> dF / Ne
#generacije 1 - 20 so generacije 40 - 60
hetQTN <- read.csv("~/Heterozygosity_QTNMean.csv")
hetNTR <- read.csv("~/Heterozygosity_NTRMean.csv")

df_QTN <- data.frame(Strategy=NA, Scenario=NA, Rep=NA, Gen=NA, F=NA, dF=NA, Ne=NA)
for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      tmp1 <- hetQTN [(hetQTN$Strategy==strategy) & (hetQTN$Scenario==scenario) & (hetQTN$Rep==rep),]
      het0 <- hetQTN$Het[(hetQTN$Strategy==strategy) & (hetQTN$Scenario==scenario) & (hetQTN$Rep==rep) & (hetQTN$Gen==1)]
      for (gen in 2:20){
        dfGen1 <- (tmp1$Hom[tmp1$Gen==gen] - tmp1$Hom[tmp1$Gen==(gen-1)]) / (1 - tmp1$Hom[tmp1$Gen==(gen-1)])
        ne <- 1 / (2 * dfGen1)
        f <- 1 - (tmp1$Het[tmp1$Gen==gen] / het0)
        df_QTN <- rbind(df_QTN, c(strategy, scenario, rep, gen+40, f, dfGen1, ne))
      }
    }
  }
}

df_NTR <- data.frame(Strategy=NA, Scenario=NA, Rep=NA, Gen=NA, F=NA, dF=NA, Ne=NA)
for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      tmp <- hetNTR[(hetNTR$Strategy==strategy) & (hetNTR$Scenario==scenario) & (hetNTR$Rep==rep),]
      het0 <- hetNTR$Het[(hetNTR$Strategy==strategy) & (hetNTR$Scenario==scenario) & (hetNTR$Rep==rep) & (hetNTR$Gen==1)]
      for (gen in 2:20){
        dfGen <- (tmp$Hom[tmp$Gen==gen] - tmp$Hom[tmp$Gen==(gen-1)]) / (1 - tmp$Hom[tmp$Gen==(gen-1)])
        ne <- 1 / (2 * dfGen)
        f <- 1 - (tmp$Het[tmp$Gen==gen] / het0)
        df_NTR <- rbind(df_NTR, c(strategy, scenario, rep, gen+40, f, dfGen, ne))
      }
    }
  }
}

df_QTN <- df_QTN[-1,]
df_NTR <- df_NTR[-1,]



#regresija F na generacijo
df_NTR$F <- as.numeric(df_NTR$F)
df_NTR$Gen <- as.numeric(df_NTR$Gen)
df_NTR$Strategy <- as.factor(df_NTR$Strategy)

df_QTN$F <- as.numeric(df_QTN$F)
df_QTN$Gen <- as.numeric(df_QTN$Gen)
df_QTN$Strategy <- as.factor(df_QTN$Strategy)

df_QTN <- df_QTN[!(is.na(df_QTN$F)),]
df_NTR <- df_NTR[!(is.na(df_NTR$F)),]
ScenarioDF <- data.frame(Marker=NA, Strategy=NA, Scenario=NA, Rep=NA, deltaF=NA, Ne=NA)
for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      #na QTNih
      tmpQ <- df_QTN[(df_QTN==strategy) & (df_QTN$Scenario == scenario) & (df_QTN$Rep == rep),]
      if (nrow(tmpQ) == 19) {
        tmpQ$t = tmpQ$Gen - min(tmpQ$Gen) + 1
        tmpQ$y = log(1 - tmpQ$F)
        fitQ = MASS:::rlm(tmpQ$y ~ tmpQ$t, maxit=2000)
        # lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
        dFQ = 1 - exp(coef(fitQ)[2])
        NeQ = 1 / (2 * dFQ)
        
        #na nevtralnih lokusih
        tmpN <- df_NTR[(df_NTR==strategy) & (df_NTR$Scenario == scenario) & (df_NTR$Rep == rep),]
        tmpN$t = tmpN$Gen - min(tmpN$Gen) + 1
        tmpN$y = log(1 - tmpN$F)
        fitN = MASS:::rlm(tmpN$y ~ tmpN$t, maxit=2000)
        dFN <-  1 - exp(coef(fitN)[2])
        NeN <-  1 / (2 * dFN)
        ScenarioDF <- rbind(ScenarioDF, c("QTN", strategy, scenario, rep, dFQ, NeQ))
        ScenarioDF <- rbind(ScenarioDF, c("NTR", strategy, scenario, rep, dFN, NeN))
      }
    }
  }
}

ScenarioDF <- ScenarioDF[-1,]
ScenarioDF$deltaF <- as.numeric(ScenarioDF$deltaF )
ScenarioDF$Ne <- as.numeric(ScenarioDF$Ne)
ScenarioDF <- ScenarioDF[ScenarioDF$Ne > 0,]
#ScenarioDF <- ScenarioDF[-which((ScenarioDF$Marker=="NTR") & (ScenarioDF$Strategy=="10K_Ref_20Rep") & (ScenarioDF$Scenario=="Gen") & (ScenarioDF$Rep==6)),]
#ScenarioDF <- ScenarioDF[-which((ScenarioDF$Marker=="NTR") & (ScenarioDF$Strategy=="10K_Ref_20Rep") & (ScenarioDF$Scenario=="Gen") & (ScenarioDF$Rep==13)),]
ScenarioDFA <- aggregate(ScenarioDF$deltaF ~ScenarioDF$Strategy + ScenarioDF$Scenario + ScenarioDF$Marker, FUN="mean")
ScenarioDFA_Ne <- aggregate(ScenarioDF$Ne ~ScenarioDF$Strategy + ScenarioDF$Scenario + ScenarioDF$Marker, FUN="mean")
colnames(ScenarioDFA) <- c( "Strategy", "Scenario", "Marker","dF")
colnames(ScenarioDFA_Ne) <- c("Strategy", "Scenario","Marker",  "Ne")
HetInb <- merge(ScenarioDFA, ScenarioDFA_Ne, by=c("Strategy", "Scenario", "Marker"))


##primerjaj z rodovniškim
pedInb <- read.csv("~/Documents/PhD/Simulaton/NEs_deltaFs.csv")
HetInb$Scenario <- as.factor(HetInb$Scenario)
levels(HetInb$Scenario)
levels(HetInb$Scenario) <- c("GS-BD", "PT", "GS", "GS-PS", "GS-C")
HetInb$Strategy <- as.factor(HetInb$Strategy)
levels(HetInb$Strategy)
levels(HetInb$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5")

pedInb <- pedInb[pedInb$Interval.gen. =="41:60",c(2,3,5,4)]
HetInbQ <- HetInb[HetInb$Marker=="QTN",c(1,2,4,5)]
HetInbN <- HetInb[HetInb$Marker=="NTR",c(1,2,4,5)]
colnames(pedInb)[c(3,4)] <- c("dF_Ped", "Ne_Ped")
colnames(HetInbQ)[c(3,4)] <- c("dF_HetQ", "Ne_HetQ")
colnames(HetInbN)[c(3,4)] <- c("dF_HetN", "Ne_HetN")
dFcompare <- merge(pedInb, HetInbQ, by=c("Strategy", "Scenario"))
dFcompare <- merge(dFcompare, HetInbN, by=c("Strategy", "Scenario"))
dFcompare[,c(3, 5, 7)] <- round(dFcompare[,c(3,5,7)], 4)
dFcompare[,c(4, 6, 8)] <- round(dFcompare[,c(4,6,8)], 1)
write.csv(dFcompare, "~/Documents/PhD/Simulaton/Compare_PedHet_Inbreeding.csv", quote=FALSE, row.names = FALSE)

#Nariši plot
#QTN
df_QTN$dF <- as.numeric(df_QTN$dF)
df_QTN$Gen <- as.numeric(df_QTN$Gen)
df_QTN$Strategy <- as.factor(df_QTN$Strategy)
levels(df_QTN$Strategy)
levels(df_QTN$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5")

df_QTN$Group <- paste0(df_QTN$Scenario, df_QTN$Rep)
ggplot(data=df_QTN, aes(x=Gen, y=dF, group = Group, colour=Scenario, fill=Group)) + 
  geom_path() + 
  xlim (c(40, 60)) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) +  
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=16), legend.position = "left",
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20),
        strip.text.y = element_text(size = 14))   + 
  facet_grid(Strategy ~ ., scales = "free_y") + theme(legend.position = "right") 

#Nariši plot
#NTR
df_NTR$dF <- as.numeric(df_NTR$dF)
df_NTR$Gen <- as.numeric(df_NTR$Gen)
df_NTR$Strategy <- as.factor(df_NTR$Strategy)
levels(df_NTR$Strategy)
levels(df_NTR$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5")

df_NTR$Group <- paste0(df_NTR$Scenario, df_NTR$Rep)
ggplot(data=df_NTR, aes(x=Gen, y=dF, group = Group, colour=Scenario, fill=Group)) + 
  geom_path() + 
  xlim (c(40, 60)) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                      labels=c("PT", "GS-PS", "GS-C", "GS-BD", "GS")) +  
  guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) +
  theme(axis.text=element_text(size=16), legend.position = "left",
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20),
        strip.text.y = element_text(size = 14))   + 
  facet_grid(Strategy ~ ., scales = "free_y") + theme(legend.position = "right") 
