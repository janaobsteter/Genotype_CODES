library(readr)
library(reshape2)



#povprečna hetero/homozigotnost --> dF / Ne
#generacije 1 - 20 so generacije 40 - 60
#Heterozygosity na QTNih
hetQTN <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_QTNsALL_14082018.csv")[, -1]
hetQTN <- read.csv("~/Documents/Projects/inProgress/GenomicStrategies_SireUSe/MeanHet_QTNsALL_OCS_11122018.csv")[, -1]
colnames(hetQTN) <-  c("Strategy",  "Scenario", "Rep", "Gen", "Het")
colnames(hetQTN) <-  c( "Scenario", "Rep", "Gen", "Het")
 hetQTN$Marker <- "QTN"
# nQTN <- 10000
# nAnim <- 8640
# hetQTN$HetCorr <- hetQTN$Het *  nQTN / nAnim
# write.csv(hetQTN, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_QTNsALL_14082018.csv", quote=FALSE, row.names=FALSE)

#popravi vrednosti!!!! NAREDILA SI VSOTO HETEROZIGOTNOST PO ZIVALI (apply - 2), NATO PA DELILA S ŠTEVILOM QTNom!!!


#Heterozygosity na nevtralnih lokusih
hetNTR <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_NeutralALL_14082018.csv")
# colnames(hetNTR) <-   c("Strategy",  "Scenario", "Rep", "Gen", "Het", "Marker", "HetCorr")
# hetNTR$Marker  <- "NTR"
# nNeutral <- 20000
# nAnim <- 8640
# hetNTR$HetCorr <- hetNTR$Het *  nNeutral / nAnim
# write.csv(hetNTR, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_NeutralALL_14082018.csv", quote=FALSE, row.names=FALSE)

#Heterozygosity na markerjih
hetM <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_MarkerALL_14082018.csv")
# colnames(hetM) <-  c("Strategy",  "Scenario", "Rep", "Gen", "Het", "Marker", "HetCorr")
# hetM$Marker  <- "M"
# nMarker <- 20000
# nAnim <- 8640
# hetM$HetCorr <- hetM$Het *  nMarker / nAnim
# write.csv(hetM, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_MarkerALL_14082018.csv", quote=FALSE, row.names=FALSE)

#združi podatke o heterozigotnost na QTNih in nevtralnih mestih
hetDF <- rbind(hetQTN, hetNTR)
hetDF <- rbind(hetDF, hetM)


'
#tukaj testiraj, ali pride do razlik z različnimi bazami (het0) ali pa različno regresijo
dF_DF <- data.frame(Strategy=NA, Scenario=NA, Rep=NA, Marker=NA, dF1=NA, Ne1=NA, dF2=NA, Ne2=NA, dF3=NA, Ne3=NA)
for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      for (marker in c("QTN", "NTR")) {
        tmp <- hetDF [(hetDF$Strategy==strategy) & (hetDF$Scenario==scenario) & (hetDF$Rep==rep) & (hetDF$marker==marker),]
        print(paste0(strategy, "_", scenario, "_", nrow(tmp)))
        het0 <- het_BurnIn$Het[het_BurnIn$marker==marker]
        het41 <- hetDF$Het[(hetDF$Strategy==strategy) & (hetDF$Scenario==scenario) & (hetDF$Rep==rep) & (hetDF$marker==marker) & (hetDF$Gen==1)]
        
        tmp$y1 = log(tmp$Het)
        fit1 = MASS:::rlm(tmp$y1 ~ tmp$Gen, maxit=2000)
        dF1 <-  1 - exp(coef(fit1)[2])
        Ne1 <-  1 / (2 * dF1)
  
        tmp$y2 = log(tmp$Het) - log(het0)
        fit2 = MASS:::rlm(tmp$y2 ~ tmp$Gen, maxit=2000)
        dF2 <-  1 - exp(coef(fit2)[2])
        Ne2 <-  1 / (2 * dF2)
        
        tmp$y3 = log(tmp$Het) - log(het41)
        fit3 = MASS:::rlm(tmp$y3 ~ tmp$Gen, maxit=2000)
        dF3 <-  1 - exp(coef(fit3)[2])
        Ne3 <-  1 / (2 * dF3)
        
        dF_DF <- rbind(dF_DF, c(strategy, scenario, rep, marker, dF1, Ne1, dF2, Ne2, dF3, Ne3))
      }
    }
  }
}
'
#ugotovili smo, da so vsi trije načini primerljivi, tako da obdržimo samo en rezultat

dF_DF <- data.frame(Strategy=NA, Scenario=NA, Rep=NA, Marker=NA, dF=NA, Ne=NA)
for (strategy in c("SU55", "SU15", "SU51")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      for (marker in c( "QTN", "NTR", "M")) {
        #print(c(strategy, scenario, rep, marker))
        tmp <- hetDF [(hetDF$Strategy==strategy) & (hetDF$Scenario==scenario) & (hetDF$Rep==rep) & (hetDF$Marker==marker),]
        if (nrow(tmp) > 1) {
  
          tmp$y1 = log(tmp$HetCorr)
          fit1 = MASS:::rlm(tmp$y1 ~ tmp$Gen, maxit=2000)
          dF <-  1 - exp(coef(fit1)[2])
          Ne <-  1 / (2 * dF)
          
          dF_DF <- rbind(dF_DF, c(strategy, scenario, rep, marker, dF, Ne))
        }
      }
    }
  }
}
dF_DF <- dF_DF[-1,]

####OCS
hetDF <- hetQTN
dF_DF <- data.frame(Strategy=NA, Scenario=NA, Rep=NA, Marker=NA, dF=NA, Ne=NA)
for (scenario in c(15, 30, 45, 60, 75)) {
  for (rep in 0:1) {
    for (marker in c( "QTN")) {
      #print(c(strategy, scenario, rep, marker))
      tmp <- hetDF [(hetDF$Scenario==scenario) & (hetDF$Rep==rep) & (hetDF$Marker==marker),]
      if (nrow(tmp) > 1) {
        
        tmp$y1 = log(tmp$HetCorr)
        fit1 = MASS:::rlm(tmp$y1 ~ tmp$Gen, maxit=2000)
        dF <-  1 - exp(coef(fit1)[2])t
        Ne <-  1 / (2 * dF)
        
        dF_DF <- rbind(dF_DF, c("OCS", scenario, rep, marker, dF, Ne))
      }
    }
  }
}
dF_DF <- dF_DF[-1,]


dF_DF$dF <- as.numeric(dF_DF$dF )
dF_DF$Ne <- as.numeric(dF_DF$Ne)
dF_DF <- dF_DF[dF_DF$Ne > 0,]



HET <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (marker in c("NTR", "M", "QTN")) {
      for (rep in 0:19) {
        #genetic gain
        base <- dF_DF$Ne[dF_DF$Scenario=="Class" & dF_DF$Strategy=="SU55" & dF_DF$Marker==marker & dF_DF$Rep==rep]
        het <- dF_DF[dF_DF$Scenario==scenario & dF_DF$Strategy==strategy & dF_DF$Marker==marker & dF_DF$Rep==rep,]
        het$per_het <- het$Ne / base
        
        HET <- rbind(HET, het)
      }
    }
  }
}

HET$per_het <- (1-HET$per_het)*100
HETavg <- aggregate(HET$dF ~ HET$Strategy + HET$Scenario + HET$Marker, FUN="mean")
colnames(HETavg) <- c("Strategy", "Scenario", "Marker", "dF")
#NAJPREJ izračunaj povprečje in SD
##################################################################
#HETa <- summarySE(HET, measurevar="Ne", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)] to je za ABSOLUTNE vrednosti
HETa <- summarySE(HET, measurevar="per_het", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)]
colnames(HETa) <- c("Strategy", "Scenario", "Marker", "per_het", "per_hetSD")
HETa$per_het <- round(HETa$per_het)
HETa$per_hetSD <- round(HETa$per_hetSD)

HETa$Strategy <- factor(HETa$Strategy, levels =c("SU55", "SU51", "SU15"))
HETa$Scenario <- factor(HETa$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
HETa$Marker <- factor(HETa$Marker, levels =c("NTR", "M", "QTN"))
HETa[order(HETa$Marker, HETa$Strategy, HETa$Scenario),]

write.csv(HETa, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NEs_HET_relative_21112018.csv", row.names=FALSE, quote=FALSE)

#####################################################################3


#preveri značilnost med QTN, M, NTR
############################################################
HET$Strategy <- factor(HET$Strategy, levels =c("SU55", "SU51", "SU15"))
HET$Scenario <- factor(HET$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))


model <- lm(per_het ~ Strategy + Scenario + Strategy + Marker + Marker: Strategy : Scenario, data=HET)
anova(model)
marginal = emmeans(model, ~ Strategy:Scenario:Marker)
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

#izračunaj značilnost
############################################################
HETQ <- HET[HET$Marker=="QTN",]
HETQ$Strategy <- factor(HETQ$Strategy, levels =c("SU55", "SU51", "SU15"))
HETQ$Scenario <- factor(HETQ$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))


model <- lm(per_het ~ Strategy + Scenario + Strategy : Scenario, data=HETQ)
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

'''
dF_DF <- dF_DF
#dF_DF <- dF_DF[-which((dF_DF$Marker=="NTR") & (dF_DF$Strategy=="10K_Ref_20Rep") & (dF_DF$Scenario=="Gen") & (dF_DF$Rep==6)),]
#dF_DF <- dF_DF[-which((dF_DF$Marker=="NTR") & (dF_DF$Strategy=="10K_Ref_20Rep") & (dF_DF$Scenario=="Gen") & (dF_DF$Rep==13)),]
dF_DFa <- aggregate(dF_DF$Ne ~ dF_DF$Strategy + dF_DF$Scenario + dF_DF$Marker, FUN="mean")
dF_DFa1 <- aggregate(dF_DF$dF ~ dF_DF$Strategy + dF_DF$Scenario + dF_DF$Marker, FUN="mean")
dF_DFa_SD <- aggregate(dF_DF$Ne ~ dF_DF$Strategy + dF_DF$Scenario + dF_DF$Marker, FUN="sd")
dF_DFa1_SD <- aggregate(dF_DF$dF ~ dF_DF$Strategy + dF_DF$Scenario + dF_DF$Marker, FUN="sd")
colnames(dF_DFa) <- c("Strategy", "Scenario","Marker",  "Ne")
colnames(dF_DFa1) <- c("Strategy", "Scenario","Marker",  "dF")
colnames(dF_DFa_SD) <- c("Strategy", "Scenario","Marker",  "Ne_sd")
colnames(dF_DFa1_SD) <- c("Strategy", "Scenario","Marker",  "dF_sd")
DFa <- merge(dF_DFa, dF_DFa1, by=c("Strategy", "Scenario","Marker"))
DFa <- merge(DFa, dF_DFa_SD, by=c("Strategy", "Scenario","Marker"))
DFa <- merge(DFa, dF_DFa1_SD, by=c("Strategy", "Scenario","Marker"))
scenarios <- c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")

DFa55 <- DFa[DFa$Strategy=="SU55",]
DFO <- DFa55[(order(match(DFa55$Scenario, scenarios))),c(1,2,3,4,6)]
DFO[,4:5] <- round(DFO[,4:5], 0)
DFO
DFa15 <- DFa[DFa$Strategy=="SU15",]
DFO <- DFa15[(order(match(DFa15$Scenario, scenarios))),c(1,2,3,4,6)]
DFO[,4:5] <- round(DFO[,4:5], 0)
DFO
DFa51 <- DFa[DFa$Strategy=="SU51",]
DFO <- DFa51[(order(match(DFa51$Scenario, scenarios))),c(1,2,3,4,6)]
DFO[,4:5] <- round(DFO[,4:5], 0)
DFO

#significance razlik v Ne
library(emmeans)
DFa$Scenario <- as.factor(DFa$Scenario)
DFa$Strategy <- as.factor(DFa$Strategy)
DFa$Marker <- as.factor(DFa$Marker)
DFa1 <- within(DFa, Scenario <- relevel(Scenario, ref = "Class"))
m1 <- lm(Ne~Marker ,data=DFa1[DFa1$Strategy=="SU55",])
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "Marker")
contrast(m1.grid, method="pairwise")
contrast(m1S, method="eff")
summary(lm(Ne~Scenario,data=DFa))

#significance razlik v Ne med strategijami
library(emmeans)
DFa$Scenario <- as.factor(DFa$Scenario)
DFa$Strategy <- as.factor(DFa$Strategy)
DFa$Marker <- as.factor(DFa$Marker)
DFa1 <- within(DFa, Strategy <- relevel(Strategy, ref = "SU55"))
m1 <- lm(Ne~ Strategy + Marker ,data=DFa1[DFa1$Strategy %in% c("SU55", "SU51") & DFa1$Scenario=="Class",])
m1.grid <- ref_grid(m1)
anova(m1)
m1S <- lsmeans(m1.grid, "Marker")
contrast(m1.grid, method="pairwise")
contrast(m1S, method="eff")
summary(lm(Ne~Scenario,data=DFa))
'''



##povprečje QTN, NTR and M - since there is no sigificant differences
DFmean <- aggregate(DFa$Ne ~DFa$Strategy + DFa$Scenario, FUN="mean")
DFmeanSD <- aggregate(DFa$Ne_sd ~DFa$Strategy + DFa$Scenario, FUN="mean")
colnames(DFmean) <- c("Strategy", "Scenario", "Ne")
colnames(DFmeanSD) <- c("Strategy", "Scenario", "Ne_sd")
DF_MEAN <- merge(DFmean, DFmeanSD, by=c("Strategy", "Scenario"))
DF_MEAN[,3:4] <- round(DF_MEAN[,3:4], 0)

write.table(DF_MEAN, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Ne_MEANMarkers.csv", quote=FALSE)


##združi z rodovniškim
INBa <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NEs_relative_21112018.csv")
INBa$Marker <- "Pedigree"
colnames(HETa)[colnames(HETa)=="per_het"] <- "per_inb"
colnames(HETa)[colnames(HETa)=="per_hetSD"] <- "per_inbSD"
INBALL <- rbind(INBa, HETa)


write.csv(INBALL, "/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PedHet_Inbreeding_Relative_21112018.csv", quote=FALSE, row.names = FALSE)
####################################################################


#izračunaj genomski F
#######################################################################

hetQTN$F <- NA
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      for (gen in 2:20) {
        #genetic gain
        hetQTN$F[hetQTN$Scenario==scenario & hetQTN$Strategy==strategy & hetQTN$Rep==rep & hetQTN$Gen==gen] <- 1 - 
          (hetQTN$HetCorr[hetQTN$Scenario==scenario & hetQTN$Strategy==strategy & hetQTN$Rep==rep & hetQTN$Gen==gen] / 
          hetQTN$HetCorr[hetQTN$Scenario==scenario & hetQTN$Strategy==strategy & hetQTN$Rep==rep & hetQTN$Gen==1])
      }
    }
  }
}

QTNa <- aggregate(hetQTN$F ~ hetQTN$Strategy + hetQTN$Scenario + hetQTN$Gen, FUN="mean")
colnames(QTNa) <- c("Strategy", "Scenario", "Generation", "F")
QTNa$Strategy <- factor(QTNa$Strategy, levels =c("SU55", "SU51", "SU15"))
QTNa$Scenario <- factor(QTNa$Scenario, levels =c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"))
QTNa$Scenario <- revalue(QTNa$Scenario, c("Class" = "PT", "GenSLO" = "GS-PS", "OtherCowsGen" = "GS-C", "BmGen" = "GS-BD", "Gen" = "GS"))

QTNa$Generation <- QTNa$Generation + 40
ggplot(data=QTNa, aes(x=Generation, y=F, group=Scenario, colour=Scenario)) + geom_path() + 
  facet_grid(Strategy ~ ., scales = "free_y")


'#tega zdj nimaš izračunanega
#agregiraj F po generacijah za plot F
df_QTN$F <- as.numeric(df_QTN$F)
hetA <- aggregate(df_QTN$F ~ df_QTN$Strategy + df_QTN$Scenario + df_QTN$Gen, FUN="mean")
hetASD <- aggregate(df_QTN$F ~ df_QTN$Strategy + df_QTN$Scenario + df_QTN$Gen, FUN="sd")

df_NTR$F <- as.numeric(df_NTR$F)
hetA1 <- aggregate(df_NTR$F ~ df_NTR$Strategy + df_NTR$Scenario + df_NTR$Gen, FUN="mean")
hetA1SD <- aggregate(df_NTR$F ~ df_NTR$Strategy + df_NTR$Scenario + df_NTR$Gen, FUN="sd")

colnames(hetA) <- c("Strategy", "Scenario", "Generation", "F")
colnames(hetASD) <- c("Strategy", "Scenario", "Generation", "sdF")
colnames(hetA1) <- c("Strategy", "Scenario", "Generation", "F")
colnames(hetA1SD) <- c("Strategy", "Scenario", "Generation", "sdF")

hetA <- merge(hetA, hetASD, by=c("Strategy", "Scenario", "Generation"))
hetA1 <- merge(hetA1, hetA1SD, by=c("Strategy", "Scenario", "Generation"))
hetA$method <- "Het_QTN"
hetA1$method <- "Het_NTR"
heta <- rbind(hetA, hetA1)

####################################################################
####################################################################
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


####################################################################
####################################################################
##primerjaj z rodovniškim - F po generacijah
inbRep$method <- "Ped"
inbRep <- inbRep[inbRep$Generation %in% 42:60,c(3,2,1,4,5,6)]


compareF <- rbind(inbRep, heta)
compareF$method <- as.factor(compareF$method)

#standardise - samo za Ped!!!
stCompare <- data.frame()
for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (method in c("Ped", "Het_QTN", "Het_NTR")) {
      tmp <- compareF[(compareF$Strategy==strategy) & (compareF$Scenario==scenario) & (compareF$method==method),]
      tmp$stF <- (tmp$F - tmp$F[1]) / tmp$sdF[1]
      stCompare <- rbind(stCompare, tmp) 
    }
  }
}


#standardiziraj - INBRIDINGA NE MOREŠ STANDARDIZIRATI, KER JE VERJETNOST!!!
#plot F po generacijah
compareF_str <- stCompare[stCompare$Strategy=="10K_Ref_20Rep",]
ggplot(data=compareF_str, aes(x=Generation, y=F, group=method, fill=method, colour=method)) + geom_path() +
  facet_grid(Scenario ~ ., scales = "free_y") + theme(legend.position = "right") +   theme(axis.text=element_text(size=16), legend.position = "left",
                                                                                           axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18), legend.title=element_text(size=20),
                                                                                           strip.text.y = element_text(size = 14)) + ggtitle("SU 5/5") 
'

##primerjaj z rodovniškim - Ne in deltaF
pedInb <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//NEs_deltaFs_14082018.csv")
HetInb$Scenario <- as.factor(HetInb$Scenario)
levels(HetInb$Scenario)
levels(HetInb$Scenario) <- c("GS-BD", "PT", "GS", "GS-PS", "GS-C")
HetInb$Strategy <- as.factor(HetInb$Strategy)
levels(HetInb$Strategy)
levels(HetInb$Strategy) <- c("SU 1/5", "SU 5/1", "SU 5/5")
pedInb$Marker <- "Pedigree"
pedInbA <- pedInb[pedInb$Interval.gen. =="41:60",c(3,2,8,4,6,5,7)]
colnames(pedInbA) <- colnames(DFa)


InbALL <- rbind(pedInbA, DFa)
InbALL$Ne <- round(InbALL$Ne, 0)
InbALL$Ne_sd <- round(InbALL$Ne_sd, 0)
InbALL$dF <- round(InbALL$dF, 5)
InbALL$dF_sd <- round(InbALL$dF_sd, 5)
write.csv(InbALL, "/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Compare_PedHet_Inbreeding.csv", quote=FALSE, row.names = FALSE)
####################################################################


