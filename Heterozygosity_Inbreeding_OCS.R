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


library(readr)
library(reshape2)


#TO JE SKRIPT ZA primerjavo OCS in osnovnih scenarijev (SU55PT, SU55GS, SU51GS)
#samo primerjava heterozigntonsto na QTNih


#povprečna hetero/homozigotnost --> dF / Ne
#generacije 1 - 20 so generacije 40 - 60
#Heterozygosity na QTNih
hetQTN <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_QTNsALL_14082018.csv")[]
hetQTN$Gen <- hetQTN$Gen + 40
hetQTN_OCS <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/MeanHet_ALL_11022019_OCS.csv")
hetQTN_OCS <- hetQTN_OCS[hetQTN_OCS$Marker == "QTN",]
hetQTN_OCS$Scenario <- hetQTN_OCS$degree
hetQTN_OCS <- hetQTN_OCS[,1:6]
#colnames(hetQTN) <-  c("Strategy",  "Scenario", "Rep", "Gen", "Het")
hetQTN_OCS$Strategy <- "OCS"
hetQTN_OCS$Marker <- "QTN"
hetQTN$Marker <- "QTN"
colnames(hetQTN_OCS)[5] <-  c("Het")
hetQTN <- hetQTN[,-5]
colnames(hetQTN)[5] <- "Het"
#hetQTN <- hetQTN[,c(colnames(hetQTN_OCS))]

#some generations are missing, some are duplicated
table(hetQTN_OCS$Rep, hetQTN_OCS$Scenario)

#remove duplicated, keep the second
hetQTN_OCS$id <- paste0(hetQTN_OCS$Scenario, hetQTN_OCS$Rep, hetQTN_OCS$Gen)
hetQTN_OCS[hetQTN_OCS$id %in% hetQTN_OCS$id[duplicated(hetQTN_OCS$id)],]
hetQTN_OCS[duplicated(hetQTN_OCS$id, fromLast=TRUE),] #odstrani prvo vrstico ponovljenih
hetQTN_OCS <- hetQTN_OCS[!(duplicated(hetQTN_OCS$id, fromLast=TRUE)),] #odstrani prvo vrstico ponovljenih
hetQTN_OCS <- hetQTN_OCS[,-7]

hetQTN_OCS$Scenario <- as.factor(hetQTN_OCS$Scenario)
hetQTN_OCS$Scenario <- as.factor(as.character(hetQTN_OCS$Scenario))
hetDF <- rbind(hetQTN, hetQTN_OCS)
table(hetDF$Scenario)


#ugotovili smo, da so vsi trije načini primerljivi, tako da obdržimo samo en rezultat

####OCS
table(hetDF$Scenario)
hetDF$Scenario <- as.character(hetDF$Scenario)
hetDF$Scenario <- revalue(hetDF$Scenario, c("Class" = "PT", "GenSLO" = "GS-PS", "OtherCowsGen" = "GS-C", "BmGen" = "GS-BD", "Gen" = "GS",
                                                "15"="15", "30"="30", "45"="45", "60"="60", "75"="75"))

hetDF$Group <- paste0(hetDF$Strategy, hetDF$Scenario)
hetDF <- hetDF[hetDF$Group %in% c("SU55PT","SU55GS", "SU51GS", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75"),]


#tukaj najprej izračunaj deltaF in NE po replikah
dF_DF <- data.frame(Strategy=NA, Scenario=NA, Rep=NA, Marker=NA, dF=NA, Ne=NA, Group=NA)
for (group in c("SU55PT","SU55GS", "SU51GS", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75")) {
  for (rep in 0:19) {
    for (marker in c( "QTN")) {
      #print(c(strategy, scenario, rep, marker))
      tmp <- hetDF [(hetDF$Group==group) & (hetDF$Rep==rep) & (hetDF$Marker==marker),]
      if (nrow(tmp) > 1) {
        
        tmp$y1 = log(tmp$Het)
        fit1 = MASS:::rlm(tmp$y1 ~ tmp$Gen, maxit=2000)
        dF <-  1 - exp(coef(fit1)[2])
        Ne <-  1 / (2 * dF)

        dF_DF <- rbind(dF_DF, c(as.character(unique(tmp$Strategy)), as.character(unique(tmp$Scenario)), rep, marker, dF, Ne, group))
      }
    }
  }
}

table(dF_DF$Scenario)
table(dF_DF$Strategy)
table(dF_DF$Group)

dF_DF$dF <- as.numeric(dF_DF$dF )
dF_DF$Ne <- as.numeric(dF_DF$Ne)
dF_DF <- dF_DF[dF_DF$Ne > 0,]
dF_DF <- dF_DF[-1,]


#tu pa sedaj izračunaj % razlike od base scenarija po replikah
HET <- data.frame()
for (group in c("SU55PT","SU55GS", "SU51GS", "OCS15", "OCS30", "OCS45", "OCS60", "OCS75")) {
  for (rep in 0:19) {
    #genetic gain
    base <- dF_DF$Ne[dF_DF$Scenario=="PT" & dF_DF$Strategy=="SU55" & dF_DF$Marker==marker & dF_DF$Rep==rep]
    basedF <- dF_DF$dF[dF_DF$Scenario=="PT" & dF_DF$Strategy=="SU55" & dF_DF$Marker==marker & dF_DF$Rep==rep]
    het <- dF_DF[dF_DF$Group==group & dF_DF$Marker==marker & dF_DF$Rep==rep,]
    het$per_Ne <- het$Ne / base
    het$per_dF <- het$dF / basedF
    
    HET <- rbind(HET, het)
  }
}



#sedaj pretvori v pravilne procente (razlika, ali je osnova 1 ali 0, ali hočeš povečanje ali zmanjšanje)
HET$per_Ne <- (HET$per_Ne)*100 - 100
HET$per_dF <- HET$per_dF*100 - 100

#HETavg <- aggregate(HET$dF ~ HET$Strategy + HET$Scenario + HET$Marker, FUN="mean")
#colnames(HETavg) <- c("Strategy", "Scenario", "Marker", "dF")

#NAJPREJ izračunaj povprečje in SD
##################################################################
#HET<- summarySE(HET, measurevar="Ne", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)] #to je za ABSOLUTNE vrednosti!
#HETabs <- summarySE(HET, measurevar="dF", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)] #to je za ABSOLUTNE vrednosti!
#to je za odstotek razlike v NE!
HETa <- summarySE(HET, measurevar="per_Ne", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)]
HETa_abs <- summarySE(HET, measurevar="Ne", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)]
#odstotek razlike zas delta F, ki je v bistvu delta coancestry (rate of coancestry) - ker je izračunano iz genomskega inbridinga
HETb <- summarySE(HET, measurevar="per_dF", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)]
HETb_abs <- summarySE(HET, measurevar="dF", groupvars=c("Strategy", "Scenario", "Marker"))[,c(1,2,3,5,6)]

colnames(HETa) <- c("Strategy", "Scenario", "Marker", "per_Ne", "per_NeSD")
colnames(HETa_abs) <- c("Strategy", "Scenario", "Marker", "Ne", "NeSD")
colnames(HETb) <- c("Strategy", "Scenario", "Marker", "per_dF", "per_dFSD")
colnames(HETb_abs) <- c("Strategy", "Scenario", "Marker", "dF", "dFSD")

#zaokroži
HETa$per_Ne <- round(HETa$per_Ne)
HETa$per_NeSD <- round(HETa$per_NeSD)
HETb$per_dF <- round(HETb$per_dF)
HETb$per_dFSD <- round(HETb$per_dFSD)

HETa <- merge(HETa, HETb, by=c("Strategy", "Scenario", "Marker"))
HETabs <- merge(HETa_abs, HETb_abs, by=c("Strategy", "Scenario", "Marker"))

HETa$Strategy <- factor(HETa$Strategy, levels =c("SU55", "SU51", "OCS"))
HETa$Scenario <- factor(HETa$Scenario, levels =c("PT", "GS", 15, 30, 45, 60, 75))
HETa <- HETa[order(HETa$Strategy, HETa$Scenario),]
HETa
HETabs$Strategy <- factor(HETabs$Strategy, levels =c("SU55", "SU51", "OCS"))
HETabs$Scenario <- factor(HETabs$Scenario, levels =c("PT", "GS", 15, 30, 45, 60, 75))
HETabs <- HETabs[order(HETabs$Strategy, HETabs$Scenario),]
HETabs

write.csv(HETa, "~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NEs_HET_relative_17012019.csv", row.names=FALSE, quote=FALSE)

#####################################################################3


#preveri značilnost med QTN, M, NTR
############################################################
HET$Group <- factor(HET$Group, level=c("SU55PT","SU55GS","SU51GS","OCS15","OCS30","OCS45","OCS60","OCS75"))

#Ne
model <- lm(per_Ne ~ Group, data=HET)
anova(model)
marginal = emmeans(model, ~ Group)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#deltaC
model <- lm(per_dF ~ Group, data=HET)
anova(model)
marginal = emmeans(model, ~ Group)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD


#abs
#Ne
model <- lm(Ne ~ Group, data=HET)
anova(model)
marginal = emmeans(model, ~ Group)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#deltaC
model <- lm(dF ~ Group, data=HET)
anova(model)
marginal = emmeans(model, ~ Group)
CLD = cld(marginal, 
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD






