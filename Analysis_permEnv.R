
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

library(ggplot2)
library(plyr)
pheno <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_pheno_29102019.csv", header=TRUE, sep=" ")
pheno$Ref <- "True"
phenoF <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_phenoFalse_29102019.csv", header=TRUE, sep=" ")
phenoF$Ref <- "False"
pheno <- rbind(pheno, phenoF)
class <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_permEnvClass_SU55_09102019.csv", header=TRUE, sep=" ")
class <- class[,-18]
class$Ref <- "Class"
pheno <- rbind(pheno, class)

pheno$scenario <- factor(pheno$scenario, levels =c("Class", "1_1", "1_2", "2_1"))
pheno$scenario <- revalue(pheno$scenario, c("Class" = "PT", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2","2_1" = "$G:$P = 2:1"))


class <- pheno[pheno$scenario == "PT",]
pheno <- pheno[pheno$scenario != "PT",]

class1 <- class
class1$scenario <- "$G:$P = 1:1"
class1$Ref <- "True"
class2 <- class
class2$scenario <- "$G:$P = 1:2"
class2$Ref <- "True"
class3 <- class
class3$scenario <- "$G:$P = 2:1"
class3$Ref <- "True"
class4 <- class
class4$scenario <- "$G:$P = 1:1"
class4$Ref <- "False"
class5 <- class
class5$scenario <- "$G:$P = 1:2"
class5$Ref <- "False"
class6 <- class
class6$scenario <- "$G:$P = 2:1"
class6$Ref <- "False"

CLASS <- rbind(class1, class2)
CLASS <- rbind(CLASS, class3)
CLASS <- rbind(CLASS, class4)
CLASS <- rbind(CLASS, class5)
CLASS <- rbind(CLASS, class6)
pheno <- rbind(pheno, CLASS)

pheno$Repeats <- as.factor(pheno$Repeats)
pheno$NoControls <- pheno$Repeats
pheno$BV <- ifelse(pheno$NoControls == 11, "Conventional", "Genomic")

pheno$Group <- paste(pheno$scenario, pheno$NoControls, sep="_")
finished <- unique(pheno$Group[pheno$Generation == 60])
pheno$scenario <- factor(pheno$scenario, levels = c("$G:$P = 2:1", "$G:$P = 1:1", "$G:$P = 1:2"))
#pheno <- pheno[pheno$Group %in% finished,]
phenoA <- summarySE(pheno, groupvars = c("Generation", "NoControls", "scenario", "Group", "BV", "Ref"), measurevar = "zMean")
phenoA$minzMean <- phenoA$zMean - phenoA$sd
phenoA$maxzMean <- phenoA$zMean + phenoA$sd

#png("/home/jana/Documents/Projects/inProgress/Phenotyping/Reallocation_plot.png", res=610, width=200, height=120, units="mm")

ggplot(data=phenoA, aes(x=Generation, y=zMean, group=NoControls, colour=NoControls, linetype = BV)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") +
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = NoControls),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario + Ref)

ggplot(data=phenoA[phenoA$Ref == "True",], aes(x=Generation, y=zMean, group=NoControls, colour=NoControls, linetype = BV)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") +
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = NoControls),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario + Ref)

#Plot genic variance
phenoAg <- summarySE(pheno, groupvars = c("Generation", "NoControls", "scenario", "Group", "BV", "Ref"), measurevar = "SDGenicSt")

ggplot(data=phenoAg, aes(x=Generation, y=SDGenicSt, group=NoControls, colour=NoControls, linetype = BV)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") +
  #geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = NoControls),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario + Ref)
#dev.off()


#naberi pheno v tabelo za excel
library(tidyr)
pheno60 <- phenoA[phenoA$Generation == 60,]
pheno60 <- pheno60[,c("NoControls", "scenario", "Ref", "zMean")] %>% spread(key = c(scenario), value = zMean)
write.table(pheno60, "/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/Results/GeneticGain.csv", quote=FALSE, row.names=FALSE, sep=";")

########################################################
########################################################
library(ggplot2)
accC <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_phenoClass_17102019.csv", header=TRUE)
accC1 <- accC
accC1$gp <- "$G:$P = 1:1"
accC1$Ref <- "True"
accC2 <- accC
accC2$gp <- "$G:$P = 1:2"
accC2$Ref <- "True"
accC3 <- accC
accC3$gp <- "$G:$P = 2:1"
accC3$Ref <- "True"
accC4 <- accC
accC4$gp <- "$G:$P = 1:1"
accC4$Ref <- "False"
accC5 <- accC
accC5$gp <- "$G:$P = 1:2"
accC5$Ref <- "False"
accC6 <- accC
accC6$gp <- "$G:$P = 2:1"
accC6$Ref <- "False"
acc <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_pheno_29102019.csv", header=TRUE)
acc$Ref <- "True"
accF <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_phenoFalse_29102019.csv", header=TRUE)
accF$Ref <- "False"
acc <- rbind(acc, accF)
acc <- rbind(acc, accC1)
acc <- rbind(acc, accC2)
acc <- rbind(acc, accC3)
acc <- rbind(acc, accC4)
acc <- rbind(acc, accC5)
acc <- rbind(acc, accC6)

acc$gp <- revalue(acc$gp, c("2_1" = "$G:$P = 2:1", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2"))
acc$Cor <- as.numeric(acc$Cor)
acc <- acc[acc$Gen %in% 40:60,]
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref"))
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref", "Gen"))
accA$NoControl <- as.factor(accA$NoControl)
accPlot <- accA[accA$AgeCat %in% c("genTest1", "telF1", "k3", "gpb4", "vhlevljeni1"),]
accPlot$AgeCat <- revalue(accPlot$AgeCat, c("genTest1" = "Male candidates", "telF1" = "Female candidates", "k3" = "Cows", "gpb4" = "gMales", "vhlevljeni1" = "Male candidates"))
accPlot$BV <- ifelse(accPlot$NoControl == 11, "Conventional", "Genomic")

write.table(accPlot[,c("AgeCat", "NoControl", "gp", "Cor", "Ref")] %>% spread(key = gp, value = Cor), 
            "/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/Results/Accuracy.csv", quote=FALSE, row.names=FALSE, sep=";")

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

ggplot(data=accPlot, aes(x=AgeCat, y=Cor, group=NoControl, colour=NoControl)) + 
  geom_point(size = 2)  + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "red")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  scale_x_discrete(breaks=unique(accPlot$AgeCat), 
                   labels=addline_format(c("Male candidates",  "Female candidates",  "Cows",  "gMales"))) + 
  facet_grid(. ~ gp + Ref)

accPlot$Gen <- as.factor(accPlot$Gen)
ggplot(data=accPlot[accPlot$Ref == "True",], aes(x=Gen, y=Cor, group=AgeCat, colour=AgeCat)) + 
  geom_line(size=1)  + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  ylab("Accuracy") +  xlab("Year") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  scale_x_discrete(labels=addline_format(c("Male candidates",  "Female candidates",  "Cows",  "gMales"))) + 
  facet_grid(rows = vars(gp), cols=vars(NoControl))



library(reshape)
refsize <- read.table("/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/ReferenceSize.txt", header=TRUE)
refsizeM <- melt(refsize, id.vars = "Generation")
refsizeM <- separate(refsizeM, col = variable, into = c("Rep", "Ref", "NoControls", "G", "P"), sep = "_")
refsizeM$gp <- paste(refsizeM$G, refsizeM$P, sep="_")
refsizeM$gp <- revalue(refsizeM$gp, c("2_1" = "$G:$P = 2:1", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2"))
refsizeM$NoControls <- as.factor(as.numeric(refsizeM$NoControls))

ggplot(data=refsizeM, aes(x=Generation, y=value, group=NoControls, colour=NoControls)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") + geom_hline(yintercept = 2000, colour="red") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ gp + Ref)

#put refsize and gain on the same plot
colnames(refsizeM)[8] <- "scenario"
colnames(refsizeM)[7] <- "TPSize"
gainref <- merge(phenoA[,c("Generation", "Ref", "NoControls", "scenario", "zMean")], 
                 refsizeM[,c("Generation", "Ref", "NoControls", "TPSize", "scenario")], by=c("Generation", "Ref", "NoControls", "scenario"))
gainrefM <- melt(gainref, id.vars = c("Generation", "Ref", "NoControls", "scenario"))
table(gainrefM$variable)
gainrefM$Group <- paste0(gainrefM$NoControls, "_", gainrefM$variable)
gainrefM$value[gainrefM$variable == "zMean"] <- gainrefM$value[gainrefM$variable == "zMean"] * 3000
gainref$NoControls <- as.factor(gainref$NoControls)
gainrefM$NoControls <- as.factor(gainrefM$NoControls)
gainrefM$variable <- revalue(gainrefM$variable, c("zMean" = "Gain", "TPSize" = "TP Size"))
gainrefM$RefValue <- paste0(gainrefM$Ref, "\n", gainrefM$variable)

#everything on one plot
ggplot(data=gainrefM, aes(x=Generation, y=value, group=NoControls, colour=NoControls)) +
    geom_line(size = 1) + 
    theme_bw(base_size=18, base_family="sans")  + 
    theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
          axis.text = element_text(size = 16)) + 
    scale_colour_manual("Number of recordings / lactation", 
                        values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", 
                                   "#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263")) + 
  ylab("Genetic mean") +
  scale_y_continuous(sec.axis = sec_axis(~ . / 3000)) + 
  #+ geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = NoControl),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(RefValue), cols=vars(scenario))


#put refsize and accuracy on the same plot
colnames(refsizeM)[8] <- "scenario"
colnames(refsizeM)[7] <- "TPSize"
colnames(accA)[5] <- "Generation"
colnames(accA)[2] <- "NoControls"
colnames(accA)[3] <- "scenario"
accref <- merge(accA[,c("Generation", "Ref", "NoControls", "scenario", "Cor", "AgeCat")], 
                refsizeM[,c("Generation", "Ref", "NoControls", "TPSize", "scenario")], by=c("Generation", "Ref", "NoControls", "scenario"))
accrefM <- melt(accref, id.vars = c("Generation", "Ref", "NoControls", "scenario", "AgeCat"))
table(accrefM$variable)
accrefM$variable <- revalue(accrefM$variable, c("Cor" = "Accuracy", "TPSize" = "TP Size"))
accrefM$Group <- paste0(accrefM$NoControls, "_", accrefM$variable)
accrefM$value[accrefM$variable == "Accuracy"] <- accrefM$value[accrefM$variable == "Accuracy"] * 25000
accref$NoControls <- as.factor(accref$NoControls)
accrefM$NoControls <- as.factor(accrefM$NoControls)

accrefM$RefValue <- paste0(accrefM$Ref, "\n", accrefM$variable)
accrefM$ScValue <- paste0(accrefM$scenario, "\n", accrefM$variable)
accrefM <- accrefM[accrefM$AgeCat %in% c("genTest1", "telF1", "k3", "gpb4"),]
accrefM$AgeCat <- revalue(accrefM$AgeCat, c("genTest1" = "Male candidates", "telF1" = "Female candidates", "k3" = "Cows", "gpb4" = "gMales"))

#everything on one plot
ggplot(data=accrefM[accrefM$Ref == "False",], aes(x=Generation, y=value, group=AgeCat, colour=AgeCat)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  ylab("Genetic mean") +
  scale_y_continuous(sec.axis = sec_axis(~ . / 25000)) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(ScValue), cols=vars(NoControls))


number = 1
plotList = list()
ggplot(data=gainref, aes(x=Generation, y=value, group=zMean, colour=Group)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", 
                                 "#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263")) + 
  ylab("Genetic mean") +
  scale_y_continuous(sec.axis = sec_axis(~ . / 3000)) + 
  #+ geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = NoControl),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(Ref), cols=vars(scenario))






mon <- read.csv("/home/jana/Documents/Projects/inProgress/Phenotyping/MoneySaved.csv")
mon$RecordingsRemoved <- as.numeric(mon$RecordingsRemoved)
ggplot(data = mon, aes(x=RecordingsRemoved, y=Genotypes, colour=P_G, group=P_G)) + geom_line(size = 1) +
  # xlab("Genic standard deviation") + ylab("Average True Genetic Value") + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16), legend.position = "top", 
        legend.text=element_text(size=16), legend.title=element_text(size=16),
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
        plot.margin = margin(t = 0, r = 10, b = 10, l = 10)) +
  scale_x_continuous(breaks=seq(1, 10, 1)) +
  guides(linetype=guide_legend(nrow=1, keyheight = unit(3, "cm"), keywidth = unit(3, "cm"), override.aes = list(alpha = 1, size=2))) 



############
#check significance
##1. GENETIC GAIN
library(emmeans)
TGV60$scenario <- as.factor(TGV60$scenario)

#ref = YES
TGV60_T <- pheno[pheno$Generation == 60 & pheno$Ref == "True",]
model_T <- lm(zMean ~ NoControls + scenario + NoControls:scenario, data=TGV60_T) #absolute
marginal = emmeans(model_T, ~ NoControls:scenario)
CLD = cld(marginal, by="NoControls",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#ref = NO
TGV60_F <- pheno[pheno$Generation == 60 & pheno$Ref == "False",]
model_F <- lm(zMean ~ NoControls + scenario + NoControls:scenario, data=TGV60_F) #absolute
marginal = emmeans(model_F, ~ NoControls:scenario)
CLD = cld(marginal, by="NoControls",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

##2. ACCURACY
library(emmeans)
acca <- summarySE()
TGV60$scenario <- as.factor(TGV60$scenario)

#ref = YES
TGV60_T <- pheno[pheno$Generation == 60 & pheno$Ref == "True",]
model_T <- lm(zMean ~ NoControls + scenario + NoControls:scenario, data=TGV60_T) #absolute
marginal = emmeans(model_T, ~ NoControls:scenario)
CLD = cld(marginal, by="NoControls",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#ref = NO
TGV60_F <- pheno[pheno$Generation == 60 & pheno$Ref == "False",]
model_F <- lm(zMean ~ NoControls + scenario + NoControls:scenario, data=TGV60_F) #absolute
marginal = emmeans(model_F, ~ NoControls:scenario)
CLD = cld(marginal, by="NoControls",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
CLD = cld(marginal, by="scenario",
          alpha   = 0.05,
          Letters = LETTERS,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
