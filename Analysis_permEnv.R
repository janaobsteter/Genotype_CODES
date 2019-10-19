
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
males <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_MaleSU55_19082019.csv", header=TRUE)
var <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_varESU55_05092019.csv", header=TRUE)
var <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_varESU55_06092019.csv", header=TRUE)
pheno <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_pheno_12092019.csv", header=TRUE)
pheno <- read.csv("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_comparisonPheno_rep0.csv", header=TRUE)

pheno <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_pheno_17102019.csv", header=TRUE, sep=" ")
class <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_permEnvClass_SU55_17102019.csv", header=TRUE, sep=" ")
class <- class[,-18]
pheno <- rbind(pheno, class)

pheno$scenario <- factor(pheno$scenario, levels =c("Class", "1_1", "1_2", "2_1"))
pheno$scenario <- revalue(pheno$scenario, c("Class" = "PT", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2","2_1" = "$G:$P = 2:1"))


class <- pheno[pheno$scenario == "PT",]
pheno <- pheno[pheno$scenario != "PT",]

class1 <- class
class1$scenario <- "$G:$P = 1:1"
class2 <- class
class2$scenario <- "$G:$P = 1:2"
class3 <- class
class3$scenario <- "$G:$P = 2:1"

CLASS <- rbind(class1, class2)
CLASS <- rbind(CLASS, class3)
pheno <- rbind(pheno, CLASS)

pheno$Repeats <- as.factor(pheno$Repeats)
pheno$NoControls <- pheno$Repeats
pheno$BV <- ifelse(pheno$NoControls == 11, "Conventional", "Genomic")

pheno$Group <- paste(pheno$scenario, pheno$NoControls, sep="_")
finished <- unique(pheno$Group[pheno$Generation == 60])
pheno$scenario <- factor(pheno$scenario, levels = c("$G:$P = 2:1", "$G:$P = 1:1", "$G:$P = 1:2"))
#pheno <- pheno[pheno$Group %in% finished,]
phenoA <- summarySE(pheno, groupvars = c("Generation", "NoControls", "scenario", "Group", "BV"), measurevar = "zMean")
phenoA$minzMean <- phenoA$zMean - phenoA$sd
phenoA$maxzMean <- phenoA$zMean + phenoA$sd

png("/home/jana/Documents/Projects/inProgress/Phenotyping/Reallocation_plot.png", res=610, width=200, height=120, units="mm")

ggplot(data=phenoA, aes(x=Generation, y=zMean, group=NoControls, colour=NoControls, linetype = BV)) +
  geom_line(size = 1.5) + 
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
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, colour = NoControls, fill = NoControls),  linetype = 0, alpha = 0.2) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario)
dev.off()


#naberi pheno v tabelo za excel
library(tidyr)
pheno60 <- phenoA[phenoA$Generation == 60,]
pheno60 <- pheno60[,c("NoControls", "scenario", "zMean")] %>% spread(key = c(scenario), value = zMean)
write.table(pheno60, "/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/Results/GeneticGain_RefYes.csv", quote=FALSE, row.names=FALSE, sep=";")

library(ggplot2)
acc <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_pheno_17102019.csv", header=TRUE)
acc <- acc[acc$Gen %in% 40:60,]
acc$gp <- revalue(acc$gp, c("2_1" = "$G:$P = 2:1", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2"))
acc$Cor <- as.numeric(acc$Cor)
acc <- acc[acc$AgeCat %in% c("genTest1", "telF1", "k3", "gpb4"),]
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Gen"))
colnames(accA)[4] <- "Generation"
colnames(accA)[2] <- "NoControls"
accA <- merge(accA, refsizeM[,c("Generation", "NoControls", "value", "gp")], by=c("Generation", "gp", "NoControls"))
accA$NoControls <- as.factor(accA$NoControls)

#accuracy against Generation by NoControls
ggplot(data=accA, aes(x=Generation, y=Cor, group=NoControls, colour=NoControls)) + 
  geom_line()  + 
  #geom_errorbar(aes(ymin = Cor - sd, ymax = Cor + sd)) + 
  theme_bw(base_size=18, base_family="sans")  + 
  # scale_colour_manual("Number of recordings / lactation", 
  #                     values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) +
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  ylab("Genetic mean") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows = vars(AgeCat), cols = vars(gp))

#accuracy against  refSize
ggplot(data=accA, aes(x=value, y=Cor, group=AgeCat, colour=AgeCat)) + 
  geom_line()  + 
  #geom_errorbar(aes(ymin = Cor - sd, ymax = Cor + sd)) + 
  theme_bw(base_size=18, base_family="sans")  + 
  # scale_colour_manual("Number of recordings / lactation", 
  #                     values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) +
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  ylab("Genetic mean") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows = vars(NoControls), cols = vars(gp))


#pheno
acc <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_pheno_17102019.csv", header=TRUE)
accC <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_phenoClass_17102019.csv", header=TRUE)
accC1 <- accC
accC1$gp <- "$G:$P = 1:1"
accC2 <- accC
accC2$gp <- "$G:$P = 1:2"
accC3<- accC
accC3$gp <- "$G:$P = 2:1"

ACCC <- rbind(accC1, accC2)
ACCC <- rbind(ACCC, accC3)
acc <- rbind(acc, ACCC)

acc <- acc[acc$Gen %in% 40:60,]
acc$gp <- revalue(acc$gp, c("2_1" = "$G:$P = 2:1", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2"))
acc$Cor <- as.numeric(acc$Cor)
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp"))
accA$NoControl <- as.factor(accA$NoControl)
accA$BV <- ifelse(accA$NoControl == 11, "Conventional", "Genomic")
accPlot <- accA[accA$AgeCat %in% c("genTest1", "telF1", "k3", "gpb4", "cak5", "pb8"),]
accPlot$AgeCat <- revalue(accPlot$AgeCat, c("genTest1" = "Male candidates","cak5" = "Male candidates", "telF1" = "Female candidates", "k3" = "Cows", "gpb4" = "gMales","pb8" = "gMales"))
accPlot$AgeCat <- factor(accPlot$AgeCat, levels = c("Male candidates","Female candidates",  "Cows", "gMales"))

write.table(accPlot[,c("AgeCat", "NoControl", "gp", "Cor")] %>% spread(key = gp, value = Cor), 
            "/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/Results/Accuracy_RefYes.csv", quote=FALSE, row.names=FALSE, sep=";")

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


ggplot(data=accPlot, aes(x=AgeCat, y=Cor, group=NoControl, colour=NoControl)) + 
  geom_point(size = 2, aes(shape = BV)  )+ 
  #geom_errorbar(aes(ymin = Cor - sd, ymax = Cor + sd)) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "red")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  scale_x_discrete( 
                   labels=addline_format(c("Male cand.",  "Female cand.",  "Cows",  "gMales")), 
                   breaks = c("Male candidates","Female candidates",  "Cows", "gMales")) + 
  facet_grid(. ~ gp)

library(reshape)
refsize <- read.table("/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/ReferenceSize.txt", header=TRUE)
refsizeM <- melt(refsize, id.vars = "Generation")
refsizeM <- separate(refsizeM, col = variable, into = c("Ref", "NoControls", "G", "P"), sep = "_")
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
  ylab("Genetic mean") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ gp)



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
