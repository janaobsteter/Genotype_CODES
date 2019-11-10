library(ggplot2)
library(plyr)
males <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_MaleSU55_19082019.csv", header=TRUE)
var <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_varESU55_05092019.csv", header=TRUE)
var <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_varESU55_06092019.csv", header=TRUE)
pheno <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_pheno_12092019.csv", header=TRUE)
pheno <- read.csv("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_comparisonPheno_rep0.csv", header=TRUE)

pheno <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_pheno_07112019.csv", header=TRUE, sep=" ")
class <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_permEnvClass_SU55_09102019.csv", header=TRUE, sep=" ")
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
#pheno <- pheno[pheno$Group %in% finished,]

png("/home/jana/Documents/Projects/inProgress/Phenotyping/Reallocation_plot.png", res=610, width=200, height=120, units="mm")

ggplot(data=pheno, aes(x=Generation, y=zMean, group=NoControls, colour=NoControls, linetype = BV)) +
  geom_line(size = 1) + theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
  values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", 
  values = c("dashed", "solid")) + 
  ylab("Genetic mean") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario)
dev.off()

levels(males$scenario)
levels(var$scenario)
males$scenario <- factor(males$scenario, levels=c("Class1_11", "Gen1_10MaleUpdate", "Gen1_11MaleNoUpdate", "Gen1_10Male", "Gen1_11"))
library(ggplot2)
var$scenario <- as.factor(var$scenario)
var$PlotGroup <- paste0(var$scenario, "_", var$Rep)
ggplot(var, aes(x=Generation, y=zMean, group=scenario, colour=scenario)) + geom_line() + 
  ylab("Genetic Mean") +    
  scale_colour_manual("Scenario", 
                      labels=c("varH=0.1", "varPE=0.2", "varH=0.5", "varPE=1","varPE=0.5",  "varPE=0.2"), 
                      values = c("red", "grey4", "red4", "grey60", "blue", "forestgreen")) +                            
 
  theme(axis.text=element_text(size=16),
  axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title=element_text(size=16)) + facet_grid(.~Rep)



acc <- read.table("~/Documents/Projects/inProgress/Phenotyping/Correlation_reference.csv", header=TRUE)[-1,]
acc$Model <- "PE"
accN <- read.table("~/Documents/Projects/inProgress/Phenotyping/Correlation_reference_noPE_permEnv.csv", header=TRUE)[-1,]
accN$Model <- "noPE"
accNew <- read.table("~/Documents/Projects/inProgress/Phenotyping/Correlation_reference_newModel.csv", header=TRUE)[-1,]
accNew$Model <- "CorrectedHY"

library(ggplot2)
acc <- read.table("~/Documents/Projects/inProgress/Phenotyping/CombinedAcc.csv", header=TRUE)
acc <- read.table("~/Documents/Projects/inProgress/Phenotyping/CombinedAcc_pheno.csv", header=TRUE)
acc <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_pheno_11102019.csv", header=TRUE)
acc$Cor <- as.numeric(acc$Cor)
acc <- acc[acc$Gen %in% 40:60,]

#pheno
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp"))
accA$NoControl <- as.factor(accA$NoControl)
accPlot <- accA[accA$AgeCat %in% c("genTest1", "telF1", "k3", "gpb4"),]
ggplot(data=accPlot, 
       aes(x=AgeCat, y=Cor, group=NoControl, colour=NoControl)) + geom_point(size = 2)  + 
  facet_grid(. ~ gp)

acc$varE <- as.factor(acc$varE)
levels(acc$varE)
ggplot(data=acc, aes(x=Gen, y=Cor, group=varE, colour=varE)) + geom_point() + 
   scale_colour_manual("Scenario", labels=c("varH=0.1", "varPE=0.2", "varH=0.5", "varPE=1","varPE=0.5",  "varPE=0.2"), values = c("red", "grey4", "red4", "grey60", "blue", "forestgreen")) +                            
   facet_grid(. ~ AgeCat + Ref)
accA <- summarySE(data=acc, measurevar = "Cor", groupvars = c("Ref", "AgeCat", "varE"))
ggplot(data=accA, aes(x=AgeCat, y=Cor, group=Ref, colour=Ref)) + geom_point()
ggplot(data=accA, aes(x=AgeCat, y=mean, group=varE, colour=varE)) + geom_point() + 
  scale_colour_manual("Scenario", labels=c("varH=0.1", "varPE=0.2", "varH=0.5", "varPE=1","varPE=0.5",  "varPE=0.2"), values = c("red", "grey4", "red4", "grey60", "blue", "forestgreen")) +                            
  facet_grid(. ~ Ref) 
levels(acc$Ref)

accA$Ref <- factor(accA$Ref, levels=c("1K", "10K"))
accA$varE <- as.factor(accA$varE)
ggplot(data=accA, aes(x=Ref, y=Cor, group=varE, colour=varE)) + geom_point() + 
  scale_colour_manual("Scenario", labels=c("c", "e", "d", "g", "f", "a","b"),  
                      values = c("black", "red", "blue", "pink", "orange", "purple", "forestgreen")) +
  facet_grid(. ~ AgeCat) +  theme(axis.text=element_text(size=16), legend.position = "top", 
                               axis.title=element_blank(), legend.text=element_text(size=16), legend.title=element_text(size=16),
                               plot.title = element_text(margin = margin(t = 0, r = 0, b = 40, l = 0), size=16, hjust=0.5),
                               plot.margin = margin(t = 0, r = 10, b = 10, l = 10))

ggplot(data=acc, aes(x=Ref, y=Cor, group=varE, colour=varE)) + geom_point() + 
  scale_colour_manual("Scenario", labels=c("c", "e", "d", "f", "a","b"), values = c("black", "red", "blue", "orange", "purple", "forestgreen")) +                            
  facet_grid(. ~ AgeCat) 

acc <- rbind(acc, accN)
acc <- rbind(acc, accNew)

acc$Rep <- as.factor(acc$Rep)
ggplot(data=acc, aes(x=nInd, y=Cor, group=Model, colour=Model)) + geom_point()
ggplot(data=acc, aes(x=nInd, y=Cor, group=Rep, colour=Rep)) + geom_point()

acc$size <- gsub("K", "", acc$Ref)
acc$size <- as.numeric(acc$size )* 1000
accA <- aggregate(acc$Cor ~ acc$size + acc$Model, FUN="mean")
colnames(accA) <- c("Ref", "Model", "Cor")

ggplot(data=accA, aes(x=Ref, y=Cor, group=Model, colour=Model)) + geom_line()
ggplot(data=accA, aes(x=Ref, y=Cor, group=Model, colour=Model)) + geom_point()
ggplot(data=acc, aes(x=size, y=Cor, group=Model, colour=Model)) + geom_point()





acc <- read.csv("~/Documents/Projects/inProgress/Phenotyping/AccuracyMale.csv", header=TRUE)
acc <- acc[acc$AgeCat %in% c("telF1", "genTest1", "potomciNP0"),]
ggplot(data=acc, aes(x=AgeCat, y=Cor, group=scenario, colour=scenario)) + geom_point()


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
