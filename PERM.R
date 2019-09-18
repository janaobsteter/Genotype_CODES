males <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_MaleSU55_19082019.csv", header=TRUE)
var <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_varESU55_05092019.csv", header=TRUE)
var <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_varESU55_06092019.csv", header=TRUE)
pheno <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_pheno_12092019.csv", header=TRUE)

pheno$Repeats <- as.factor(pheno$Repeats)
pheno$NoControls <- pheno$Repeats
ggplot(data=pheno, aes(x=Generation, y=zMean, group=NoControls, colour=NoControls)) +
  geom_line() + 
  facet_grid(. ~ scenario)


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
acc$Cor <- as.numeric(acc$Cor)
acc <- acc[acc$Gen %in% 40:60,]

#pheno
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp"))
accA$NoControl <- as.factor(accA$NoControl)
ggplot(data=accA, aes(x=AgeCat, y=Cor, group=NoControl, colour=NoControl)) + geom_point()  + facet_grid(. ~ gp)

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


