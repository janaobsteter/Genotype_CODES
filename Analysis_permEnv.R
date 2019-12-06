
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
pheno <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_pheno_15112019.csv", header=TRUE, sep=" ")
pheno$Ref <- "True"
phenoF <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_phenoFalse_15112019.csv", header=TRUE, sep=" ")
phenoF$Ref <- "False"
pheno <- rbind(pheno, phenoF)
table(pheno$Rep)
class <- read.csv("~/Documents/PhD/Projects/inProgress/Phenotyping/TGVsAll_permEnvClass_SU55_17102019.csv", header=TRUE, sep=" ")
class <- class[,-18]
class$Ref <- "Class"
pheno <- rbind(pheno, class)

pheno$scenario <- factor(pheno$scenario, levels =c("Class", "1_1", "1_2", "2_1"))
pheno$scenario <- revalue(pheno$scenario, c("Class" = "PT", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2","2_1" = "$G:$P = 2:1"))
table(pheno$Rep)

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
phenoA$Ref <- revalue(phenoA$Ref, c("True" = "With initial TP", "False" = "Without initial TP"))
phenoA$RealSc <- ifelse(phenoA$NoControl != 11, paste0("G", phenoA$NoControl), paste0("C", phenoA$NoControl))
phenoA$RealSc <- factor(phenoA$RealSc, levels = c("G1", "G2", "G5", "G8", "G9", "G10", "C11"))


ggplot(data=phenoA, aes(x=Generation, y=zMean, group=RealSc, colour=RealSc, linetype = BV)) +
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
  facet_grid(rows = vars(scenario), cols = vars(Ref))

#z začetno referenco
ggplot(data=phenoA[phenoA$Ref == "With initial TP",], aes(x=Generation, y=zMean, group=RealSc, colour=RealSc, linetype = BV)) +
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
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=2, keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  guides(colour=guide_legend(nrow=2, override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario)

#z začetno referenco _ SLO
ggplot(data=phenoA[phenoA$Ref == "With initial TP",], aes(x=Generation, y=zMean, group=RealSc, colour=RealSc, linetype = BV)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Scenarij", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Scenarij", 
                    values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", labels = c("Klasična", "Genomska"),
                        values = c("dashed", "solid")) + 
  ylab("Genetski napredek") + xlab("Generacija") + 
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=2, keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  guides(colour=guide_legend(nrow=2, override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario)


#brez začetne reference
ggplot(data=phenoA[phenoA$Ref == "Without initial TP",], aes(x=Generation, y=zMean, group=RealSc, colour=RealSc, linetype = BV)) +
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
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario + Ref)

#brez začetne reference - SLO
ggplot(data=phenoA[phenoA$Ref == "Without initial TP",], aes(x=Generation, y=zMean, group=RealSc, colour=RealSc, linetype = BV)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Scenarij", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Scenarij", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", labels = c("Klasična", "Genomska"),
                        values = c("dashed", "solid")) + 
  ylab("Genetski napredek") + xlab("Generacija") +
  geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ scenario)

#plot just gain in generation 60
ggplot(data=phenoA[phenoA$Generation == 60,], aes(x=NoControls, y=zMean, group=scenario, colour=scenario)) +
  geom_point(size = 2) + geom_line() + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("red", "blue", "green")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") +
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ Ref)


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

pheno60$zMean <- round(pheno60$zMean, 2)
pheno60$sd <- round(pheno60$sd, 2)
pheno60[pheno60$Ref =="With initial TP" & pheno60$scenario == "$G:$P = 1:1",]
pheno60[pheno60$Ref =="With initial TP" & pheno60$scenario == "$G:$P = 1:2",]
pheno60[pheno60$Ref =="With initial TP" & pheno60$scenario == "$G:$P = 2:1",]
pheno60[pheno60$Ref =="Without initial TP" & pheno60$scenario == "$G:$P = 1:1",]
pheno60[pheno60$Ref =="Without initial TP" & pheno60$scenario == "$G:$P = 1:2",]
pheno60[pheno60$Ref =="Without initial TP" & pheno60$scenario == "$G:$P = 2:1",]

ggplot(data=phPar, aes(x=Total, y=zMean, group=scenario, colour=scenario)) +
  geom_point(size = 2) + geom_line() + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("red", "blue", "green")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Genetic mean") +
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ Ref)

#################################################
#################################################
pheno60 <- phenoA[phenoA$Generation == 60,]
par <- read.csv("Genotipi/Genotipi_CODES/ParameterFile_Simulation.csv")
par <- par[,c(2, 3, 4, 5, 6, 7, 8)]
par$Males <- par$PotomciNPPerYear + par$telMPerYear
library(reshape)
library(ggplot2)
PAR <- par[,c(2,3,6,7,8)]
PARm <- melt(PAR, id.vars = c("G_P", "NoControl"))
PARm$NoControl <- as.numeric(PARm$NoControl)
PARm$NoControl <- as.numeric(PARm$NoControl)
PARm$value <- as.numeric(PARm$value)
PARm$variable <- factor(PARm$variable, c("Total", "FemalePerYear", "Males"))
library(plyr)
PARm$G_P <- revalue(PARm$G_P, c("1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2", "2_1" = "$G:$P = 2:1"))
PARm$G_P <- factor(PARm$G_P, c("$G:$P = 2:1", "$G:$P = 1:1", "$G:$P = 1:2"))

ggplot(PARm, aes(x = NoControl, y = value, colour=variable)) + geom_point(size = 2) + geom_line(size = 1) + 
  ylab("Število genotipiziranih živali") + theme_bw(base_size = 20) + 
  theme(axis.text = element_text(size = 18), legend.text = element_text(size = 18)) + 
  scale_colour_manual("", labels = c("Skupno", "Ženske", "Moški"), values = c("black", "red", "skyblue3")) + 
  scale_x_reverse("Število meritev / laktacijo", breaks = c(10, 9, 8, 5, 2, 1)) + 
  facet_grid(.~G_P)


colnames(par)[1:3] <- c("Ref", "scenario", "NoControls")
par$scenario <- revalue(par$scenario, c("2_1" = "$G:$P = 2:1", "1_1" = "$G:$P = 1:1", "1_2" = "$G:$P = 1:2"))
phPar <- merge(pheno60, par, by=c("Ref", "scenario", "NoControls"))
phPar$MalePerYear <- phPar$PotomciNPPerYear + phPar$telMPerYear


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

acc <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_pheno_15112019.csv", header=TRUE)
acc$Ref <- "True"
accF <- read.table("~/Documents/PhD/Projects/inProgress/Phenotyping/CombinedAcc_phenoFalse_15112019.csv", header=TRUE)
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
#accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref"))
#accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref", "Gen"))
acc$NoControl <- as.factor(acc$NoControl)
accPlot <- acc[acc$AgeCat %in% c("genTest1", "telF1",  "k3", "k4", "k5", "k6", "pBM4", "pBM5", "pBM6", "bm7", "gpb2", "gpb3", "gpb4", "gpb5", "gpb6", "pb6", "pb7", "pb8", "pbv9", "vhlevljeni1"),]
accPlot$AgeCat <- revalue(accPlot$AgeCat, c("genTest1" = "Male candidates", "telF1" = "Female candidates", 
                                            "k3" = "Mothers", 
                                            "k4" = "Mothers", 
                                            "k5" = "Mothers",
                                            "k6" = "Mothers", 
                                            "pBM4" = "Mothers", 
                                            "pBM5" = "Mothers", 
                                            "pBM6" = "Mothers", 
                                            "bm7" = "Mothers", 
                                            "gpb2" = "Fathers", 
                                            "gpb3" = "Fathers", 
                                            "gpb4" = "Fathers", 
                                            "gpb5" = "Fathers", 
                                            "gpb6" = "Fathers", 
                                            "pb6" = "Fathers", 
                                            "pb7" = "Fathers", 
                                            "pb8" = "Fathers", 
                                            "pb9" = "Fathers", 
                                            "vhlevljeni1" = "Male candidates"))
accA <- summarySE(data = accPlot, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref", "Rep"))
accA <- summarySE(data = accA, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref"))
accPlotA <- accA
accPlotA$BV <- ifelse(accPlotA$NoControl == 11, "Conventional", "Genomic")

write.table(accPlotA[,c("AgeCat", "NoControl", "gp", "Cor", "Ref")] %>% spread(key = gp, value = Cor), 
            "/home/jana/Documents/PhD/Projects/inProgress/Phenotyping/Results/Accuracy.csv", quote=FALSE, row.names=FALSE, sep=";")

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

ggplot(data=accPlotA, aes(x=AgeCat, y=Cor, group=NoControl, colour=NoControl)) + 
  geom_point(size = 2)  + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Number of recordings / lactation", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "red")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Accuracy") + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  scale_x_discrete(breaks=unique(accPlot$AgeCat), 
                   labels=addline_format(c("Male candidates",  "Female candidates",  "Cows",  "gMales"))) + 
  facet_grid(. ~ gp + Ref)

#the plot with lines
accPlotA$AgeCat <- factor(accPlotA$AgeCat, levels =  c("Male candidates",  "Female candidates",   "Fathers", "Mothers" ))
accPlotA$Ref <- revalue(accPlotA$Ref, c("True" = "With initial TP", "False" = "Without initial TP"))

accPlotA$Cor <- round(accPlotA$Cor, 2)
accPlotA$sd <- round(accPlotA$sd, 2)
accPlotA[accPlotA$Ref  == "With initial TP" & accPlotA$gp == "$G:$P = 1:1",]
accPlotA[accPlotA$Ref  == "With initial TP" & accPlotA$gp == "$G:$P = 2:1",]
accPlotA[accPlotA$Ref  == "With initial TP" & accPlotA$gp == "$G:$P = 1:2",]
accPlotA[accPlotA$Ref  == "Without initial TP" & accPlotA$gp == "$G:$P = 2:1",]
accPlotA[accPlotA$Ref  == "Without initial TP" & accPlotA$gp == "$G:$P = 1:1",]
accPlotA[accPlotA$Ref  == "Without initial TP" & accPlotA$gp == "$G:$P = 1:2",]

accPlotA$RealSc <- ifelse(accPlotA$NoControl != 11, paste0("G", accPlotA$NoControl), paste0("C", accPlotA$NoControl))
accPlotA$RealSc <- factor(accPlotA$RealSc, levels = c("G1", "G2", "G5", "G8", "G9", "G10", "C11"))

#just for GP 1:1
# & accPlot$NoControl != 11,
ggplot() + #
  geom_point(data=accPlotA[accPlotA$gp == "$G:$P = 1:1",], aes(x=RealSc, y=Cor, group=AgeCat, colour=AgeCat, linetype = AgeCat), size = 3)  +  
  geom_line(data=accPlotA[accPlotA$gp == "$G:$P = 1:1" & accPlotA$NoControl !=11,], aes(x=RealSc, y=Cor, group=AgeCat, colour=AgeCat, linetype = AgeCat), size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
  #geom_errorbar(aes(ymin = Cor - sd, ymax = Cor + sd), alpha = 0.5) + 
  scale_colour_manual("Animal group", 
                      values = c("skyblue2", "#e57691", "#0a488e", "#8e0a2a"),
                      breaks =  c("Male candidates",  "Female candidates",   "Fathers", "Mothers" ),
                      labels = c("Male candidates", "Female candidates", "Sires", "Dams")) + 
  scale_linetype_manual("Animal group", values=c("solid", "dashed", "solid", "dashed"),
                        breaks =  c("Male candidates",  "Female candidates",   "Fathers", "Mothers" ),
                        labels = c("Male candidates", "Female candidates", "Sires", "Dams")) +
  ylab("Accuracy") + xlab("Scenario") + 
  guides(colour=guide_legend(nrow=2, keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~  Ref + gp) +   theme(strip.text = element_text(size = 16))

#just proven males with errorbars
ggplot(data=accPlotA[accPlotA$AgeCat == "gMales",], aes(x=NoControl, y=Cor, group=AgeCat, colour=AgeCat, linetype = AgeCat)) + 
  geom_point(size = 3)  +  geom_line() + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
  geom_errorbar(aes(ymin = Cor - sd, ymax = Cor + sd)) + 
  scale_colour_manual("Animal group", 
                      values = c("skyblue2", "#e57691", "#0a488e", "#8e0a2a"),
                      breaks =  c("Male candidates",  "Female candidates",   "gMales", "Cows" ),
                      labels = c("Male candidates", "Female candidates", "Sires", "Dams")) + 
  scale_linetype_manual("Animal group", values=c("solid", "dashed", "solid", "dashed"),
                        breaks =  c("Male candidates",  "Female candidates",   "gMales", "Cows" ),
                        labels = c("Male candidates", "Female candidates", "Sires", "Dams")) +
  ylab("Genetic mean") + xlab("Number of phenotype recordings") + 
  guides(colour=guide_legend(nrow=2, keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~  Ref + gp)

accGen <- summarySE(data = accPlot, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref", "Gen"))
accGen$Gen <- as.factor(accGen$Gen)
ggplot(data=accGen[accGen$Ref == "True",], aes(x=Gen, y=Cor, group=AgeCat, colour=AgeCat)) + 
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
refsizeM$RealSc <- ifelse(refsizeM$NoControl != 11, paste0("G", refsizeM$NoControl), paste0("C", refsizeM$NoControl))
refsizeM$RealSc <- factor(refsizeM$RealSc, levels = c("G1", "G2", "G5", "G8", "G9", "G10", "C11"))

ggplot(data=refsizeM[refsizeM$Ref  == "True",], aes(x=Generation, y=value, group=RealSc, colour=RealSc)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16)) + 
  scale_colour_manual("Scenarij", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", 
                        values = c("dashed", "solid")) + 
  ylab("Število živali v referenci") + xlab("Generacija") +
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(. ~ gp)

ggplot(data=refsizeM[refsizeM$Ref  == "True",], aes(x=Generation, y=value, group=RealSc, colour=RealSc)) +
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
refsizeM$Ref <- revalue(refsizeM$Ref, c("True" = "With initial TP", "False" = "Without initial TP"))

colnames(refsizeM)[8] <- "scenario"
colnames(refsizeM)[7] <- "TPSize"
gainref <- merge(phenoA[,c("Generation", "Ref", "NoControls", "scenario", "zMean")], 
                 refsizeM[,c("Generation", "Ref", "NoControls", "TPSize", "scenario")], by=c("Generation", "Ref", "NoControls", "scenario"), all.x = TRUE)
gainrefM <- melt(gainref, id.vars = c("Generation", "Ref", "NoControls", "scenario"))
gainrefSD <- phenoA[,c("Generation", "Ref", "NoControls", "scenario", "sd")] 
gainrefSD$variable ="zMean"
gainrefSDM <- melt(gainrefSD, id.vars = c("Generation", "Ref", "NoControls", "scenario", "variable"))
gainrefSDM <- gainrefSDM[,-6]
colnames(gainrefSDM)[6] <- "SD"

table(gainrefM$variable)
table(gainrefSDM$variable)
gainrefM <- merge(gainrefM, gainrefSDM, by=c("Generation", "Ref", "NoControls", "scenario", "variable"), all.x = TRUE)
table(gainrefM$variable)
gainrefM$Group <- paste0(gainrefM$NoControls, "_", gainrefM$variable)
gainrefM$value[gainrefM$variable == "zMean"] <- gainrefM$value[gainrefM$variable == "zMean"] 
gainref$NoControls <- as.factor(gainref$NoControls)
gainrefM$NoControls <- as.factor(gainrefM$NoControls)
gainrefM$variable <- revalue(gainrefM$variable, c("zMean" = "Genetic mean", "TPSize" = "TP Size"))
gainrefM$RefValue <- paste0(gainrefM$Ref, "\n", gainrefM$variable)
gainrefM$SD[gainrefM$variable == "TPSize"] <- 0

#everything on one plot
#ref == yes
gainrefM$BV <- ifelse(gainrefM$NoControls == 11,"Conventional",  "Genomic")
gainrefM$RealSc <- ifelse(gainrefM$NoControl != 11, paste0("G", gainrefM$NoControl), paste0("C", gainrefM$NoControl))
gainrefM$RealSc <- factor(gainrefM$RealSc, levels = c("G1", "G2", "G5", "G8", "G9", "G10", "C11"))

g1 <- ggplot(data=gainrefM[gainrefM$Ref == "With initial TP",], aes(x=Generation, y=value, group=RealSc, colour=RealSc, linetype = BV)) +
    geom_line(size = 1) + 
    theme_bw(base_size=18, base_family="sans")  + 
    theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
          axis.text = element_text(size = 16), axis.title.y = element_blank()) + 
    scale_colour_manual("Scenario", 
                        values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Scenario", 
                        values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", values = c("dashed", "solid")) + xlab("Year") + 
  geom_ribbon(aes(ymin = value - SD, ymax = value + SD, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=2, keyheight = unit(.8, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  guides(colour=guide_legend(byrow = TRUE, nrow=2, keywidth = unit(1.2, "cm"),override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(variable), cols=vars(scenario), scales = "free_y")+
  theme(strip.text = element_text(size = 16))
#SLO
gainrefM$variable <- revalue(gainrefM$variable, c("TP Size" = "Velikost RP", "Genetic mean" = "Genetski napredek"))
g1 <- ggplot(data=gainrefM[gainrefM$Ref == "With initial TP",], aes(x=Generation, y=value, group=RealSc, colour=RealSc, linetype = BV)) +
    geom_line(size = 1) + 
    theme_bw(base_size=18, base_family="sans")  + 
    theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
          axis.text = element_text(size = 16), axis.title.y = element_blank()) + 
    scale_colour_manual("Scenarij", 
                        values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Scenarij", 
                        values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", labels = c("Klasična", "Genomska"), values = c("dashed", "solid")) + xlab("Leto") + 
  geom_ribbon(aes(ymin = value - SD, ymax = value + SD, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=2, keyheight = unit(.8, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  guides(colour=guide_legend(byrow = TRUE, nrow=2, keywidth = unit(1.2, "cm"),override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(variable), cols=vars(scenario), scales = "free_y")+
  theme(strip.text = element_text(size = 16))
G1 <- ggplotGrob(g1)
G1$heights
G1$heights[[10]] <- unit(5, "null")
G1$heights[[12]] <- unit(3, "null")
library(grid)
grid.draw(G1)


#ref == no
gainrefM$value[gainrefM$Ref == "Without initial TP" & gainrefM$Generation == 40 & gainrefM$variable == "TP Size"] <- 0
g2 <- ggplot(data=gainrefM[gainrefM$Ref == "Without initial TP" & gainrefM$scenario == "$G:$P = 1:1",], aes(x=Generation, y=value, group=RealSc, colour=RealSc, linetype = BV)) + 
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16),  axis.title.y = element_blank()) + 
  scale_colour_manual("Scenario", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Scenario", 
                    values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", values = c("dashed", "solid")) + 
  ylab("Genetic mean") +xlab("Year") + 
  geom_hline(data = gainrefM[gainrefM$Ref == "Without initial TP" & gainrefM$variable == "TP Size" & gainrefM$scenario == "$G:$P = 1:1",], aes(yintercept = 2000), colour = "red", size = 1) + 
  geom_ribbon(aes(ymin = value - SD, ymax = value + SD, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(label.position = "top", nrow=2, keywidth = unit(0.7, "cm"), keyheight = unit(0.7, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  guides(colour=guide_legend(byrow = TRUE,label.position = "top", nrow=2, keywidth = unit(1.2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(variable), cols = vars(scenario),  scales = "free_y") +
  theme(strip.text = element_text(size = 16))
#SLO
gainrefM$value[gainrefM$Ref == "Without initial TP" & gainrefM$Generation == 40 & gainrefM$variable == "Velikost RP"] <- 0

gainrefM$variable <- revalue(gainrefM$variable, c("TP Size" = "Velikost RP", "Genetic mean" = "Genetski napredek"))
g2 <- ggplot(data=gainrefM[gainrefM$Ref == "Without initial TP",], aes(x=Generation, y=value, group=RealSc, colour=RealSc, linetype = BV)) +
  geom_line(size = 1) + 
  theme_bw(base_size=18, base_family="sans")  + 
  theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
        axis.text = element_text(size = 16), axis.title.y = element_blank()) + 
  scale_colour_manual("Scenarij", 
                      values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_fill_manual("Scenarij", 
                    values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual("", labels = c("Klasična", "Genomska"), values = c("dashed", "solid")) + xlab("Leto") + 
  geom_ribbon(aes(ymin = value - SD, ymax = value + SD, fill = RealSc),  linetype = 0, alpha = 0.3) + 
  geom_hline(data = gainrefM[gainrefM$Ref == "Without initial TP" & gainrefM$variable == "Velikost RP" ,], aes(yintercept = 2000), colour = "red", size = 1) + 
  
  guides(linetype=guide_legend(nrow=2, keyheight = unit(.8, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  guides(colour=guide_legend(byrow = TRUE, nrow=2, keywidth = unit(1.2, "cm"),override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(variable), cols=vars(scenario), scales = "free_y")+
  theme(strip.text = element_text(size = 16))
G2<- ggplotGrob(g2)
G2$heights
G2$heights[[10]] <- unit(5, "null")
G2$heights[[12]] <- unit(3, "null")
#library(grid)
grid.draw(G2)

#all
ggplot(data=gainrefM, aes(x=Generation, y=value, group=NoControls, colour=NoControls, linetype = BV)) +
    geom_line(size = 1) + 
    theme_bw(base_size=18, base_family="sans")  + 
    theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=16),
          axis.text = element_text(size = 16)) + 
    scale_colour_manual("Number of recordings / lactation", 
                        values = c("#7dcbf5", "#62bff0", "#41b4f0", "#17a3eb", "#076ea3", "#024263", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  ylab("Genetic mean") +
  scale_y_continuous(sec.axis = sec_axis(~ . / 3000)) + 
  #+ geom_ribbon(aes(ymin = minzMean, ymax = maxzMean, fill = NoControl),  linetype = 0, alpha = 0.3) + 
  guides(linetype=guide_legend(nrow=1, keyheight = unit(1.2, "cm"), keywidth = unit(2, "cm"), override.aes = list(alpha = 1, size=1.2))) +
  facet_grid(rows=vars(RefValue), cols=vars(scenario))


#put refsize and accuracy on the same plot
accA <- summarySE(data = acc, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Ref", "Gen"))
accA$Ref <- revalue(accA$Ref, c("True" = "With initial TP", "False" = "Without initial TP"))

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
ggplot(data=accrefM[accrefM$Ref == "Without initial TP",], aes(x=Generation, y=value, group=AgeCat, colour=AgeCat)) +
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
pheno$scenario <- as.factor(pheno$scenario)
pheno$NoControls <- as.factor(pheno$NoControls)

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
          alpha   = 0.05,sort=FALSE,reverse = TRUE,
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
accPlot$gp <- as.factor(accPlot$gp)

#ref = YES
accPlotS <- accPlot[accPlot$Ref == "True" & accPlot$gp == "$G:$P = 2:1",]
accPlotS <- accPlot[accPlot$Ref == "True" & accPlot$gp == "$G:$P = 1:1",]
accPlotS <- accPlot[accPlot$Ref == "True" & accPlot$gp == "$G:$P = 1:2",]
accPlotS <- accPlot[accPlot$Ref == "False" & accPlot$gp == "$G:$P = 2:1",]
accPlotS <- accPlot[accPlot$Ref == "False" & accPlot$gp == "$G:$P = 1:1",]
accPlotS <- accPlot[accPlot$Ref == "False" & accPlot$gp == "$G:$P = 1:2",]
accPlotS <- summarySE(accPlotS, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "Rep"))

accPlotS <- accPlot[accPlot$Ref == "True" & accPlot$AgeCat == "Male candidates",]
accPlotS <- accPlot[accPlot$Ref == "True" & accPlot$AgeCat == "Male candidates",]
accPlotS <- summarySE(accPlotS, measurevar = "Cor", groupvars = c("NoControl", "Rep", "gp"))

#accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:1",]


model_T <- lm(Cor ~ NoControl +AgeCat + AgeCat:NoControl, data=accPlotS) #absolute
marginal = emmeans(model_T, ~ AgeCat:NoControl)
CLD = cld(marginal, by="AgeCat", reverse = TRUE,
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

accPlotS$gp <- factor(accPlotS$gp, levels = c("$G:$P = 2:1","$G:$P = 1:1","$G:$P = 1:2"))
model_T <- lm(Cor ~ NoControl + gp + NoControl:gp, data=accPlotS) #absolute
marginal = emmeans(model_T, ~ NoControl:gp)
CLD = cld(marginal, by="NoControl", 
          alpha   = 0.05, sort = FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#by Ref
accPlotS <- accPlot[accPlot$gp == "$G:$P = 2:1" & accPlot$AgeCat == "Male candidates",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 2:1" & accPlot$AgeCat == "Female candidates",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 2:1" & accPlot$AgeCat == "Mothers",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 2:1" & accPlot$AgeCat == "Fathers",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:1" & accPlot$AgeCat == "Male candidates",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:1" & accPlot$AgeCat == "Female candidates",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:1" & accPlot$AgeCat == "Mothers",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:1" & accPlot$AgeCat == "Fathers",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:2" & accPlot$AgeCat == "Male candidates",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:2" & accPlot$AgeCat == "Female candidates",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:2" & accPlot$AgeCat == "Mothers",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:2" & accPlot$AgeCat == "Fathers",]

accPlotS <- summarySE(accPlotS, measurevar = "Cor", groupvars = c("NoControl", "Rep", "Ref"))

model_T <- lm(Cor ~ NoControl + Ref + NoControl:Ref, data=accPlotS) #absolute
marginal = emmeans(model_T, ~ NoControl:Ref)
CLD = cld(marginal, by="NoControl", 
          alpha   = 0.05, sort = FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD


#ref yes - no
accPlotS <- accPlotA[accPlotA$gp == "$G:$P = 2:1",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:1",]
accPlotS <- accPlot[accPlot$gp == "$G:$P = 1:2",]
accPlotS <- summarySE(accPlotS, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "Ref", "Rep"))

model_T <- lm(Cor ~ NoControl +AgeCat + Ref + AgeCat:NoControl:Ref, data=accPlotA) #absolute
marginal = emmeans(model_T, ~ AgeCat:NoControl:Ref)
CLD = cld(marginal, by="NoControl", reverse = TRUE,
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD


#ref = NO
accPlot <- accPlot[accPlot$Ref == "False",]
model_T <- lm(Cor ~ NoControl + gp + NoControl:gp, data=accPlot) #absolute
marginal = emmeans(model_T, ~ NoControl:gp)
CLD = cld(marginal, by="NoControl",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD

#between gps
accPlotS <- accPlot[accPlot$Ref == "True",]
accPlotS <- summarySE(accPlotS, measurevar = "Cor", groupvars = c("AgeCat", "NoControl", "gp", "Rep"))
model_T <- lm(Cor ~ AgeCat + NoControl + gp + AgeCat:NoControl:gp, data=accPlotS) #absolute
marginal = emmeans(model_T, ~ AgeCat:NoControl:gp)
CLD = cld(marginal, by="AgeCat",
          alpha   = 0.05, sort=FALSE,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "tukey") 
CLD
