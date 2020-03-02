
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
library(Rmisc)
library(ggplot2)
library(reshape)
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/MSE.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Interaction/MSE.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Intercept/MSE.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/NewModel/MSE.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/NewModel/Reversed/MSE.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/SimpleSimulation//MSE.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Interaction/MSE_noS.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/MSE_newModel.csv", header=TRUE)[-1,]
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/MSE_noS.csv", header=TRUE)[-1,]

table(h2$H2)
h2$Trait <- as.factor(h2$Trait)
levels(h2$Trait)
levels(h2$Trait) <- c("Trait1", "Trait2")
levels(h2$Program) <- c("Program1", "Program2")

mseA <- summarySE(data=h2, groupvars = c("H2", "Program", "Trait", "Path", "Population"), measurevar = "MSE")
mseA <- summarySE(data=h2, groupvars = c("H2", "Program", "Trait", "Path"), measurevar = "MSE")
##both traits same heritabiity
ggplot(data = mseA[mseA$Program == "Program1",], aes(x=H2, y = MSE, group = Path, colour = Path)) + 
  geom_point() + geom_line() + 
  facet_grid(. ~ Program + Trait + Population )
ggplot(data = mseA, aes(x=H2, y = MSE, group = Path, colour = Path)) + 
  geom_point() + geom_line() 


library(reshape)

#READ IN FILES FROM THE CHOSEN MODEL

#####intercept model
part1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Intercept//PartitionPN1.csv", header=TRUE)[-1,]
colnames(part1)[10] <- "BV"
table(part1$h2)
part1 <- part1[,-2]

part2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Intercept//PartitionPN2.csv", header=TRUE)[-1,]
colnames(part2)[10] <- "BV"
part2 <- part2[,-2]
table(part2$h2)
##############################################

#####fixed effects model (newModel), no interaction

part1 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart//NewModel/PartitionPN1.csv", header=TRUE)
colnames(part1)[10] <- "BV"
table(part1$h2)
table(part1$rep)
table(part1$Generation)
part1 <- part1[,-2]

part2 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/NewModel/PartitionPN2.csv", header=TRUE)
colnames(part2)[10] <- "BV"
part2 <- part2[,-2]
table(part2$h2)
table(part2$Generation)
table(part2$rep)
############################################################################################

#####fixed effects model (newModel), no interaction - SAME NUMBERS IN THE MULTIPLIER AND NUCLEUS

part1 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart//NewModel//PartitionPN1.csv", header=TRUE)
colnames(part1)[10] <- "BV"
table(part1$h2)
table(part1$rep)
table(part1$Generation)
part1 <- part1[,-2]

part2 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/NewModel//PartitionPN2.csv", header=TRUE)
colnames(part2)[10] <- "BV"
part2 <- part2[,-2]
table(part2$h2)
table(part2$Generation)
table(part2$rep)
##############################################
'''
#####fixed effects model (newModel), no interaction, REVERSED
part1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/NewModel/Reversed/PartitionPN1.csv", header=TRUE)[-1,]
colnames(part1)[10] <- "BV"
table(part1$h2)
table(part1$rep)
part1 <- part1[,-2]

##############################################

#####interaction model
part1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Interaction/PartitionPN1.csv", header=TRUE)[-1,]
part1 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/Interaction/PartitionPN1.csv", header=TRUE)[-1,]
colnames(part1)[10] <- "BV"
table(part1$h2)
part1 <- part1[,-2]

part2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Interaction/PartitionPN2.csv", header=TRUE)[-1,]
part2 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/Interaction/PartitionPN2.csv", header=TRUE)[-1,]
colnames(part2)[10] <- "BV"
part2 <- part2[,-2]
table(part2$h2)
#####################################################################
#####################################################################
'''


library(reshape)
part1M <- melt(part1, id.vars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "rep"))
#part1M <- melt(part1, id.vars = c("Generation", "Program", "Trait", "BV", "h2", "Population"))
head(part1M)
part1M$BV <- factor(part1M$BV, levels = c("Tbv", "Ebv"))
part1M$Trait <- factor(part1M$Trait, levels = c("T1", "T2", "I"))

part2M <- melt(part2, id.vars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "rep"))
head(part2M)
part2M$BV <- factor(part2M$BV, levels = c("Tbv", "Ebv"))
part2M$Trait <- factor(part2M$Trait, levels = c("T1", "T2", "I"))

library(ggplot2)
# ggplot(data = part1M[part1M$h2 == 0.25,], aes(x=Generation, y = value, group = variable, colour = variable)) + 
#   geom_line() + ggtitle("PN1") + 
#   facet_wrap(. ~ rep + Population + Trait + BV, nrow=3)

# ggplot(data = part1M, aes(x=Generation, y = value, group = variable, colour = variable)) + 
#   geom_line() + ggtitle("PN1") + 
#   facet_wrap(. ~ rep +  BV, nrow=3)

#explore what is going on with male contributions
part1Copy <- part1
part1Copy$Generation[part1Copy$Population == "PN1"] <- part1Copy$Generation[part1Copy$Population == "PN1"] - 1
gen40pn <- part1Copy[part1Copy$Generation == 40 & part1Copy$Trait == "T1" & part1Copy$BV == "Tbv",]

for (rep in 0:9) {
  gen40pn$MaleDiff[gen40pn$rep == rep] <- gen40pn$GN.M[gen40pn$rep == rep & gen40pn$Population == "GN1"] - 
                      gen40pn$GN.M[gen40pn$rep == rep & gen40pn$Population == "PN1"]
  gen40pn$FemaleDiff[gen40pn$rep == rep] <- gen40pn$GN.F[gen40pn$rep == rep & gen40pn$Population == "GN1"] - 
                      gen40pn$GN.F[gen40pn$rep == rep & gen40pn$Population == "PN1"]
}
#to je kul, PN males imajo višje
gen40pn[gen40pn$Population == "PN1",]


part2Copy <- part2
part2Copy$Generation[part2Copy$Population == "PN2"] <- part2Copy$Generation[part2Copy$Population == "PN2"] - 1
gen40pn2 <- part2Copy[part2Copy$Generation == 40 & part2Copy$Trait == "T1" & part2Copy$BV == "Tbv",]

for (rep in 0:9) {
  gen40pn2$MaleDiff[gen40pn2$rep == rep] <- gen40pn2$GN.M[gen40pn2$rep == rep & gen40pn2$Population == "GN2"] - 
    gen40pn2$GN.M[gen40pn2$rep == rep & gen40pn2$Population == "PN2"]
  gen40pn2$FemaleDiff[gen40pn2$rep == rep] <- gen40pn2$GN.F[gen40pn2$rep == rep & gen40pn2$Population == "GN2"] - 
    gen40pn2$GN.F[gen40pn2$rep == rep & gen40pn2$Population == "PN2"]
}
gen40pn2[gen40pn2$Population == "PN2",]

#significance
gen40_t1 <- part1[part1$Generation == 40 & part1$Trait == "T1" & part1$BV == "Tbv",]

gen40_t2 <- part1[part1$Generation == 40 & part1$Trait == "T2" & part1$BV == "Tbv",]

part1Ma <- summarySE(data=part1M, measurevar = "value", groupvars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "variable"))
part2Ma <- summarySE(data=part2M, measurevar = "value", groupvars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "variable"))


#caluclate SD
part1Ma$min <- part1Ma$value - part1Ma$sd
part1Ma$max <- part1Ma$value + part1Ma$sd
part2Ma$min <- part2Ma$value - part2Ma$sd
part2Ma$max <- part2Ma$value + part2Ma$sd



library(ggplot2)
plotFunction <- function (h2) { ggplot(data = part1Ma[part1Ma$h2 == h2,], aes(x=Generation, y = value, group = variable, colour = variable)) + 
  geom_line() + ggtitle("PN2") + ggtitle(h2) +
  facet_grid(. ~ Population + Trait + BV )}
p1PLot <- ggplot(data = part1Ma[part1Ma$h2 == 0.25,], aes(x=Generation, y = value, group = variable, colour = variable)) + 
  geom_line() + ggtitle("PN1") + 
  facet_grid(. ~ Population + Trait + BV )
p2Plot <- ggplot(data = part2Ma[part2Ma$h2 == 0.25,], aes(x=Generation, y = value, group = variable, colour = variable)) + 
  geom_line() +  ggtitle ("PN2") +
  facet_grid(. ~ Population + Trait + BV )

part1Ma$Population <- as.character(part1Ma$Population)
part1Ma$Population[part1Ma$Population == "PN1"] <- "Multiplier"
part1Ma$Population[part1Ma$Population == "GN1"] <- "Nucleus"

part1Ma$variable <- factor(part1Ma$variable, levels = c("Sum", "GN.F", "GN.M", "PN1.F", "PN1.M"))

part2Ma$Population <- as.character(part2Ma$Population)
part2Ma$Population[part2Ma$Population == "PN2"] <- "Multiplier"
part2Ma$Population[part2Ma$Population == "GN2"] <- "Nucleus"

part2Ma$variable <- factor(part1Ma$variable, levels = c("Sum", "GN.F", "GN.M", "PN1.F", "PN1.M"))

part1Ma$Trait <- revalue(part1Ma$Trait, c("T1" = "Trait 1", "T2" = "Trait 2", "I" = "Index"))
part2Ma$Trait <- revalue(part2Ma$Trait, c("T1" = "Trait 1", "T2" = "Trait 2", "I" = "Index"))

part1Ma$Population <- factor(part1Ma$Population, levels = c("Nucleus", "Multiplier"))
part2Ma$Population <- factor(part2Ma$Population, levels = c("Nucleus", "Multiplier"))

tiff("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart//Figures/Obsteter_2.tiff", res=610, width=175, height=100, units="mm")
ggplot(data = part1Ma[(part1Ma$BV == "Tbv") ,],  #$& (part1Ma$variable != "PN1.M")
    aes(x=Generation, y = value, colour = variable, linetype = variable)) + 
    geom_line(size = 1, aes(linetype = variable)) + 
    #ggtitle("Program 1") +  
    geom_ribbon(data = part1Ma[part1Ma$BV == "Tbv",], aes(ymin = min, ymax = max, x = Generation, fill = variable), linetype = 0, alpha = 0.3) + 
    scale_colour_manual("\n\nSelection path", values=c("black", "#bd0b58", "#3ea4ed", "#bd0b58", "#3ea4ed"), 
                        labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) +
    scale_fill_manual("\n\nSelection path", values=c("black", "#bd0b58", "#3ea4ed", "#bd0b58", "#3ea4ed"), 
                      labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) + 
    scale_linetype_manual("\n\nSelection path", values = c("solid", "solid", "solid", "twodash", "twodash"), 
                          labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) +
    ylab("Genetic trend") + 
    theme_bw(base_size=10, base_family="arial") + theme(legend.position="top", legend.text=element_text(size=10), legend.title=element_text(size=12), 
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    strip.text = element_text(size = 10)) + scale_y_continuous(breaks = seq(0, 12, 3)) + 
    guides(colour = guide_legend(keywidth = unit(2.5, "cm"), nrow = 1, byrow = TRUE, label.position =  "top")) +
    facet_grid(. ~ Population + Trait)
dev.off()

tiff("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart//Figures/Obsteter_3.tiff", res=610, width=175, height=100, units="mm")
ggplot(data = part2Ma[(part2Ma$BV == "Tbv"),], 
    aes(x=Generation, y = value, colour = variable, linetype = variable)) + 
    geom_line(size = 1, aes(linetype = variable)) + 
    #ggtitle("Program 2") +  
    geom_ribbon(data = part2Ma[(part2Ma$BV == "Tbv"),], aes(ymin = min, ymax = max, x = Generation, fill = variable), linetype = 0, alpha = 0.3) + 
  scale_colour_manual("\n\nSelection path", values=c("black", "#bd0b58", "#3ea4ed", "#bd0b58", "#3ea4ed"), 
                      labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) +
  scale_fill_manual("\n\nSelection path", values=c("black", "#bd0b58", "#3ea4ed", "#bd0b58", "#3ea4ed"), 
                    labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) + 
  scale_linetype_manual("\n\nSelection path", values = c("solid", "solid", "solid", "twodash", "twodash"), 
                        labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) +
    ylab("Genetic trend") + 
  theme_bw(base_size=10, base_family="Arial") + theme(legend.position="top", legend.text=element_text(size=10), legend.title=element_text(size=12), 
                                                      axis.text=element_text(size=10),
                                                      axis.title=element_text(size=12), strip.text = element_text(size = 10)) + 
    scale_y_continuous(breaks = seq(0, 12, 3)) + 
    guides(colour = guide_legend(keywidth = unit(2.5, "cm"), label.position =  "top")) +
    facet_grid(. ~ Population + Trait)
dev.off()

tiff("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart//Figures/ObsteterLarge_3.tiff", res=610, width=300, height=200, units="mm")
ggplot(data = part2Ma[(part2Ma$BV == "Tbv"),], 
    aes(x=Generation, y = value, colour = variable, linetype = variable)) + 
    geom_line(size = 1, aes(linetype = variable)) + 
    #ggtitle("Program 2") +  
    geom_ribbon(data = part2Ma[(part2Ma$BV == "Tbv"),], aes(ymin = min, ymax = max, x = Generation, fill = variable), linetype = 0, alpha = 0.3) + 
  scale_colour_manual("\n\nSelection path", values=c("black", "#bd0b58", "#3ea4ed", "#bd0b58", "#3ea4ed"), 
                      labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) +
  scale_fill_manual("\n\nSelection path", values=c("black", "#bd0b58", "#3ea4ed", "#bd0b58", "#3ea4ed"), 
                    labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) + 
  scale_linetype_manual("\n\nSelection path", values = c("solid", "solid", "solid", "twodash", "twodash"), 
                        labels = c("Total", "Nucleus\nfemales", "Nucleus\nmales", "Multiplier\nfemales", "Multiplier\nmales")) +
  ylab("Genetic trend") + 
  theme_bw(base_size=18, base_family="Arial") + theme(legend.position="top", legend.text=element_text(size=16), legend.title=element_text(size=18), 
                                                      axis.text=element_text(size=16),
                                                      axis.title=element_text(size=18), strip.text = element_text(size = 16)) + 
    scale_y_continuous(breaks = seq(0, 12, 3)) + 
    guides(colour = guide_legend(keywidth = unit(2.5, "cm"), label.position =  "top")) +
    facet_grid(. ~ Population + Trait)
dev.off()


part1Ma[part1Ma$Generation %in% 40:41 & part1Ma$Trait == "Trait 1" & part1Ma$variable == "Sum" & part1Ma$BV == "Tbv",]
part1Ma[part1Ma$Generation  %in% 40:41 & part1Ma$Trait == "Trait 1" & part1Ma$BV == "Tbv",]
part1Ma[part1Ma$Generation  %in% 40:41 & part1Ma$Trait == "Trait 2" & part1Ma$variable == "Sum" & part1Ma$BV == "Tbv",]
part1Ma[part1Ma$Generation  %in% 40:41 & part1Ma$Trait == "Trait 2"  & part1Ma$BV == "Tbv",]

part2Ma[part2Ma$Generation  %in% 40:41 & part2Ma$Trait == "Trait 1" & part2Ma$variable == "Sum" & part2Ma$BV == "Tbv",]
part2Ma[part2Ma$Generation  %in% 40:41 & part2Ma$Trait == "Trait 1" & part2Ma$BV == "Tbv",]
part2Ma[part2Ma$Generation  %in% 40:41 & part2Ma$Trait == "Trait 2" & part2Ma$variable == "Sum" & part2Ma$BV == "Tbv",]
part2Ma[part2Ma$Generation  %in% 40:41 & part2Ma$Trait == "Trait 2" &  part2Ma$BV == "Tbv",]

part2Ma[part2Ma$Generation == 40 & part2Ma$Trait == "T1",]
part2Ma[part2Ma$Generation == 40 & part2Ma$Trait == "T1",]
part2Ma[part2Ma$Generation == 40 & part2Ma$Trait == "T2",]
part2Ma[part2Ma$Generation == 40 & part2Ma$Trait == "I",]

PART <- data.frame()
for (rep in 0:9) {
  for (generation in 21:41) {
      for (trait in c("T1", "T2")) {
        for (BV in c("Ebv", "Tbv")) {
          for (population in c("GN1", "PN1")) {
            tmp <- part1M[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV & part1M$Population == population & part1M$Generation == generation,]
            tmp$PerSum <- tmp$value / tmp$value[tmp$variable == "Sum"]
            PART <- rbind(PART, tmp)
        }
      }
    }
  }
}

mean(PART$PerSum[PART$Population == "PN1" & PART$variable == "PN1.M" & PART$BV == "Tbv" & PART$Trait == "T1"])
mean(PART$PerSum[PART$Population == "PN1" & PART$variable == "PN1.F" & PART$BV == "Tbv" & PART$Trait == "T1"])
PARTA1 <- summarySE(PART, measurevar = "PerSum", groupvars = c("Program", "Trait", "BV", "Population", "variable", "Generation"))
PARTA1[PARTA1$Generation >= 40 & PARTA1$BV == "Tbv" & PARTA1$Population == "GN1",]
PARTA1[PARTA1$Generation == 41 & PARTA1$BV == "Tbv" & PARTA1$Population == "PN1",]
summary(PARTA1$PerSum[PARTA1$Population == "PN1" & PARTA1$BV == "Tbv" & PARTA1$variable == "PN1.F" & PARTA2$Trait == "T1"])


PARTA1 <- PARTA1[PARTA1$Generation > 23,]
PARTA1 <- summarySE(PARTA1, measurevar = "PerSum", groupvars = c("Program", "Trait", "BV", "Population", "variable"))

PART2 <- data.frame()
  for (rep in 0:9) {
  for (generation in 21:41) {
      for (trait in c("T1", "T2")) {
        for (BV in c("Ebv", "Tbv")) {
          for (population in c("GN2", "PN2")) {
            tmp <- part2M[part2M$rep == rep & part2M$Trait == trait & part2M$BV == BV & part2M$Population == population & part2M$Generation == generation,]
            tmp$PerSum <- tmp$value / tmp$value[tmp$variable == "Sum"]
            PART2 <- rbind(PART2, tmp)
          }
        }
      }
    }
  }

PARTA2 <- summarySE(PART2, measurevar = "PerSum", groupvars = c("Program", "Trait", "BV", "Population", "variable", "Generation"))
summary(PARTA2$PerSum[PARTA2$Population == "PN2" & PARTA2$BV == "Tbv" & PARTA2$variable == "PN2.F" & PARTA2$Trait == "T1"])
summary(PARTA2$PerSum[PARTA2$Population == "PN2" & PARTA2$BV == "Tbv" & PARTA2$variable == "PN2.M" & PARTA2$Trait == "T1"])
PARTA2[PARTA2$Generation >= 40 & PARTA2$BV == "Tbv" & PARTA2$Population == "GN2",]
PARTA2[PARTA2$Generation == 41 & PARTA2$BV == "Tbv" & PARTA2$Population == "PN2",]
PARTA2 <- PARTA2[PARTA2$Generation > 23,]
PARTA2 <- summarySE(PARTA2, measurevar = "PerSum", groupvars = c("Program", "Trait", "BV", "Population", "variable"))
PARTA <- rbind(PARTA1, PARTA2)
PARTA$PerSum <- round(PARTA$PerSum*100, 1)

#ggplot(PARTA[PARTA$BV == "Tbv" & PARTA$Program == "PN1" & PARTA$Population == "PN1",], aes(x=Generation, y = PerSum, group = variable, colour=variable)) + geom_point()

PARTA[PARTA$Program == "PN1" & PARTA$Trait == "T1" & PARTA$BV == "Tbv",]
PARTA[PARTA$Program == "PN1" & PARTA$Trait == "T2" & PARTA$BV == "Tbv",]
PARTA[PARTA$Program == "PN2" & PARTA$Trait == "T1" & PARTA$BV == "Tbv",]
PARTA[PARTA$Program == "PN2" & PARTA$Trait == "T2" & PARTA$BV == "Tbv",]



#naredi razliko med GN and PN v vseh generacijah

part1M$diff <- NA
part1M$diffPer <- NA
for (rep in 0:9) {
  for (generation in 21:40) {
    for (trait in c("T1", "T2")) {
      for (BV in c("Ebv", "Tbv")) {
        part1M$diff[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV &  part1M$Generation == (generation + 1) & part1M$Population == "PN1"]  <- 
        part1M$value[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV &  part1M$Generation == (generation + 1) & part1M$Population == "PN1"] -
        part1M$value[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV &  part1M$Generation == generation & part1M$Population == "GN1"] 
        
        part1M$diffPer[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV &  part1M$Generation == (generation + 1) & part1M$Population == "PN1"]  <- 
        part1M$diff[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV &  part1M$Generation == (generation + 1) & part1M$Population == "PN1"]  /
        part1M$diff[part1M$rep == rep &  part1M$Trait == trait & part1M$BV == BV &  part1M$Generation == (generation + 1) & part1M$Population == "PN1" & part1M$variable == "Sum"]  
      }
    }
  }
}

DIFF <- part1M[part1M$Generation > 23,]
diffA <- summarySE(DIFF[DIFF$Population == "PN1",], measurevar = "diffPer", groupvars = c("Trait", "BV", "Population", "variable", "Generation"), na.rm=TRUE)
diffA <- summarySE(DIFF[DIFF$Population == "PN1",], measurevar = "diffPer", groupvars = c("Trait", "BV", "Population", "variable"), na.rm=TRUE)
diffA[diffA$Generation == 41,]
diffA$diffPer <- round((diffA$diffPer * 100), 1)
diffA

part2M$diff <- NA
part2M$diffPer <- NA
for (rep in 0:9) {
  for (generation in 21:40) {
    for (trait in c("T1", "T2")) {
      for (BV in c("Ebv", "Tbv")) {
        part2M$diff[part2M$rep == rep &  part2M$Trait == trait & part2M$BV == BV &  part2M$Generation == (generation + 1) & part2M$Population == "PN2"]  <- 
        part2M$value[part2M$rep == rep &  part2M$Trait == trait & part2M$BV == BV &  part2M$Generation == (generation + 1) & part2M$Population == "PN2"] -
        part2M$value[part2M$rep == rep &  part2M$Trait == trait & part2M$BV == BV &  part2M$Generation == generation & part2M$Population == "GN2"] 
        
        part2M$diffPer[part2M$rep == rep &  part2M$Trait == trait & part2M$BV == BV &  part2M$Generation == (generation + 1) & part2M$Population == "PN2"]  <- 
        part2M$diff[part2M$rep == rep &  part2M$Trait == trait & part2M$BV == BV &  part2M$Generation == (generation + 1) & part2M$Population == "PN2"]  /
        part2M$diff[part2M$rep == rep &  part2M$Trait == trait & part2M$BV == BV &  part2M$Generation == (generation + 1) & part2M$Population == "PN2" & part2M$variable == "Sum"]  
      }
    }
  }
}

DIFF2 <- part2M[part2M$Generation > 23,]
DIFF2[DIFF2$Generation == 41,]
diffA2 <- summarySE(DIFF2[DIFF2$Population == "PN2",], measurevar = "diffPer", groupvars = c("Trait", "BV", "Population", "variable", "Generation"), na.rm=TRUE)
diffA2[diffA2$Generation == 41,]
diffA2 <- summarySE(DIFF2[DIFF2$Population == "PN2",], measurevar = "diffPer", groupvars = c("Trait", "BV", "Population", "variable"), na.rm=TRUE)
diffA2$diffPer <- round((diffA2$diffPer * 100), 1)
diffA2

 ##################
#stats
part1Ma[part1Ma$Generation == 20 & part1Ma$BV == "Ebv" & part1Ma$Trait == "T1",]
part1Ma[part1Ma$Generation == 20 & part1Ma$BV == "Ebv" & part1Ma$Trait == "T2",]

part20 <- part1Ma[part1Ma$Generation == 20,]
part20$PercentTotal <- NA
for (pop in c("GN", "PN")) {
  for (trait in c("T1", "T2")) {
    part20$PercentTotal[part20$Population == pop & part20$Trait == trait & part20$BV == "Ebv"] <- part20$value[part20$Population == pop & part20$Trait == trait & part20$BV == "Ebv"] / part20$value[part20$Population == pop & part20$Trait == trait & part20$variable == "Sum" & part20$BV == "Ebv"]
  }
}

part20$PercentTotal <- round(part20$PercentTotal, 2)
part20[part20$BV == "Ebv" & part20$Trait == "T1",]
part20[part20$BV == "Ebv" & part20$Trait == "T2",]

for (gen in 1:20) {
  for (pop in c("GN", "PN")) {
    for (trait in c("T1", "T2")) {
      part1Ma$PercentTotal[part1Ma$Population == pop & part1Ma$Trait == trait & part1Ma$BV == "Ebv" & part1Ma$Generation == gen] <- 
        part1Ma$value[part1Ma$Population == pop & part1Ma$Trait == trait & part1Ma$BV == "Ebv" & part1Ma$Generation == gen] / 
        part1Ma$value[part20$Population == pop & part1Ma$Trait == trait & part1Ma$variable == "Sum" & part1Ma$BV == "Ebv" & part1Ma$Generation == gen]
    }
  }
}


part1Ma[part1Ma$PercentTotal == max(part1Ma$PercentTotal[part1Ma$variable == "PN1.F"], na.rm=TRUE) & part1Ma$variable == "PN1.F",]

p2Plot <- ggplot(data = part2Ma, aes(x=Generation, y = value, group = variable, colour = variable)) + 
  geom_line() +  ggtitle ("PN2") +
  facet_grid(. ~ Population + Trait + BV)


library(Rmisc)
multiplot(p1PLot, p2Plot)
multiplot(plotFunction(0.1), plotFunction(0.25),  plotFunction(0.5), plotFunction(0.99))




##################################################################
##################################################################
#simple simulation
acc <- read.csv("~/Documents/Projects/inProgress/AlphaPart/SimpleSimulation/Acc.csv")
mst <- read.csv("~/Documents/Projects/inProgress/AlphaPart/SimpleSimulation/MST.csv")

accA <- summarySE(data = acc, measurevar = "rEBV", groupvars = c("H2", "Generation", "Gender"))
mstA <- summarySE(data = mst, measurevar = "rMST", groupvars = c("H2", "Generation", "Gender"))

ggplot(data=accA, aes(x=Generation, y=rEBV, colour=Gender, group=Gender)) +
  geom_line() +
  facet_grid(. ~ H2)
acc$Group <- paste0(acc$Gender, acc$Rep)
ggplot(data=acc, aes(x=Generation, y=rEBV, colour=Group, group=Group)) +
  geom_line() +
  facet_grid(. ~ H2)

ggplot(data=mstA, aes(x=Generation, y=rMST, colour=Gender, group=Gender)) +
  geom_line() +
  facet_grid(. ~ H2)
mst$Group <- paste0(mst$Gender, mst$Rep)
ggplot(data=mst, aes(x=Generation, y=rMST, colour=Group, group=Group)) +
  geom_line() +
  facet_grid(. ~ H2)
###############################################################
##############################################################

acc <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/Interaction/Accuracies.csv", header=TRUE)[-1,]
acc <- acc[!is.na(acc$Program),]
accA <- summarySE(data=acc, groupvars = c("Program", "Trait", "h2"), measurevar = "Cor")
accA$Trait <- as.factor(accA$Trait)

ggplot(data = accA, aes(x=Trait, y = Cor, group = Program, colour = Program)) + 
  geom_point() + 
  ylim(c(0, 1))

ggplot(data = accA, aes(x=h2, y = Cor, group = Trait, colour = Trait)) + 
  geom_line() + 
  ylim(c(0, 1)) + 
  facet_grid(. ~ Program  )



acc1 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/NewModel/Accuracies_PN1.csv", header=TRUE)[-1,]
acc1A <- summarySE(data=acc1, measurevar = "Cor", groupvars = c("Program", "Trait", "h2", "Generation"))
acc1A$Trait <- as.factor(acc1A$Trait)
ggplot(data = acc1A, aes(x=Generation, y = Cor, group = Trait, colour = Trait)) + 
  geom_line() + ggtitle("PN1") + 
  ylim(c(0, 1)) + 
  facet_grid(. ~ Program + h2  )


acc1A <- summarySE(data=acc1, measurevar = "Cor", groupvars = c("Program", "Trait", "h2"))
acc1A
acc1A$Trait <- as.factor(acc1A$Trait)
acc1Plot <- ggplot(data = acc1A, aes(x=h2, y = Cor, group = Trait, colour = Trait)) + 
  geom_line() + ggtitle("PN1") + 
  ylim(c(0, 1)) + 
  facet_grid(. ~ Program  )

acc2 <- read.table("~/Documents/PhD/Projects/inProgress/AlphaPart/NewModel/Accuracies_PN2.csv", header=TRUE)[-1,]
acc2A <- summarySE(data=acc2, measurevar = "Cor", groupvars = c("Program", "Trait", "h2", "Generation"))
acc2A <- summarySE(data=acc2, measurevar = "Cor", groupvars = c("Program", "Trait", "h2"))
acc2A
acc2A$Trait <- as.factor(acc2A$Trait)
ggplot(data = acc2A, aes(x=Generation, y = Cor, group = Trait, colour = Trait)) + 
  geom_line() + ggtitle("PN2") + 
  ylim(c(0, 1)) + 
  facet_grid(. ~ Program + h2  )
acc2Plot <- ggplot(data = acc2A, aes(x=h2, y = mean, group = Trait, colour = Trait)) + 
  geom_line() + ggtitle("PN2") + 
  ylim(c(0, 1)) + 
  facet_grid(. ~ Program  )


multiplot(acc1Plot, acc2Plot)

##accuracy of partitiones EBVs
partAcc <- read.table("~/Documents/Projects/inProgress/AlphaPart/AccuracyPart.csv", header=TRUE)[-1,]

plot1 <- ggplot(data=partAcc[partAcc$Program == "PN1", ], aes(x=H2, y=Cor, group=Path, colour=Path)) + 
  geom_line() + 
  facet_grid(. ~ Trait)
plot2 <- ggplot(data=partAcc[partAcc$Program == "PN2", ], aes(x=H2, y=Cor, group=Path, colour=Path)) + 
  geom_line() + 
  facet_grid(. ~ Trait)

##post corelation
postCor <- read.table("~/Documents/Projects/inProgress/AlphaPart/PostCorrelation.csv", header=TRUE)

ggplot(data=postCor, aes(x=h2, y=COR, grouop=ProgramGender, colour=ProgramGender)) + 
  geom_line() + 
  facet_grid(. ~ Program + Trait)
  
  


#####standardization factors
stF <- read.table("Documents/Projects/inProgress/AlphaPart/Interaction/StandFactors.csv", header=TRUE)[-1,]
colnames(stF)[1] <- "rep"

ggplot(data=stF[stF$Program == "Program1" & stF$rep == 1,], aes(x=Variable, y=Value, group=BV, colour=BV)) + geom_point() + 
  facet_grid(.~h2 + Trait)

stF <- spread(stF, key="Variable", value="Value")
levels(stF$Trait) <-  c("T1", "T2")

stMean <- stF[,1:6]
levels(stMean$Program) <- c("PN1", "PN2")

#calculate difference between EBV and TBV Mean

stFF <- spread(data = stMean, key = "BV", value= "Mean")
stFF$DIFF <- stFF$Tbv - stFF$Ebv  
part1SD <- merge(part1[part1$Trait != "I" & part1$BV == "Ebv",], stFF, by=c("rep", "h2", "Program", "Trait"))
head(part1SD)




#distribution of EBVs
PedEval1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/SimpleSimulation/PedEval_0.99.csv", header=TRUE)
PedEval1 <- read.table("~/PedEval1_0.5_NewModel.csv", header=TRUE)
PedEval1 <- read.table("~/PedEval1_0.5_Old.csv", header=TRUE)
PedEval1 <- read.table("~/PedEval1_0.25.csv", header=TRUE)


gen0 <- PedEval1[PedEval1$Generation == 0,]
plot(density(gen0$TbvT1))
sd(gen0$TbvT1)
mean(gen0$TbvT1)
plot(density(gen0$EbvT1))
sd(gen0$EbvT1)
mean(gen0$EbvT1)

tbv1 <- ggplot(data=gen0, aes(x=TbvT1)) + geom_density() + ggtitle("TbvT1")
ebv1 <- ggplot(data=gen0, aes(x=EbvT1)) + geom_density() + ggtitle("EbvT1")
tbv2 <- ggplot(data=gen0, aes(x=TbvT2)) + geom_density() + ggtitle("TbvT2")
ebv2 <- ggplot(data=gen0, aes(x=EbvT2)) + geom_density() + ggtitle("EbvT2")
multiplot(tbv1, ebv1, tbv2, ebv2, cols=2)
plot(density(gen0$TbvT2))
sd(gen0$TbvT2)
mean(gen0$TbvT2)
plot(density(gen0$EbvT2))
sd(gen0$EbvT2)
mean(gen0$EbvT2)

#PedEval1 <- read.table("~/PedEval1_0.05.csv", header=TRUE)
#EBV
PedEval1$EbvT1_s <- (PedEval1$EbvT1 - mean(PedEval1$EbvT1[PedEval1$Generation == 0])) / sd(PedEval1$EbvT1[PedEval1$Generation == 0])
PedEval1$EbvT2_s <- (PedEval1$EbvT2 - mean(PedEval1$EbvT2[PedEval1$Generation == 0])) / sd(PedEval1$EbvT2[PedEval1$Generation == 0])
PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)
PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")

PedEval1$GenerationProgram <- paste(PedEval1$Generation, PedEval1$Program, sep="-")



# ---- Partitioning the trend PN1 ----
PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
PedEval1$TbvT1_s <- (PedEval1$TbvT1 - mean(PedEval1$TbvT1[PedEval1$Generation == 0])) / sd(PedEval1$TbvT1[PedEval1$Generation == 0])
PedEval1$TbvT2_s <- (PedEval1$TbvT2 - mean(PedEval1$TbvT2[PedEval1$Generation == 0])) / sd(PedEval1$TbvT2[PedEval1$Generation == 0])
PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)

PedEval1$PhenoI = 0.5 * (PedEval1$PhenoT1 + PedEval1$PhenoT2)
PedEval1$PhenoT1_s <- (PedEval1$PhenoT1 - mean(PedEval1$PhenoT1[PedEval1$Generation == 0])) / sd(PedEval1$PhenoT1[PedEval1$Generation == 0])
PedEval1$PhenoT2_s <- (PedEval1$PhenoT2 - mean(PedEval1$PhenoT2[PedEval1$Generation == 0])) / sd(PedEval1$PhenoT2[PedEval1$Generation == 0])
PedEval1$PhenoI_s <- 0.5 * (PedEval1$PhenoT1_s + PedEval1$PhenoT2_s)



PedEval1 <- read.table("/home/jana/PedEval1_0.25.csv", header=TRUE)
PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")

library(AlphaPart)
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1", "EbvT2"))

Part1g = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                   colId = "IId", colFid = "FId", colMid = "MId",
                   colPath = "ProgramGender", colAGV = c("TbvT1", "TbvT2"))

Part1Summary = summary(object = Part1, by = "Generation")
Part1gSummary = summary(object = Part1g, by = "Generation")


p1 <- plot(Part1Summary)
p1g <- plot(Part1gSummary)

head(PedEval)
gainAe <- summarySE(data=PedEval1, measurevar = "EbvT1_s", groupvars = c("Generation", "Program"))
gainAe$BV <- "EBV"
colnames(gainAe)[4] <- "value"
gainAt <- summarySE(data=PedEval1, measurevar = "TbvT1_s", groupvars = c("Generation", "Program"))
gainAt$BV <- "TBV"
colnames(gainAt)[4] <- "value"
gainAp <- summarySE(data=PedEval1, measurevar = "PhenoT1_s", groupvars = c("Generation", "Program"))
gainAp$BV <- "Pheno"
colnames(gainAp)[4] <- "value"
gainA <- rbind(gainAe, gainAt)
gainA <- rbind(gainA, gainAp)
gainA$Generation <- as.factor(gainA$Generation)
ggplot(data = gainA, aes(x=Generation, y=value, group=BV, colour=BV)) + geom_line() + facet_grid(.  ~ Program)





head(gainA)
gainAm <- melt(gainA, id.vars = c("Generation", "Program", "BV"), measure.vars = c("value", "sd"))

ggplot(data = gainAm, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~ Program + BV)


PedEval1[PedEval1$Generation == 1,]




##trends
de1 <- read.csv("~/Documents/PhD/Projects/inProgress/AlphaPart/NewModel/GeneticTrendsPN1.csv")
de2 <- read.csv("~/Documents/PhD/Projects/inProgress/AlphaPart/NewModel/GeneticTrendsPN2.csv")
#de1 <- read.csv("~/Documents/PhD/Projects/inProgress/AlphaPart/Interaction/GeneticTrendsPN1.csv")
#de2 <- read.csv("~/Documents/PhD/Projects/inProgress/AlphaPart/Interaction/GeneticTrendsPN2.csv")
#de1 <- read.csv("~/Documents/PhD/Projects/inProgress/AlphaPart/Intercept/GeneticTrendsPN1.csv")
#de2 <- read.csv("~/Documents/Projects/inProgress/AlphaPart/Intercept/GeneticTrendsPN2.csv")
#de1 <- read.csv("~/Documents/Projects/inProgress/AlphaPart/NewModel/Reversed/GeneticTrendsPN1.csv")
head(de1)
table(de1$h2)

de1A <- summarySE(data=de1, groupvars = c("Generation", "Program", "trait"), measurevar = "mean")
de1A$P <- "P1"
#de1A <- summarySE(data=de1, groupvars = c("Generation",  "trait", "BV", "h2"), measurevar = "mean")
de2A <- summarySE(data=de2, groupvars = c("Generation", "Program", "trait"), measurevar = "mean")
de2A$P <- "P2"
deA <- rbind(de1A, de2A)
deA[deA$Generation == 40,]
deA[deA$Generation == 21,]

ggplot(dat = de1A[de1A$trait == "T2",], aes(x=Generation, y=mean, colour = trait, group=trait))  + geom_line() + 
  facet_grid(. ~ Program)
ggplot(dat = de2A, aes(x=Generation, y=mean, colour = trait, group=trait))  + geom_line() + 
facet_grid(. ~ Program)

ggplot(dat = deA, aes(x=Generation, y=mean, colour = trait, group=trait))  + geom_line() + 
  theme_bw(base_size = 18) + 
  facet_grid(. ~ P + Program )

ggplot(dat = de1A, aes(x=Generation, y=mean, colour = BV, group=BV))  + geom_line() + 
  theme_bw(base_size = 18) + 
  facet_grid(. ~ Program + h2 +  trait)

ggplot(dat = deA, aes(x=Generation, y=mean, colour = BV, group=BV))  + geom_line() + 
  facet_grid(. ~ Program + h2 +  trait) 

##correlations
corr <- read.csv("~/Documents/Projects/inProgress/AlphaPart/NewModel/CorrelationBVs.csv")[-1,]
ggplot(data=corr, aes(x=H2, y=Cor, group=Trait, colour=Trait)) + geom_point() +
  facet_grid(. ~ BVCor + Program)
library(ggplot2)



fathers <-  unique(PedEval1$FId)
mothers <-  unique(PedEval1$MId)
pedFather <- PedEval1[PedEval1$IId %in% fathers,c("IId", "TbvT1")]
pedMother <- PedEval1[PedEval1$IId %in% mothers,c("IId", "TbvT1")]
colnames(pedFather) <- c("FId", "gvF")
colnames(pedMother) <- c("MId", "gvM")

ped1 <- merge(PedEval1, pedFather, by="FId", all.x=TRUE)
nrow(ped1)
ped1 <- merge(ped1, pedMother, by="MId", all.x=TRUE)
nrow(ped1)
ped1$MST <- ped1$TbvT1 -  ((ped1$gvF +  ped1$gvM) / 2)


fathers <-  unique(PedEval1$FId)
mothers <-  unique(PedEval1$MId)
pedFather <- PedEval1[PedEval1$IId %in% fathers,c("IId", "EbvT1")]
pedMother <- PedEval1[PedEval1$IId %in% mothers,c("IId", "EbvT1")]
colnames(pedFather) <- c("FId", "ebvF")
colnames(pedMother) <- c("MId", "ebvM")

ped1 <- merge(ped1, pedFather, by="FId", all.x=TRUE)
nrow(ped1)
ped1 <- merge(ped1, pedMother, by="MId", all.x=TRUE)
nrow(ped1)
ped1$MSTe <- ped1$EbvT1 -  ((ped1$ebvF +  ped1$ebvM) / 2)



as.data.frame(ped1 %>% 
  group_by(Generation, Program, Gender) %>%
  summarize(COR=cor(MST, MSTe)))


ped <- read.table("/home/jana/PedEval1_0.25.csv", header=TRUE)
ped <- ped[ped$Generation > 20,]
library(kinship2)
pedig <- fixParents(ped$IId, ped$FId, ped$MId, ped$Gender)
pedig <- pedig[1:100,]
Ped <- fixParents(pedig$id, pedig$dadid, pedig$momid, pedig$sex)
pedig <- Ped
Ped <- pedigree(pedig$id, pedig$dadid, pedig$momid, pedig$sex)
plot(Ped)

install.packages("pedantics")
library(pedantics)
colnames(pedig) <- c("id", "dam", "sire")
pedig <- ped[,1:3]
colnames(pedig) <- c("id", "dam", "sire")
drawPedigree(pedig[,1:3])


#distribution of true breeding values
ped2 <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart/NewModel//PedEval2_0.25.csv", header=TRUE, sep=" ")
ped2$P <- "Program 2"
ped1 <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart/NewModel//PedEval1_0.25.csv", header=TRUE, sep=" ")
ped1$P <- "Program 1"
##check the number of fathers
gen40_2 <- ped2[ped2$Generation == 40,]
length(unique(gen40_2$FId))
length(unique(gen40_2$MId))
father <- ped2[ped2$IId %in% gen40_2$FId,]
table(father$Program)

ped2$TbvT1_s <- (ped2$TbvT1 - mean(ped2$TbvT1[ped2$Generation == 20])) / sd(ped2$TbvT1[ped2$Generation == 20])
ped2$TbvT2_s <- (ped2$TbvT2 - mean(ped2$TbvT2[ped2$Generation == 20])) / sd(ped2$TbvT2[ped2$Generation == 20])
ped1$TbvT1_s <- (ped1$TbvT1 - mean(ped1$TbvT1[ped1$Generation == 20])) / sd(ped1$TbvT1[ped1$Generation == 20])
ped1$TbvT2_s <- (ped1$TbvT2 - mean(ped1$TbvT2[ped1$Generation == 20])) / sd(ped1$TbvT2[ped1$Generation == 20])

ped <- rbind(ped1, ped2)
table(ped$P)
table(ped$Program)

table(ped$Program)
pedL <- rbind(ped[ped$Generation == 40 & ped$Program == "GN",],
              ped[ped$Generation == 41 & ped$Program == "PN1",],
              ped[ped$Generation == 41 & ped$Program == "PN2",])
pedL$TbvI_s <- (pedL$TbvT1_s + pedL$TbvT2_s) / 2
pedL$Program <- as.character(pedL$Program)
pedL$Program[pedL$Program == "GN" & pedL$P == "Program 2"] <- "GN2"
table(pedL$Program)
pedL <- pedL[,c("Program", "P", "TbvT1_s", "TbvT2_s", "TbvI_s")]
pedLm <- melt(pedL, id.vars = c("Program", "P"))
head(pedLm)
pedLm$variable <- revalue(pedLm$variable, c("TbvT1_s" = "Trait 1", "TbvT2_s" = "Trait 2", "TbvI_s" = "Index"))
pedLm$Program <- revalue(pedLm$Program, c("GN1" = "GN", "PN1" = "PN", "GN2" = "GN", "PN2" = "PN"))
pedLm$ProgramPaper <- ifelse(pedLm$P == "Program 1", "MaleFlow100", "MaleFlow52")
#če imaš samo en program
#tiff("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart//Figures/Obsteter_1.tiff", res=610, width=175, height=100, units="mm")
#če imaš oba programa (dve vrsti)
tiff("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart//Figures/Obsteter_1.tiff", res=610, width=175, height=120, units="mm")
ggplot(pedLm, aes(x = value, group = Program, linetype = Program)) + geom_density() +  #, ((..count..)/sum(..count..)))
  facet_grid(rows = vars(ProgramPaper), cols = vars(variable)) + 
  ylab("Density") + 
  xlab("True breeding value") + 
  scale_linetype_manual("Tier", labels = c("Nucleus", "Multiplier"), values = c("solid", "dashed")) + 
  theme_bw(base_size=10, base_family="arial") + 
  theme(legend.position="top", legend.text=element_text(size=10), legend.title=element_text(size=12), 
                                                      axis.text=element_text(size=10),
                                                      axis.title=element_text(size=12),
                                                      strip.text = element_text(size = 10))
dev.off()
table(pedLm$Program)

ggplot(pedLm[pedLm$variable == "Index",], aes(x = value, group = Program, linetype = Program)) + geom_freqpoly() 
  

