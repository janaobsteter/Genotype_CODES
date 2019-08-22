library(Rmisc)
library(ggplot2)
library(reshape)
h2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/MSE.csv", header=TRUE)[-1,]

h2$Trait <- as.factor(h2$Trait)
levels(h2$Trait)
levels(h2$Trait) <- c("Trait1", "Trait2")
levels(h2$Program) <- c("Program1", "Program2")


##both traits same heritabiity
ggplot(data = h2[h2$Program == "Program1",], aes(x=H2, y = MSE, group = Path, colour = Path)) + 
  geom_point() + geom_line() + 
  facet_grid(. ~ Program + Trait + Population )


library(reshape)
part1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/PartitionPN1.csv", header=TRUE)[-1,]
#part1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/PartitionPN1_oneRenum.csv", header=TRUE)[-1,]
colnames(part1)[10] <- "BV"
part1 <- part1[,-2]
part2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/PartitionPN2.csv", header=TRUE)[-1,]
colnames(part2)[10] <- "BV"
part2 <- part2[,-2]


part1M <- melt(part1, id.vars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "rep"))
#part1M <- melt(part1, id.vars = c("Generation", "Program", "Trait", "BV", "h2", "Population"))
head(part1M)
part1M$BV <- factor(part1M$BV, levels = c("Tbv", "Ebv"))
part1M$Trait <- factor(part1M$Trait, levels = c("T1", "T2", "I"))
part2M <- melt(part2, id.vars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "rep"))
head(part2M)
part2M$BV <- factor(part2M$BV, levels = c("Tbv", "Ebv"))
part2M$Trait <- factor(part2M$Trait, levels = c("T1", "T2", "I"))

part1Ma <- summarySE(data=part1M, measurevar = "value", groupvars = c("Generation", "Program", "Trait", "BV", "h2", "Population", "variable"))
library(ggplot2)
p1PLot <- ggplot(data = part1M[part1M$h2 == 0.5,], aes(x=Generation, y = value, group = variable, colour = variable)) + 
  geom_line() + ggtitle("PN1") + 
  facet_grid(. ~ Population + Trait + BV )
p2Plot <- ggplot(data = part2M[part2M$h2 == 0.5,], aes(x=Generation, y = value, group = variable, colour = variable)) + 
  geom_line() +  ggtitle ("PN2") +
  facet_grid(. ~ Population + Trait + BV )

library(Rmisc)
multiplot(p1PLot, p2Plot)
multiplot(p1PLot, p1PLotOne)


p1PLotOne <- p1PLot
#multiplot(p1PLot, p1PLOT)
##accuracy
acc1 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Accuracies_PN1.csv", header=TRUE)[-1,]
acc1A <- summarySE(data=acc1, measurevar = "Cor", groupvars = c("Program", "Trait", "h2"))

acc1A$Trait <- as.factor(acc1A$Trait)
acc1Plot <- ggplot(data = acc1A, aes(x=h2, y = mean, group = Trait, colour = Trait)) + 
  geom_line() + ggtitle("PN1") + 
  ylim(c(0, 1)) + 
  facet_grid(. ~ Program  )

acc2 <- read.table("~/Documents/Projects/inProgress/AlphaPart/Accuracies_PN2.csv", header=TRUE)[-1,]
acc2A <- summarySE(data=acc2, measurevar = "Cor", groupvars = c("Program", "Trait", "h2"))

acc2A$Trait <- as.factor(acc2A$Trait)
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
  
  
