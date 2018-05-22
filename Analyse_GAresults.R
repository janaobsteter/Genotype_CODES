acc <- read.csv("~/Documents/PhD/CompBio/TestingGBLUP/REP10_Success/AccuraciesRep.txt", header=TRUE)
library(ggplot2)
library(reshape)
colnames(acc)[1] <- "Rep"

accM <- melt(acc, id.vars = "Rep")
accM$variable <- as.factor(accM$variable)


as.data.frame(colSums(acc) / 10)


ggplot(data=accM, aes(x="Rep", y="value")) +  geom_bar(group="variable", position = "dodge", stat="identity")

# Grouped
ggplot(accM, aes(fill=accM$variable, y=accM$value, x=accM$Rep)) + 
  geom_bar(position="dodge", stat="identity") + xlab("Ponovitev") + ylab("Točnost gPV") + 
  scale_fill_manual(breaks = c("Opt", "Random", "RandomHerd"),"Metoda", 
                      values = c("red", "forestgreen", "blue3"),
                      labels = c("Optimizirana", "Naključna izbira krav", "Naključna izbira čred")) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14))

acc <- read.csv("~/Documents/PhD/CompBio/TestingGBLUP/ACCURACIES.txt", header=TRUE)
acc <- acc[acc$RandomHerd != "RandomHerd",]
colnames(acc) <- c("Rep", "Opt", "Random", "RandomHerd")

rel <- read.csv("~/Documents/PhD/CompBio/TestingGBLUP/RELATIONS.csv", header=TRUE)
rel <- rel[rel$Score != "Score",]
rel$NoAnimals <- as.numeric(as.character(rel$NoAnimals))
rel$NoHerds <- as.numeric(as.character(rel$NoHerds))
aggregate(rel$NoAnimals ~ rel$Way, FUN="mean")
aggregate(rel$NoHerds ~ rel$Way, FUN="mean")

accM <- melt(acc, id.vars = "Rep")
accM$variable <- as.factor(accM$variable)
accM$value <- as.numeric(accM$value)
aggregate(accM$value ~ accM$variable, FUN="mean")
acc$Opt <- as.numeric(as.character(acc$Opt))
acc$Random <- as.numeric(as.character(acc$Random))
acc$RandomHerd <- as.numeric(as.character(acc$RandomHerd))
colSums(acc[,2:4])

as.data.frame(colSums(acc) / 10)

write.table(acc, "~/Accuracies.csv", quote=FALSE)
ggplot(data=accM, aes(x="Rep", y="value")) +  geom_bar(group="variable", position = "dodge", stat="identity")

# Grouped
ggplot(accM, aes(fill=accM$variable, y=accM$value, x=accM$Rep)) + 
  geom_bar(position="dodge", stat="identity") + xlab("Ponovitev") + ylab("Točnost gPV") + 
  scale_x_continuous(breaks =  c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), labels=c(1,2,3,4,5,6,7,8,9,10)) + 
  scale_fill_manual(breaks = c("Opt", "Random", "RandomHerd"),"Metoda", 
                    values = c("red", "forestgreen", "blue3"),
                    labels = c("Optimizirana", "Naključna izbira krav", "Naključna izbira čred")) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14))


