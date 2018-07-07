acc <- read.csv("~/Documents/PhD/CompBio/TestingGBLUP//ACCURACIES.txt", header=TRUE)
library(ggplot2)
library(reshape)
colnames(acc)[1] <- "Rep"

accM <- melt(acc, id.vars = "Rep")
accM$variable <- as.factor(accM$variable)


as.data.frame(colSums(acc) / 10)


ggplot(data=accM, aes(x="Rep", y="value")) +  geom_bar(group="variable", position = "dodge", stat="identity")
accM$value <- as.numeric(as.character(accM$value))
# Grouped
ggplot(accM, aes(fill=accM$variable, y=accM$value, x=accM$Rep)) + 
  geom_bar(position="dodge", stat="identity") + xlab("Ponovitev") + ylab("Točnost gPV") + 
  scale_fill_manual(breaks = c("Opt", "Random", "RandomHerd"),"Metoda", 
                      values = c("red", "forestgreen", "blue3"),
                      labels = c("Optimizirana", "Naključna izbira krav", "Naključna izbira čred")) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14))


# Random as a line
acc$Opt <- as.numeric(as.character(acc$Opt))
acc$Random <- as.numeric(as.character(acc$Random))
acc$RandomHerd <- as.numeric(as.character(acc$RandomHerd))
ggplot(data = acc, aes(x=acc$Rep, y=acc$Opt)) + 
  geom_bar(position="dodge", stat="identity", fill="orangered3") + xlab("Ponovitev") + ylab("Točnost gPV") + 
  geom_text(aes(label=round(acc$Opt,2)), vjust=-0.3, size=5) +
  geom_hline(yintercept = mean(acc$Random), color="forestgreen", size=2) + 
  geom_hline(yintercept = mean(acc$RandomHerd), color="blue", size=2) + 
  geom_text(aes(0,mean(acc$Random),label = round(mean(acc$Random),2), vjust = -1, hjust=+2), color="forestgreen", size=5) +
  geom_text(aes(0,mean(acc$RandomHerd),label = round(mean(acc$RandomHerd),2), vjust = -1, hjust=+2), color="blue", size=5) +
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


