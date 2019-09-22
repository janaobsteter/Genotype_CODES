acc <- read.table("~/Documents/Projects/inProgress/Phenotyping/CombinedAcc.csv", header=TRUE)
acc$Cor <- as.numeric(acc$Cor)
acc <- acc[acc$Gen %in% 40:60,]
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
ggplot(data=accA, aes(x=Ref, y=Cor, group=varE, colour=varE)) + geom_point() + 
  scale_colour_manual("Scenario", labels=c("c", "e", "d", "f", "a","b"),  values = c("black", "red", "blue", "orange", "purple", "forestgreen")) +
  facet_grid(. ~ AgeCat) 

ggplot(data=acc, aes(x=Ref, y=Cor, group=varE, colour=varE)) + geom_point() + 
  scale_colour_manual("Scenario", labels=c("c", "e", "d", "f", "a","b"), values = c("black", "red", "blue", "orange", "purple", "forestgreen")) +                            
  facet_grid(. ~ AgeCat) 