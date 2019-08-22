males <- read.table("~/Documents/Projects/inProgress/Phenotyping/TGVsAll_permEnv_MaleSU55_16082019.csv", header=TRUE)

library(ggplot2)
ggplot(males, aes(x=Generation, y=zMean, group=scenario, colour=scenario)) + geom_line()



acc <- read.table("~/Documents/Projects/inProgress/Phenotyping/Correlation_reference.csv", header=TRUE)[-1,]
acc$size <- gsub("K", "", acc$Ref)
acc$size <- as.numeric(acc$size )* 1000
accA <- aggregate(acc$Cor ~ acc$size + acc$Pop, FUN="mean")
colnames(accA) <- c("Ref", "Pop", "Cor")

ggplot(data=accA[accA$Pop != "All",], aes(x=Ref, y=Cor, group=Pop, colour=Pop)) + geom_line()





acc <- read.csv("~/Documents/Projects/inProgress/Phenotyping/AccuracyMale.csv", header=TRUE)
acc <- acc[acc$AgeCat %in% c("telF1", "genTest1", "potomciNP0"),]
ggplot(data=acc, aes(x=AgeCat, y=Cor, group=scenario, colour=scenario)) + geom_point()


ped25 <- 