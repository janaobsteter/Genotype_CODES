library(ggplot2)
# acc <- read.table("~/Documents/Projects/inProgress/Phenotyping/Correlation_reference.csv", header=TRUE)[-1,]
ggplot(data=acc, aes(x=nInd, y=Cor)) + geom_point() + scale_y_continuous(breaks = c(seq(0, 0.9, by=0.1)))

acc$size <- gsub("K", "", acc$Ref)
acc$size <- as.numeric(acc$size )* 1000
accA <- aggregate(acc$Cor ~ acc$size + acc$Pop, FUN="mean")
colnames(accA) <- c("Ref", "Pop", "Cor")

ggplot(data=accA[accA$Pop != "All",], aes(x=Ref, y=Cor, group=Pop, colour=Pop)) + geom_point() ## + geom_smooth()





Quantile1 <- read.csv("CombinedQuantile_T1.csv")[-1,]
Quantile2 <- read.csv("CombinedQuantile_T2.csv")[-1,]
Quantile <- rbind(Quantile1, Quantile2)

part11 <- read.csv("/home/jana/CombinedPartitionALL_PN1_T1.csv")
part11$Trait <- "T1"
part11$Program <- "PN1"
ggplot(data=part, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line()

part12 <- read.csv("/home/jana/CombinedPartition_PN1_T2.csv")
part12$Trait <- "T2"
part12$Program <- "PN1"

ggplot(data=part12, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line()

part21 <- read.csv("/home/jana/CombinedPartition_PN2_T1.csv")
part21$Trait <- "T1"
part21$Program <- "PN2"

ggplot(data=part21, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line()

part22 <- read.csv("/home/jana/CombinedPartition_PN2_T2.csv")
part22$Trait <- "T2"
part22$Program <- "PN2"

ggplot(data=part22, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_line()

PART <- rbind(part11, part12)
PART <- rbind(PART, part21)
PART <- rbind(PART, part22)
colnames(PART)[2] <- "Pop"
PART <- merge(PART, Quantile, by=c("Generation", "Pop", "Trait", "Program"))
ggplot(data=PART[PART$Program == "PN1",], aes(x=Generation, y=value, group=Pop, colour=Pop)) + 
  geom_line() +
  geom_ribbon(aes(ymin=X2.5., ymax=X97.5.,x=Generation, fill=Pop),  alpha=0.2, colour=NA) + 
  facet_grid(. ~ Program + Trait)
