ped <- read.csv("/home/jana/TGVsAll_permEnv_SU55_28072019.csv", sep=" ")

head(ped)

library(ggplot2)
ped$PlotGroup <- paste0(ped$scenario, ped$Repeats, "_", ped$Rep)
ped$PlotGroup1 <- paste0(ped$scenario, ped$Repeats)
head(ped)
ped$Repeats <- as.factor(ped$Repeats)
ggplot(data = ped, aes(x=Generation, y=zMean, colour=Repeats, linetype = scenario, group = PlotGroup)) + geom_line()

pedA <- aggregate(ped$zMean ~ ped$Generation + ped$scenario + ped$Repeats, FUN="mean")
head(pedA)
colnames(pedA) <- c("Generation", "Scenario", "Repeats", "TGV")
pedA$PlotGroup <- paste0(pedA$Scenario, pedA$Repeats)

ggplot(data=pedA, aes(x=Generation, y = TGV, group=PlotGroup, colour=Repeats, linetype=Scenario)) + geom_line()
