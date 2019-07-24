ped <- read.table("PedImport100.txt", header=TRUE)
head(ped)

popS <- read.csv("PopSplit100.txt")
colnames(popS) <- c("Group", "Indiv")

ped <- merge(ped, popS, by="Indiv")
ped$import <- 100

#for 75
ped75 <- read.table("PedImport75.txt", header=TRUE)
head(ped)

popS75 <- read.csv("PopSplit75.txt")
colnames(popS75) <- c("Group", "Indiv")

ped75 <- merge(ped75, popS75, by="Indiv")
ped75$import <- 75

ped <- rbind(ped, ped75)

##############
pedA <- summarySE(data=ped, measurevar = "gvNormUnres1", groupvar = c("Generation", "Group", "import")  )
head(pedA)
pedA$PlotGroup <- paste0(pedA$Group, pedA$import)

library(ggplot2)
ggplot(data=pedA, aes(x=Generation, y=gvNormUnres1, group=PlotGroup, colour=PlotGroup)) + geom_line()

#standardise by group
PEDG <- data.frame()
pedA <- pedA[pedA$Generation %in% 40:60,]

for (group in c("home", "import")) {
  for (imp in unique(pedA$import)) {
  pedG <- pedA[pedA$Group == group & pedA$import == imp  ,]
  pedG$zMean <- (pedG$gvNormUnres1 - pedG$gvNormUnres1[pedG == 40]) / pedG$sd[pedG == 40]
  PEDG <- rbind(PEDG, pedG)
  }
}

PEDG$PlotGroup <- paste0(PEDG$Group, PEDG$import)


head(PEDG)
ggplot(data=PEDG, aes(x=Generation, y=zMean, group=PlotGroup, colour=PlotGroup)) + geom_line()
