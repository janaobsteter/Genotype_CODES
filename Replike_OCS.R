ped <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results//TGVsAll_OCS_19032019.csv", sep=" ")


ggplot(data=ped, aes(x=Generation, y=zMean, group=degree, colour=degree)) + geom_path()
