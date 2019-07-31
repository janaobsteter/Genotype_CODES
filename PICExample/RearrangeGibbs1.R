args = commandArgs(trailingOnly=TRUE)

gibbsFile = args[1]
sampleNum = args[2]


g <- read.table(gibbsFile)
g$samples <- rep(1:sampleNum, length(unique(g$V2)))
gS <- spread(g, samples, V3)
write.csv(gS, "Gibbs1.csv", row.names=FALSE, quote=FALSE)

