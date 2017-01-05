#read genotyped individuals

GenInd <- read.csv ("GenotypedInd", sep=" ", header=F)
colnames(GenInd) <- c("ID_ZIVALI", "GEN_CHIP", "GEN_DATUM")
seq <- read.csv(IDSeqFile, sep=",", header=T)
GenInd1 <- merge(GenInd, seq, by="ID_ZIVALI")
write.table(GenInd1, "Genotyped_individuals.csv", row.names=F, quote=F, sep=",", na="")
