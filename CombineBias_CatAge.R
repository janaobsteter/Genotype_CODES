Cat <- data.frame()
CatAge <- data.frame()

for (strategy in c("SU55", "SU51", "SU15")) {
	cat <- read.csv(paste0("Bias_Cat", strategy, ".csv"))
	catage <- read.csv(paste0("Bias_CatAge", strategy, ".csv"))
	Cat <- rbind(Cat, cat)
	CatAge <- rbind(CatAge, catage)
}

write.csv(Cat, "Bias_Cat.csv", quote=FALSE, row.names=FALSE)
write.csv(CatAge, "Bias_CatAge.csv", quote=FALSE,	row.names=FALSE)

