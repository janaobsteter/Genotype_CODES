Cat <- data.frame()
CatAge <- data.frame()

for (rep in c(1,2,5,8,9,11)) {
	cat <- read.csv(paste0("Accuracy_Cat", rep, ".csv"))
	catage <- read.csv(paste0("Accuracy_CatAge", rep, ".csv"))
	Cat <- rbind(Cat, cat)
	CatAge <- rbind(CatAge, catage)
}

write.csv(Cat, "Accuracy_Cat_permEnv.csv", quote=FALSE, row.names=FALSE)
write.csv(CatAge, "Accuracy_CatAge_permEnv.csv", quote=FALSE,	row.names=FALSE)

