ref <- data.frame()

for (start in c(seq(0, 4545, 505))) {
	tmp <- read.csv(paste0("RefADF_mean",  start, ".csv"))
	ref <- rbind(ref, tmp)
}

write.csv(ref, "RefADF_mean.csv", quote=FALSE, row.names=FALSE)
