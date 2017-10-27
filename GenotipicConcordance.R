WorkingDir = "/home/jana/Genotipi/Genotipi_WORK/ImputationOntoIDBv03"
imputed <- read.table("/home/jana/Genotipi/Genotipi_WORK/ImputationOntoIDBv03/IMPUTEDRAWFILE.raw", header=TRUE)
imputed <- imputed[,-c(1,3,4,5,6)]

orig <- read.table("/home/jana/Genotipi/Genotipi_WORK/ImputationOntoIDBv03/ORIGINALRAWFILE.raw", header=TRUE)
orig <- orig[,-c(1,3,4,5,6)]

cor <- c()
for (row in 1:nrow(imputed)) {
  imp <- as.data.frame(t(imputed[row,]))
  imp$snps <- row.names(imp)
  sample <- as.character(imp[1,1])
  colnames(imp) <- c(as.character(imp[1,1]), "Snps")
  imp <- imp[-1,]
  imp[,1] <- as.numeric(imp[,1])
  or <- as.data.frame(t(orig[row,]))
  or$snps <- row.names(or)
  colnames(or) <- c(as.character(or[1,1]), "Snps")
  or <- or[-1,]
  or[,1] <- as.numeric(or[,1])  
  comp <- merge(imp, or, by="Snps")
  corC <- cor(comp[,2], comp[,3], use="pairwise.complete.obs")
  cor <- c(cor, c(sample, corC))
}

write.table(cor, paste0(WorkingDir, "/GenCorr_samples", "NUMBEROFMASKING.txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)


imputed <- imputed[,-1]
orig <- orig[,-1]
corSNP <- data.frame(matrix(ncol= 2, nrow=ncol(imputed)))
colnames(corSNP) <- c("SNP", "cor")
for (col in 1:ncol(imputed)) {
  imputed[,col] <- as.numeric (imputed[,col])
  snp <- colnames(imputed)[col]
  orig[,col] <- as.numeric(orig[,col])
  cor <- cor (imputed[,col], orig[,col], use="complete.obs" )
  corSNP$SNP[col] <- snp
  corSNP$cor[col] <- cor
}

write.table(corSNP, paste0(WorkingDir, "/GenCorr_snps", "NUMBEROFMASKING.txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)
