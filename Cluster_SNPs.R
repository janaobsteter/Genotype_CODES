#############################################################################3
#10x CV for smaller chip --> larger chip
###############################################################################
#setwd("/home/jana/Genotipi/Genotipi_CODES")

############################################################################################################3
#1)READ IN MAP file of the CONC CHIP - in CONCPLINK_REF are only SNPs that are also on REF chip
Conc_Map <- read.table("CONCPLINK_REF.map")
Conc_SNPs <- data.frame(SNPs=Conc_Map[,2])


Conc_SNPs$Clust <- NA
#assign SNPs to 10 clusters for cross-validation
CVseq1 <- seq(from=1, to=nrow(Conc_SNPs), by=10) #get row number for each 10th SNP
Conc_SNPs$Clust[CVseq1] <- "C1" #add C1-C10 cluster names to Clust column
CVset1 <- Conc_SNPs[CVseq1,] #create a list of SNPs in order to extract and compare the concordance using PLINK
CVseq2 <- seq(from=2, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq2] <- "C2"
CVset2 <-  Conc_SNPs[CVseq2,] 
CVseq3 <- seq(from=3, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq3] <- "C3"
CVset3 <-  Conc_SNPs[CVseq3,] 
CVseq4 <- seq(from=4, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq4] <- "C4"
CVset4 <-  Conc_SNPs[CVseq4,] 
CVseq5 <- seq(from=5,to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq5] <- "C5"
CVset5 <-  Conc_SNPs[CVseq5,] 
CVseq6 <- seq(from=6,to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq6] <- "C6"
CVset6 <-  Conc_SNPs[CVseq6,] 
CVseq7 <- seq(from=7, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq7] <- "C7"
CVset7 <-  Conc_SNPs[CVseq7,] 
CVseq8 <- seq(from=8, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq8] <- "C8"
CVset8 <-  Conc_SNPs[CVseq8,] 
CVseq9 <- seq(from=9, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq9] <- "C9"
CVset9 <-  Conc_SNPs[CVseq9,] 
CVseq10 <- seq(from=10, to=nrow(Conc_SNPs), by=10)
Conc_SNPs$Clust[CVseq10] <- "C10"
CVset10 <-  Conc_SNPs[CVseq10,] 

#write cluster 10xCV for small-large common SNPs
write.table(Conc_SNPs, "SNP_Cluster.txt", row.names=F, col.names=F, quote=F, sep=" ")

#write SNP set separately

write.table(CVset1, "CVSNPset1.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset2, "CVSNPset2.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset3, "CVSNPset3.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset4, "CVSNPset4.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset5, "CVSNPset5.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset6, "CVSNPset6.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset7, "CVSNPset7.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset8, "CVSNPset8.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset9, "CVSNPset9.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset10, "CVSNPset10.txt", sep="\n", quote=F, row.names=F, col.names=F)

