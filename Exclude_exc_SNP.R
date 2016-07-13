#############################################################################3
#10x CV for smaller chip --> larger chip
###############################################################################
setwd("~/Genotipi/Genotipi1_12042016/PLINK_genotypeFiles/")

############################################################################################################3
#1)READ IN MAP FILES AND FIND THE SNPS COMMON TO ALL CHIPS
GGP_map <- read.csv("GGP/OUTPUT/PLINK_MERGED.map", sep=" ", header=F)
GGP3_map <- read.csv("GGPv03/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
GP3_map <- read.csv("GP3v02/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
GP4_map <- read.csv("GP4/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
HD_map <- read.csv("HD/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
K50_map <- read.csv("50K/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)

################maps from StepImputation
GP3_Imp_map <- read.csv("/home/janao/Genotipi/Genotipi1_12042016/StepImputation/Step1/26K.map", sep="\t", header=F)
##############################################################################################################

######################################
####################################
#CHOSE SNPS FOR IMPUTATION STEP
################################
Small_map <- GP3_Imp_map
Large_map <- GP4_map
#################################
#################################

#setwd to current working directory
setwd("/home/janao/Genotipi/Genotipi1_12042016/StepImputation/Step2")

##################################################################################
#get a list of SNPs exclusively on small chip
#exclude them prior to masking and imputation using PLINK
###################################################################################
SmallExcSNP <- which (!(Small_map$V2 %in% Large_map$V2)) #SNPs only on the smaller chip
LargeExcSNP <- which (!(Large_map$V2 %in% Small_map$V2)) # SNPs only on the larger chip
SmallComSNPs <- data.frame(Small_map$V2[c(SmallExcSNP)]) #keep only the common SNPs on the smaller chip
length(which(!(SmallComSNPs$V2 %in% Large_map$V2))) #check if there is any exclusive SNPs left
write.table(SmallComSNPs, "Exclude_from_Small.txt", row.names=F, col.names=F, quote=F)

############################################################################################3
#get common SNPs for 10x cross validation
######################################################################################

#which SNPs in common on smaller and larger chip - get row numbers
common <- which (Small_map$V2 %in% Large_map$V2) 
#get a dataframe from Small_map with the common SNPs
commonSmall <- Small_map[c(common),]


#exclude sex SNPs and 0 chr SNPs --> No need for this in subsequent steps since you read in the file from imputation, these SNPs already removed
#SNP0 <- which(commonSmall$V1==0)
#SNPX <- which(commonSmall$V1=="X")
#SNPY <- which(commonSmall$V1=="Y")
#SNPMT <- which(commonSmall$V1=="MT")
#ExcludeSNPs <- c(SNP0, SNPX, SNPY, SNPMT)
#commonSmall <- commonSmall[-c(ExcludeSNPs),]


#create a cluster file with the list of 10,000 SNPs and 10 clusters
commonClust <- data.frame(commonSmall$V2)
commonClust$Clust <- NA

#assign SNPs to 5 clusters for cross-validation
CVseq1 <- seq(from=1, to=nrow(commonClust), by=10) #get row number for each 10th SNP
commonClust$Clust[CVseq1] <- "C1" #add C1-C10 cluster names to Clust column
CVset1 <- commonClust$commonSmall.V2[CVseq1] #create a list of SNPs in order to extract and compare the concordance using PLINK
CVseq2 <- seq(from=2, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq2] <- "C2"
CVset2 <-  commonClust$commonSmall.V2[CVseq2] 
CVseq3 <- seq(from=3, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq3] <- "C3"
CVset3 <-  commonClust$commonSmall.V2[CVseq3] 
CVseq4 <- seq(from=4, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq4] <- "C4"
CVset4 <-  commonClust$commonSmall.V2[CVseq4] 
CVseq5 <- seq(from=5,to=nrow(commonClust), by=10)
commonClust$Clust[CVseq5] <- "C5"
CVset5 <-  commonClust$commonSmall.V2[CVseq5] 
CVseq6 <- seq(from=6,to=nrow(commonClust), by=10)
commonClust$Clust[CVseq6] <- "C6"
CVset6 <-  commonClust$commonSmall.V2[CVseq6] 
CVseq7 <- seq(from=7, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq7] <- "C7"
CVset7 <-  commonClust$commonSmall.V2[CVseq7] 
CVseq8 <- seq(from=8, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq8] <- "C8"
CVset8 <-  commonClust$commonSmall.V2[CVseq8] 
CVseq9 <- seq(from=9, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq9] <- "C9"
CVset9 <-  commonClust$commonSmall.V2[CVseq9] 
CVseq10 <- seq(from=10, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq10] <- "C10"
CVset10 <-  commonClust$commonSmall.V2[CVseq10] 

#write cluster 10xCV for small-large common SNPs
write.table(commonClust, "SNP_Cluster.txt", row.names=F, col.names=F, quote=F, sep=" ")

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

