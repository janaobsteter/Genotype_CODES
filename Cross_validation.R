###########################################################################3
#create a script to cross-validate imputation using FIMPUTE
#############################################################################
setwd("~/Genotipi/Genotipi1_12042016/PLINK_genotypeFiles/")
#1) READ IN MAP FILES AND FIND THE SNPS COMMON TO ALL CHIPS
GGP_map <- read.csv("GGP/OUTPUT/PLINK_MERGED.map", sep=" ", header=F)
GGP3_map <- read.csv("GGPv03/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
GP3_map <- read.csv("GP3v02/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
GP4_map <- read.csv("GP4/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
HD_map <- read.csv("HD/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
K50_map <- read.csv("50K/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
HD138K_map <- read.csv("/home/janao/Genotipi/Genotipi03062016/PLINK_genotypeFiles/HD138K/OUTPUT/PLINK_MERGED.map", sep=" ", header=F)
min <- read.csv("~/Genotipi/MinSNP_Set.txt", header=F)

common <- intersect(min$V1, HD138K_map$V2)

common <- intersect (GGP_map$V2, GGP3_map$V2)
common <- intersect (GGP_map$V2, GP3_map$V2)
common <- intersect (common, GP3_map$V2)
common <- intersect (common, GP4_map$V2)
common <- intersect (common, HD_map$V2)
common <- intersect (common, HD138K_map$V2)
common <- intersect (HD_map$V2, GP4_map$V2)
common <- intersect (GGP_map$V2, GP3_map$V2)
common <- intersect (GP4_map$V2, GP3_map$V2)
common2 <- intersect (GP4_map$V2, GGP_map$V2)
common <- intersect (K50_map$V2, GP4_map$V2)
common <- intersect (K50_map$V2, GP3_map$V2)
common1 <- intersect(GP4_map$V2, HD138K_map$V2)
common <- which(HD138K_map$V2 %in% GP4_map$V2)

#####################
#read in excluded maps
#######################

setwd("~/Genotipi/Genotipi1_12042016/StepImputation/")
#1) READ IN MAP FILES AND FIND THE SNPS COMMON TO ALL CHIPS
GGP_map <- read.csv("Step1/GGP_exc.map", sep="\t", header=F)
GP3_map <- read.csv("Step2/GP3_exc.map", sep="\t", header=F)
GP4_map <- read.csv("Step2/GP41_29.map", sep="\t", header=F)
HD_map <- read.csv("Step3/HD_exc.map", sep="\t", header=F)
HD138K_map <- read.csv("StepAllInOne/HD138K_exc.map", sep="\t", header=F)
min <- read.csv("~/Genotipi/MinSNP_Set.txt", header=F)

common1 <- 

"#create a cluster file
commonClust <- data.frame (common)
commonClust$Clust <- NA

#assign SNPs to 5 clusters for cross-validation
CV5seq1 <- seq(from=1, to=length(common), by=5)
commonClust$Clust[CV5seq1] <- "C1"
CVset1 <- common[CV5seq1] 
CV5seq2 <- seq(from=2, to=length(common), by=5)
commonClust$Clust[CV5seq2] <- "C2"
CVset2 <- common[CV5seq2] 
CV5seq3 <- seq(from=3, to=length(common), by=5)
commonClust$Clust[CV5seq3] <- "C3"
CVset3 <- common[CV5seq3] 
CV5seq4 <- seq(from=4, to=length(common), by=5)
commonClust$Clust[CV5seq4] <- "C4"
CVset4 <- common[CV5seq4] 
CV5seq5 <- seq(from=5, to=length(common), by=5)
commonClust$Clust[CV5seq5] <- "C5"
CVset5 <- common[CV5seq5] 

write.table(commonClust, "/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/ClusterSNPs_CV5_AllChips.txt", row.names=F, col.names=F, quote=F, sep=" ")
write.table(commonClust, "/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/ClusterSNPs_CV5_GSChips.txt", row.names=F, col.names=F, quote=F, sep=" ")
write.table(commonClust, "/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/ClusterSNPs_CV5_GS30KChips.txt", row.names=F, col.names=F, quote=F, sep=" ")

#create within file for PLINK, individuals and clusters
IndClust <- read.csv ("/home/janao/Genotipi/Genotipi1_12042016/GS30KChips/OUTPUT/FIDandID.txt", sep=" ", header=F)
setwd("/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/GS30KChips")
IndClust$V3 <- NA
IndClust$V3 <- "C1" # cluster for C1
IndClust$V3[which(IndClust$V2=="SI14515440")] <- "C2"
write.table (IndClust, "IndClus1.txt", quote=F, row.names=F, col.names=F, sep=" ")
IndClust$V3 <- NA
IndClust$V3 <- "C2" # cluster for C2
IndClust$V3[which(IndClust$V2=="SI14515440")] <- "C1"
write.table (IndClust, "IndClus2.txt", quote=F, row.names=F, col.names=F, sep=" ")
IndClust$V3 <- NA
IndClust$V3 <- "C3" # cluster for C3
IndClust$V3[which(IndClust$V2=="SI14515440")] <- "C1"
write.table (IndClust, "IndClus3.txt", quote=F, row.names=F, col.names=F, sep=" ")
IndClust$V3 <- NA
IndClust$V3 <- "C4" # cluster for C4
IndClust$V3[which(IndClust$V2=="SI14515440")] <- "C1"
write.table (IndClust, "IndClus4.txt", quote=F, row.names=F, col.names=F, sep=" ")
IndClust$V3 <- NA
IndClust$V3 <- "C5" # cluster for C5
IndClust$V3[which(IndClust$V2=="SI14515440")] <- "C1"
write.table (IndClust, "IndClus5.txt", quote=F, row.names=F, col.names=F, sep=" ")



#write SNP set separately
write.table(CVset1, "CVSNPset1.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset2, "CVSNPset2.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset3, "CVSNPset3.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset4, "CVSNPset4.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset5, "CVSNPset5.txt", sep="\n", quote=F, row.names=F, col.names=F)

#################################################################
#find distinctive SNPs, combined
####################################################################
GGP <- as.vector(GGP_map$V2)
GGP3 <- as.vector(GGP3_map$V2)
GP4 <- as.vector(GP4_map$V2)
GP3 <- as.vector(GP3_map$V2)
HD <- as.vector(HD_map$V2)
K50 <- as.vector(K50_map$V2)

combine <- c (GGP, GGP3)
combine <- c (combine, GP3)
combine <- c (combine, GP4)
combine <- c (combine, HD)
combine <- c (combine, K50)
length(unique(combine))

#############################################################################3
#10x CV for 19K --> 26K
###############################################################################

#which SNPs in GGP also in GP3 - get row numbers
common <- which (GGP_map$V2 %in% GGP3_map$V2) #19,279 SNPs in common, 445 / 441 excluded
#get a dataframe from GGP_map with the common SNPs
commonGGP <- GGP_map[c(common),]
#generate a sample of 10000 randomly chosen numbers for rows from 19279 GGPcommon for cross validation
#K10Sample <- sample(1:length(common), 10000, replace=FALSE)
#extract those 10K SNPs
#common10K <- commonGGP[c(K10Sample),]

#exclude sex SNPs and 0 chr SNPs
SNP0 <- which(commonGGP$V1==0)
SNPX <- which(commonGGP$V1=="X")
SNPY <- which(commonGGP$V1=="Y")
SNPMT <- which(commonGGP$V1=="MT")
ExcludeSNPs <- c(SNP0, SNPX, SNPY, SNPMT)
commonGGP <- commonGGP[-c(ExcludeSNPs),]


#create a cluster file with the list of 10,000 SNPs and 10 clusters
commonClust <- data.frame(commonGGP$V2)
commonClust$Clust <- NA

#assign SNPs to 5 clusters for cross-validation
CVseq1 <- seq(from=1, to=nrow(commonClust), by=10) #get row number for each 10th SNP
commonClust$Clust[CVseq1] <- "C1" #add C1-C10 cluster names to Clust column
CVset1 <- commonClust$commonGGP.V2[CVseq1] #create a list of SNPs in order to extract and compare the concordance using PLINK
CVseq2 <- seq(from=2, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq2] <- "C2"
CVset2 <-  commonClust$commonGGP.V2[CVseq2] 
CVseq3 <- seq(from=3, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq3] <- "C3"
CVset3 <-  commonClust$commonGGP.V2[CVseq3] 
CVseq4 <- seq(from=4, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq4] <- "C4"
CVset4 <-  commonClust$commonGGP.V2[CVseq4] 
CVseq5 <- seq(from=5,to=nrow(commonClust), by=10)
commonClust$Clust[CVseq5] <- "C5"
CVset5 <-  commonClust$commonGGP.V2[CVseq5] 
CVseq6 <- seq(from=6,to=nrow(commonClust), by=10)
commonClust$Clust[CVseq6] <- "C6"
CVset6 <-  commonClust$commonGGP.V2[CVseq6] 
CVseq7 <- seq(from=7, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq7] <- "C7"
CVset7 <-  commonClust$commonGGP.V2[CVseq7] 
CVseq8 <- seq(from=8, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq8] <- "C8"
CVset8 <-  commonClust$commonGGP.V2[CVseq8] 
CVseq9 <- seq(from=9, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq9] <- "C9"
CVset9 <-  commonClust$commonGGP.V2[CVseq9] 
CVseq10 <- seq(from=10, to=nrow(commonClust), by=10)
commonClust$Clust[CVseq10] <- "C10"
CVset10 <-  commonClust$commonGGP.V2[CVseq10] 

#write cluster 10xCV for 19K-26K common SNPs
write.table(commonClust, "/home/janao/Genotipi/Genotipi1_12042016/StepImputation/19K26K_Cluster.txt", row.names=F, col.names=F, quote=F, sep=" ")

#write SNP set separately
setwd("/home/janao/Genotipi/Genotipi1_12042016/StepImputation")
write.table(CVset1, "19K26K_CVSNPset1.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset2, "19K26K_CVSNPset2.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset3, "19K26K_CVSNPset3.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset4, "19K26K_CVSNPset4.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset5, "19K26K_CVSNPset5.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset6, "19K26K_CVSNPset6.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset7, "19K26K_CVSNPset7.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset8, "19K26K_CVSNPset8.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset9, "19K26K_CVSNPset9.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(CVset10, "19K26K_CVSNPset10.txt", sep="\n", quote=F, row.names=F, col.names=F)

