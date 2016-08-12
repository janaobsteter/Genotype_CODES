#preberi ref za kromosom 1 in 
setwd("~/Genotipi/MS_impute_phased_Ref+Marker_files/")
chr1 <- read.csv("p_All_1kb+MS_BT_ref_chr1.txt", sep=" ", header=F) #header is false due to replicated ids
colnames (chr1)[2] <- "Name"
chr1GP4 <- read.csv ("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGP4v01/Beagle_phased/PedGP4v01chr1.PedGP4v01chr1.bgl.phased", sep=" ", header=F)
colnames(chr1GP4)[2] <- "Name"
chr1GGP <- read.csv ("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGGP/Beagle_phased/PedGGPchr1.PedGGPchr1.bgl.phased", sep=" ", header=F)
colnames(chr1GGP)[2] <- "Name"
chr1GP3 <- read.csv ("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGP3v02/Beagle_phased/PedGP3v02chr1.PedGP3v02chr1.bgl.phased", sep=" ", header=F)
colnames(chr1GP3)[2] <- "Name"
chr1HD <- read.csv ("~/Genotipi/Genotipi1_12042016/Combined_AB/PedHD/Beagle_phased/PedHDchr1.PedHDchr1.bgl.phased", sep=" ", header=F)
colnames(chr1HD)[2] <- "Name"

#dodaj phased GP4 živali za kromosom 1 v referenco
Chr1_ref <- merge (chr1, chr1GP4, all=T, by="Name")
Chr1_ref <- merge (Chr1_ref, chr1GGP, all=T, by="Name")
Chr1_ref <- merge (Chr1_ref, chr1GP3, all=T, by="Name")

Chr1_ref_O <- Chr1_ref [order(Chr1_ref$V1.x),]
Chr1_ref_order <- Chr1_ref [c(which(Chr1_ref$V1.x=="I"),
which(Chr1_ref$V1.x=="FID"),
which(Chr1_ref$V1.x=="MID"),
which(Chr1_ref$V1.x=="S"),
which(Chr1_ref$V1.x=="A"),
c(which(Chr1_ref$V1.x=="M"))), ]
Chr1_ref_order <- Chr1_ref_order[, c(2,1, 3:ncol(Chr1_ref_order))]
                                                                                                                              

