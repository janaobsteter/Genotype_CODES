Chr2markers <- read.csv ("~/Genotipi/MS_impute_phased_Ref+Marker_files/MinSNP+MS_chr2.txt", sep=" ", header=F)
colnames(Chr2markers) <- c("Marker", "Position", "Alele1", "Alele2")

#combined Imputed
#read in bgl file - header is false due to duplicated ID (colnames)
Chr2 <- read.csv("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/NoHD/OUTPUT/Beagle_imputedMS/Chr2.bgl", header=F, sep=" ")
Chr2 <- read.csv("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT/Chr2.bgl", header=F, sep=" ")
colnames (Chr2) [2] <- ("Marker")
Chr2 <- merge(Chr2markers, Chr2, by="Marker") #merge McClure marker file for chr2 and ped-->bgl gneotype file
Chr2 <- Chr2 [, 1:(ncol(Chr2)-2)] # remove empty columns
Chr2 <- Chr2 [order(Chr2$Position),] # order file according to McClure chr2 marker file
Chr2 <- Chr2 [, -c(4,3,2)] # remove merged columns
Chr2 <- Chr2 [,c(2,1,3:ncol(Chr2))] #reorder
write.table (Chr2, "/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT/Beagle_100/Chr2_sorted.bgl", sep=" ", row.names=F, col.names=F, quote=F)
#write the table, then go into excel and add first two rows with id and A
#read the table back into and save it in proper format
CHR <- read.csv ("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT//Chr2_sorted.bgl", sep=",", header=F)
write.table(CHR, "/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT//Chr2_sorted.bgl", sep="\t", row.names=F, col.names=F, quote=F)

#GP4 chip
GP4chr2 <- read.csv("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGP4v01/PedGP4v01chr2.bgl", header=F, sep=" ")
colnames (GP4chr2) [2] <- ("Marker")
Chr2 <- merge(Chr2markers, GP4chr2, by="Marker")
Chr2 <- Chr2 [, 1:71]
Chr2 <- Chr2 [order(Chr2$Position),]
Chr2 <- Chr2 [, -c(4,3,2)]
Chr2 <- Chr2 [,c(2,1,3:68)]

write.table (Chr2, "GP4Chr2_sorted.bgl", sep=" ")
CHR <- read.csv ("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGP4v01/GP4Chr2_sorted.csv", sep=",", header=F)
write.table(CHR, "Chr2.csv", sep="\t", row.names=F, col.names=F)


#GGP chip
setwd("/home/janao/Genotipi/Genotipi1_12042016/Combined_AB/PedGGP")
GGPchr2 <- read.csv("/home/janao/Genotipi/Genotipi1_12042016/Combined_AB/PedGGP/PedGGPchr2.bgl", header=F, sep=" ")
colnames (GGPchr2) [2] <- ("Marker")
Chr2 <- merge(Chr2markers, GGPchr2, by="Marker")
Chr2 <- Chr2 [, 1:17]
Chr2 <- Chr2 [order(Chr2$Position),]
Chr2 <- Chr2 [, -c(4,3,2)]
Chr2 <- Chr2 [,c(2,1,3:14)]

write.table (Chr2, "GGPChr2_sorted.bgl", sep=" ")

#GP3v02 chip
setwd("/home/janao/Genotipi/Genotipi1_12042016/Combined_AB/PedGP3v02")
GP3chr2 <- read.csv("PedGP3v02chr2.bgl", header=F, sep=" ")
colnames (GP3chr2) [2] <- ("Marker")
Chr2 <- merge(Chr2markers, GP3chr2, by="Marker")
Chr2 <- Chr2 [, 1:35]
Chr2 <- Chr2 [order(Chr2$Position),]
Chr2 <- Chr2 [, -c(4,3,2)]
Chr2 <- Chr2 [,c(2,1,3:32)]

write.table (Chr2, "GP3Chr2_sorted.bgl", sep=" ")

#HD chip
setwd("/home/janao/Genotipi/Genotipi1_12042016/Combined_AB/PedHD")
HDchr2 <- read.csv("PedHDchr2.bgl", header=F, sep=" ")
colnames (HDchr2) [2] <- ("Marker")
Chr2 <- merge(Chr2markers, HDchr2, by="Marker")
Chr2 <- Chr2 [, 1:67]
Chr2 <- Chr2 [order(Chr2$Position),]
Chr2 <- Chr2 [, -c(4,3,2)]
Chr2 <- Chr2 [,c(2,1,3:64)]

write.table (Chr2, "HDChr2_sorted.bgl", sep=" ")
