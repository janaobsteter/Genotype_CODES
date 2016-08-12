setwd("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGP4v01/Beagle_phased/")
setwd("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGGP/Beagle_phased/")

setwd("~/Genotipi/Genotipi1_12042016/Combined_AB/PedHD/Beagle_phased/")
setwd("~/Genotipi/Genotipi1_12042016/Combined_AB/Ped50K//Beagle_phased/")
setwd("~/Genotipi/Genotipi1_12042016/Combined_AB/PedGP3v02///Beagle_phased/")


###########################################################################
#poskus branja v loopu, neuspešen
##########################################################################
MS_chr <- read.csv ("~/Genotipi/MS_byChr.csv", sep=",", header=F)
MS_chr_2 <- MS_chr[-2,]
chr <- c(1,3,5,9,15,16,18,19,20,21,23)
Chr_MS1 <- t(as.matrix(colnames(lala)))
colnames (Chr_MS1) <- Chr_MS1[1,]
Chr_MS1[1,] <- NA

for (i in 1:11) {
  nam <- paste("PedGP4v01", MS_chr_2[i,1], ".PedGP4v01", MS_chr_2[i,1], ".bgl.phased", sep="")
  chr_ms <- read.csv (nam, sep=" ", header=T)
  Chr_MS1 <- chr_ms [(chr_ms$id==(MS_chr_2[i,2])),]
  Chr_MS <- rbind(Chr_MS, chr_ms)
}

for (i in chr) {
  nam <- paste( chr, "a.txt", sep="")
  chr_ms <- read.csv (nam, header=T, sep=" ")
  Chr_MS1 <- rbind(Chr_MS1, chr_ms)
}


#####################################################################################################
#deli na roke
######################################################################################################
BM1824 <- read.csv ("1a.txt", header=T, sep=" ")
BM2113 <- read.csv ("2a.txt", header=T, sep=" ")
INRA023 <- read.csv ("3a.txt", header=T, sep=" ")
ETH10 <- read.csv ("5a.txt", header=T, sep=" ")
ETH225 <- read.csv ("9a.txt", header=T, sep=" ")
SPS115 <- read.csv ("15a.txt", header=T, sep=" ")
TGLA53 <- read.csv ("16a.txt", header=T, sep=" ")
TGLA227 <- read.csv ("18a.txt", header=T, sep=" ")
ETH3 <- read.csv ("19a.txt", header=T, sep=" ")
TGLA126 <- read.csv ("20a.txt", header=T, sep=" ")
TGLA122 <- read.csv ("21a.txt", header=T, sep=" ")
BM1818 <- read.csv ("23a.txt", header=T, sep=" ")

MSs <- rbind(BM1824, INRA023)
MSs <-rbind(MSs,BM2113)
MSs <- rbind(MSs, ETH10)
MSs <- rbind(MSs, ETH225)
MSs <- rbind(MSs, SPS115)
MSs <- rbind(MSs, TGLA53)
MSs <- rbind(MSs, TGLA227)
MSs <- rbind(MSs, ETH3)
MSs <- rbind(MSs, TGLA126)
MSs <- rbind(MSs, TGLA122)
MSs <- rbind(MSs, BM1818)


##extract odd and columns
##!!! prilagodi število stolpcev glede na število genotuipiziranih posameznikov na čipu!!!
soda <- seq(from=4, to=68, by=2)
liha <- seq(from=3, to=68, by=2)
MS_alel1 <- MSs[,c(1,2,liha)]
MS_alel2 <- MSs[,c(1,2,soda)]

MS_alel1 <- t(MS_alel1)
colnames(MS_alel1) <- MS_alel1[2,]
MS_alel1 <- MS_alel1[-c(1,2),]

MS_alel2 <- t(MS_alel2)
colnames(MS_alel2) <- MS_alel2[2,]
MS_alel2 <- MS_alel2[-c(1,2),]

write.table(MS_alel1, "Alel1.csv", sep=",") ##pejt v excel pa poprav imena ID in spremeni X v SI
MS_alel1 <- read.table("Alel1.csv", sep=",", header=T)
colnames(MS_alel1)[1] <- "ID"
rownames(MS_alel1) <- MS_alel1[,1]

write.table(MS_alel2, "Alel2.csv", sep=",")
MS_alel2 <- read.table("Alel2.csv", sep=",", header=T)
colnames(MS_alel2) <- c("ID", "BM2113_1","BM1824_1", "INRA023_1", "ETH10_1", "ETH225_1", "SPS115_1", "TGLA53_1", "TGLA227_1", "ETH3_1", "TGLA126_1", "TGLA122_1", "BM1818_1")
rownames(MS_alel2) <- MS_alel2[,1]
#"BM2113"

##merge back together
MS_alel <- cbind (MS_alel1, MS_alel2)
MS_alel <- MS_alel[,c(1,3,16,2,15,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26)]
colnames(MS_alel)[1] <- "ID"
write.table(MS_alel, "MS_imp_alleles.csv", sep=",")
#############################
MS_lab <- read.csv("~/Documents/jana_MS.csv", sep=",", header=T)
Lab_and_imp <- intersect(MS_lab$ID, MS_alel$ID)
Lab_and_imp <- intersect(MS_lab$ID, MS_alel$ID)
Lab_and_imp <- intersect(MS_lab$ID, rownames(MS_alel))

#7 for GP4 chip
MS_lab_7 <- MS_lab[(MS_lab$ID %in% Lab_and_imp),]
MS_alel_7 <- MS_alel[(MS_alel$ID %in% Lab_and_imp),]
MS_lab_7 <- MS_lab_7[order(MS_lab_7$ID),]
MS_alel_7 <- MS_alel_7[order(MS_alel_7$ID),]

#2 for GGP chip 
MS_lab_2 <- MS_lab[(MS_lab$ID %in% Lab_and_imp),]
MS_alel_2 <- MS_alel[(MS_alel$ID %in% Lab_and_imp),]
MS_lab_2 <- MS_lab_2[order(MS_lab_2$ID),]
MS_alel_2<- MS_alel_2[order(MS_alel_2$ID),]

#15 for HD chip
MS_lab_15 <- MS_lab[(MS_lab$ID %in% Lab_and_imp),]
MS_alel_15 <- MS_alel[(MS_alel$ID %in% Lab_and_imp),]
MS_lab_15 <- MS_lab_15[order(MS_lab_15$ID),]
MS_alel_15<- MS_alel_15[order(MS_alel_15$ID),]


#57 for Illumina50K
MS_lab_57 <- MS_lab[(MS_lab$ID %in% Lab_and_imp),]
MS_alel_57 <- MS_alel[(MS_alel$ID %in% Lab_and_imp),]
MS_lab_57 <- MS_lab_57[order(MS_lab_57$ID),]
MS_alel_57<- MS_alel_57[order(MS_alel_57$ID),]

#4 for GeneSeekGp3v02
MS_lab_4 <- MS_lab[(MS_lab$ID %in% Lab_and_imp),]
MS_alel_4 <- MS_alel[(MS_alel$ID %in% Lab_and_imp),]
MS_lab_4 <- MS_lab_4[order(MS_lab_4$ID),]
MS_alel_4<- MS_alel_4[order(MS_alel_4$ID),]

#!!!!!##remove non-imputed MS from lab file
#MS_lab_7 <- MS_lab_7[,-c(4,5)]
#MS_lab_2 <- MS_lab_2[,-c(4,5)]
#MS_lab_15 <- MS_lab_15[,-c(4,5)]
#MS_lab_4 <- MS_lab_4[,-c(4,5)]

write.table(MS_alel_7, "GP4_imp.csv", sep=",")
write.table(MS_lab_7, "GP4_lab.csv",sep=",")
write.table(MS_alel_2, "GGP_imp.csv", sep=",")
write.table(MS_lab_2, "GGP_lab.csv", sep=",")
write.table(MS_alel_15, "HD_imp.csv",sep=",")
write.table(MS_lab_15, "HD_lab.csv",sep=",")
write.table(MS_alel_57, "50K_imp.csv",sep=",")
write.table(MS_lab_57, "50K_lab.csv",sep=",")
write.table(MS_alel_4, "GP3v02_imp.csv",sep=",")
write.table(MS_lab_4, "GP3v02_lab.csv",sep=",")

## change imp in excell to match lab
MS_imp_7 <- read.csv("GP4_imp.csv", sep=",", header=T)
MS_imp_2 <- read.csv("GGP_imp.csv", sep=",", header=T)
MS_imp_15 <- read.csv("HD_imp.csv", sep=",", header=T)
MS_imp_57 <- read.csv("50K_imp.csv", sep=",", header=T)
MS_imp_4 <- read.csv("GP3v02_imp.csv", sep=",", header=T)



###check intersect
Check <- matrix (nrow=7, ncol=12)
rownames(Check) <- MS_imp_7$ID

colnames(Check) <- c("BM1824","BM2113","ETH10","ETH225","ETH3","INRA23","SPS115","TGLA122","TGLA126","TGLA227","TGLA53","BM1818")
#

MS_col <- seq(from=2, to=25, by=2) #imena mikrosatelitov

for (i in MS_col) {
  for (l in 1:7) {
    a <- c(as.character(MS_lab_7[l, i]), as.character(MS_lab_7[l, i+1]))
    b <- c(as.character(MS_imp_7[l, i]), as.character(MS_imp_7[l, i+1]))
    c <- (a %in% b)
    d <- length(which(c=="TRUE"))
    Check[l,(which(MS_col==i))] <- d
  }
}

per <- rowSums(Check)*100/24
Check <- cbind(Check,per)

write.table(Check, "GP4_check_concordance.csv", sep=",")
write.table(Check, "HD_check_concordance.csv", sep=",")
write.table(Check, "GP_check_concordance.csv", sep=",")
write.table(Check, "50K_check_concordance.csv", sep=",")
write.table(Check, "GP3v02_check_concordance.csv", sep=",")

  

