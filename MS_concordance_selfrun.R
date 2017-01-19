#####################################################################################################
#imputation accuracy for imputed data with GP4 reference
######################################################################################################
#setwd ("~/Genotipi/Genotipi1_12042016/Imputation_GP4ref/OUTPUT/Beagle_imputedMS/")
#setwd ("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT/Beagle_imputedMS")
#setwd("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GGPv03ref/OUTPUT/Beagle_imputedMS")
#setwd("/home/janao/Genotipi/Genotipi1_12042016/Imputation_HDref/OUTPUT/Beagle_imputedMS")
#setwd("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/NoHD/OUTPUT/Beagle_imputedMS")
#setwd("/home/janao/Genotipi/MS_Imputation/MS_03062016/BeagleImputed10it")
#setwd("/home/janao/Genotipi/MS_Imputation/MS_03062016/BeagleImputed100it")
#setwd("/home/janao/Genotipi/MS_Imputation/MS_03062016/BeagleImputed1000it")
#setwd("~/Genotipi/MS_Imputation/MS_03062016/AddedRefImputation_10it")
setwd("~/Genotipi/MS_Imputation/MS_05012017/10it")


#read in MS file, each one with additional ID header
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

#merge by IDs
#, BM2113, BM1824, ETH10, ETH225, ETH3, INRA23, SPS115, TGLA122, TGLA126, TGLA227, TGLA53, BM1818
MSs <- rbind(BM2113, BM1824)
MSs <-rbind(MSs,ETH10)
MSs <- rbind(MSs, ETH225)
MSs <- rbind(MSs, ETH3)
MSs <- rbind(MSs, INRA023)
MSs <- rbind(MSs, SPS115)
MSs <- rbind(MSs, TGLA122)
MSs <- rbind(MSs, TGLA126)
MSs <- rbind(MSs, TGLA227)
MSs <- rbind(MSs, TGLA53)
MSs <- rbind(MSs, BM1818)


##extract odd and columns
##!!! prilagodi število stolpcev glede na število genotuipiziranih posameznikov na čipu!!!
soda <- seq(from=4, to=ncol(MSs), by=2)
liha <- seq(from=3, to=ncol(MSs), by=2)
MS_alel1 <- MSs[,c(1,2,liha)]
MS_alel2 <- MSs[,c(1,2,soda)]

MS_alel1 <- t(MS_alel1)
colnames(MS_alel1) <- MS_alel1[2,]
MS_alel1 <- MS_alel1[-c(1,2),]

MS_alel2 <- t(MS_alel2)
colnames(MS_alel2) <- MS_alel2[2,]
MS_alel2 <- MS_alel2[-c(1,2),]
rownames(MS_alel2) <- rownames(MS_alel1)

#write.table(MS_alel1, "Alel1.csv", sep=",") ##pejt v excel pa poprav imena ID in spremeni X v SI
#MS_alel1 <- read.table("Alel1.csv", sep=",", header=T)
#colnames(MS_alel1)[1] <- "ID"
#rownames(MS_alel1) <- MS_alel1[,1]

#write.table(MS_alel2, "Alel2.csv", sep=",")
#MS_alel2 <- read.table("Alel2.csv", sep=",", header=T)
colnames(MS_alel2) <- c("BM2113_1","BM1824_1", "ETH10_1", "ETH225_1", "ETH3_1", "INRA023_1", "SPS115_1", "TGLA122_1", "TGLA126_1", "TGLA227_1","TGLA53_1", "BM1818_1")
#rownames(MS_alel2) <- MS_alel2[,1]
#"BM2113"

##merge back together
MS_alel <- cbind (MS_alel1, MS_alel2)
MS_alel <- as.table(MS_alel[,c(1,13,2,14,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22, 11, 23, 12, 24)])
#colnames(MS_alel)[1] <- "ID"
#write.table(MS_alel, "MS_imp_alleles.csv", sep=",")
write.table(MS_alel, "MS_imputed.csv", sep=",", quote=F)
MS_alel <- read.table ("MS_imputed.csv", sep=",", header=T)


MS_alel$ID <- rownames(MS_alel)
#remove1 <- which(MS_alel$ID=="SI24273440")
#remove2 <- which(MS_alel$ID=="SI14514119")
#MS_alel <- MS_alel[-c(remove1, remove2),]


#####################################################
#read in lab MSs
#######################################################
#MS_alel <- read.table ("MS_imputed.csv", sep=",", header=T)
MS_lab <- read.csv("/home/janao/Genotipi/MS_Imputation/MS_05012017/LabMS05012017.csv", sep=",", header=T)
Lab_and_imp <- intersect(MS_lab$ID, MS_alel$ID)


#extract 68 animals with lab MS, 67for GeneSeek imputation, 51 for no HD
MS_lab_68<- MS_lab[(MS_lab$ID %in% Lab_and_imp),]
MS_alel_68<- MS_alel[(MS_alel$ID %in% Lab_and_imp),]
#write.table(MS_lab_68, "MS_lab_subset.csv", sep=",")
MS_lab_68 <- MS_lab_68[order(MS_lab_68$ID),]
MS_alel_68 <- MS_alel_68[order(MS_alel_68$ID),]
#write.table(MS_alel_68, "MS_imp_subset.csv", sep=",")
#remane to get the same column names and then reorder to get the same column order
colnames(MS_alel_68)[c(11,12)] <- c("INRA23", "INRA23_1")
MS_alel_68 <- MS_alel_68[, c(3,4,1,2,5:25)]
#add rownames SInumbers and remove the first column with IDs
rownames(MS_lab_68) <- MS_lab_68[,1]
MS_lab_68 <- MS_lab_68[,-1]




###check intersect
#remove last ID column from imputed table
MS_alel_68 <- MS_alel_68 [,1:24]
#substitute empty spaces in TGLA227 columns
MS_alel_68$TGLA227 <- gsub(" ", "", MS_alel_68$TGLA227)
MS_alel_68$TGLA227_1 <- gsub(" ", "", MS_alel_68$TGLA227_1)

#set matrix and variables for the concordance check
MS_col <- seq(from=1, to=23, by=2) #imena mikrosatelitov
Check <- as.data.frame(matrix (nrow=nrow(MS_lab_68), ncol=12))
rownames(Check) <- rownames(MS_alel_68)#imena živali - oba fajla morata biti order by ID!!!
colnames(Check) <- c("BM1824","BM2113","ETH10","ETH225","ETH3","INRA23","SPS115","TGLA122","TGLA126","TGLA227","TGLA53","BM1818")
for (i in MS_col) {
  for (l in 1:nrow(MS_lab_68)) {
    a <- c(as.character(MS_lab_68[l, i]), as.character(MS_lab_68[l, i+1]))
    b <- c(as.character(MS_alel_68[l, i]), as.character(MS_alel_68[l, i+1]))
    if (a[1]==a[2]) {
      if (b[1]==b[2]) {
        if (a[1]==b[1]) {
          Check[l,(which(MS_col==i))] <- 2  
        }
        else {
          Check[l,(which(MS_col==i))] <- 0 
        }
      }
      else if (b[1]!=b[2]) {
        c <- length(intersect(a,b))
        Check[l,(which(MS_col==i))] <- c
      }
    }
    else if (a[1]!=a[2]) {
      c <- (a %in% b)
      d <- sum(c=="TRUE")
      Check[l,(which(MS_col==i))] <- d }
  }
}



per <- rowSums(Check)/24
Check <- cbind(Check,per)
ConcByMS <- colSums(Check) / (nrow(Check)*2)


write.table(Check, "Concordance_Imputed_NewrefBase.csv", sep=",")


######################################
#tabela za Andrejo
########################################
#MSalel od prej - 25 stolpcev - ali preberi "MS_imputed.csv"
#MS_alel <- read.csv("MS_imputed.csv", sep=",")
MS_Andreja <- data.frame(nrow=1, ncol=1)
colnames(MS_Andreja) <- ("MS_len")

#loop - postavi stolpce eden pod drugega --> samo en stolpec, ponovljeno ime živali
MS_col_all <- matrix(ncol=3, nrow=1)
colnames(MS_col_all) <- c("ID", "MS_len", "MS_code")
for (i in 1:(ncol(MS_alel)-1)) {
  MS_col <- data.frame(MS_alel[,c(25,i)])
  colnames(MS_col) <- c("ID", "MS_len")
  MS_col$MS_code <- colnames(MS_alel[i])
  MS_col_all <- rbind(MS_col_all, MS_col)
}

MS_col_all <- MS_col_all[-1,] # izbirši prvo vrstico - nastavljena za začetek loopa
#order po ID živali
Ani_MS <- seq(from=1, to=(nrow(MS_alel)*24), by=24) #--> nastavi variable vsaka 24 vrstica za začetek pisanja MS kode (24 alelov)



MS_list <- read.csv("~/Genotipi/MS_Imputation/MS_names.csv", sep=",", header=F) #imena MS
head(MS_list)
MS_col_all$MS_code <- gsub("_1", "", MS_col_all$MS_code) #zamenjej _1 pri alelih
MS_sorted <- MS_col_all[order(MS_col_all$ID, MS_col_all$MS_code, MS_col_all$MS_len),]

#write.table(MS_col_all, "MS_col_all.csv", sep=" ", row.names=F) #--> pojdi v excell in sortiraj po ID, kodi in dolžini!!!
#MS_sorted <- read.table("MS_col_all.csv", sep=" ", header=T)
MS_sorted$Code <- NA #stolpec za MS kodo


Codes <- c(27,28,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26) #kode za MS kot sortirane po abecedi
MS_sorted$Code <- NA
Ani_MS
length(Codes)

#loop - pojdi v vsako 24 vrstico in dodaj kode po vrsti v naslednjih 24 vrstic - MS morajo biti sortirani!!!!
for (i in Ani_MS) {
  for (l in 0:23) {
    MS_sorted$Code[i+l] <- Codes[l+1]
}
}
head(MS_sorted)

MS_sorted <- MS_sorted[,c(1,4,2)]

#dodaj sekvence živali
seq <- read.csv("~/Genotipi/Genotipi_CODES/Rjave_seq_ID.csv", sep=",")
colnames(seq) [1] <- "ID"
nrow(MS_sorted)
MS_sorted <- merge (MS_sorted, seq, by="ID")
head(MS_sorted)
nrow(MS_sorted)
MS_sorted <- MS_sorted[,c(4,2,3)]
colnames(MS_sorted) <- c("ZGP_ZIV_ID_SEQ","ZGP_SGPL_SIFRA","ZGP_VREDNOST")
head(MS_sorted)
write.table(MS_sorted, "Imputed_MS_Andreja.csv", sep=" ", row.names=F, quote=F)


IDs <- as.data.frame(rownames(MS_alel_68))
colnames(IDs) <- "ID"
ID <- merge(IDs, seq, by="ID")
