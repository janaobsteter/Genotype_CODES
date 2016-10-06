Hol800 <- read.table("~/Genotipi/Genotipi_DATA/Holstein/Genotipi_09082016/GenotypedInd")
pedigreCB <- read.csv("~/Genotipi/Genotipi_DATA/Holstein/Pedigre_All.csv", header=T)




#when you strip at 8
Hol800$V3 <- substring(Hol800$V1, 8)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol8 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#when you strip at 9
Hol800$V3 <- substring(Hol800$V1, 9)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol9 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#when you strip at 10
Hol800$V3 <- substring(Hol800$V1, 10)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol10 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#when you strip at 11
Hol800$V3 <- substring(Hol800$V1, 11)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol11 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#when you strip at 12
Hol800$V3 <- substring(Hol800$V1, 12)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol12 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]


#when you strip at 13
Hol800$V3 <- substring(Hol800$V1, 13)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol13 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#when you strip at 14
Hol800$V3 <- substring(Hol800$V1, 14)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol14 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]


#when you strip at 15
Hol800$V3 <- substring(Hol800$V1, 15)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol15 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#when you strip at 16
Hol800$V3 <- substring(Hol800$V1, 16)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol16 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]

#up to 16 - now we have all 192 animals
TillNow <- rbind(Hol8, Hol9, Hol10, Hol11, Hol12, Hol13, Hol14, Hol15, Hol16)
nrow(TillNow)
length(unique(TillNow$V1))
TillNowU <- TillNow[unique(TillNow$V1),]

#merge to get sequences
TillNowU <- TillNowU[,c(1,3)]
colnames(TillNowU) <- c("Interbull", "STEV_ORIG_ZIVAL") 
pedigreCB <- pedigreCB[,c(1,5,6)]
CB_seq_num <- merge(TillNowU, pedigreCB, by="STEV_ORIG_ZIVAL")
write.table(CB_seq_num, "~/Genotipi/Genotipi_DATA/Holstein/InterbullID_seq.csv", row.names=F, quote=F)




##########################################################################################
##########################################################################################
#Lisaste
##############################################################################################
##############################################################################################
Ids <- read.csv("/home/janao/Genotipi/Genotipi_DATA/Lisaste/LsID_IB.csv", header=F)
pedigreLS <- read.csv("/home/janao/Genotipi/Genotipi_DATA/Lisaste/Pedigree.csv")
pedigreLS$ID_ZIVALI <- as.character(pedigreLS$ID_ZIVALI)
 

#when you strip at 12
Ids$V7 <- substring(Ids$V1, 12) #--> 97/99 animals
length(intersect(Ids$V7, pedigreLS$ID_ZIVALI))
Ls12A <- Ids[which(Ids$V7 %in% pedigreLS$ID_ZIVALI),]
Ls12B <- pedigreLS[which(pedigreLS$ID_ZIVALI %in% Ids$V7 ),]
Ls12 <- merge(Ls12A, Ls12B, by="ID_ZIVALI")
colnames(Ls12A) <- c("IB_ID", "ID_ZIVALI")
Ls12 <- Ls12[,c(2,4)]
colnames(Ls12) <- c("ID_ZIVALI", "SEQZIVALI")

write.table(Ls12, "~/Genotipi/Genotipi_DATA/Lisaste/InterbullID_seq.csv", row.names=F, quote=F, sep=",")
