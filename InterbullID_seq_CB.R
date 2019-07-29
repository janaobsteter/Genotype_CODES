Hol800 <- read.table("~/Genotipi/Genotipi_DATA/Genotipi_latest/Crnobela/IndividualsHolstein_50Kv01.txt")
pedigreCB <- read.csv("~/Genotipi/TransformGeno/Crnobela_seq_ID.csv", header=T)
pedigreCB$STEV_ORIG_ZIVAL <- substr(pedigreCB$ID_ZIVALI, start =3, stop=10)




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
Hol14 <- Hol14[1:32,]


#when you strip at 15
Hol800$V3 <- substring(Hol800$V1, 15)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol15 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]
Hol15 <- Hol15[c(2,3,6),]

#when you strip at 16
Hol800$V3 <- substring(Hol800$V1, 16)
length(intersect(Hol800$V3, pedigreCB$STEV_ORIG_ZIVAL))
Hol16 <- Hol800[which(Hol800$V3 %in% pedigreCB$STEV_ORIG_ZIVAL),]
Hol16 <- Hol16[3,]

#up to 16 - now we have all 192 animals
TillNow <- rbind(Hol8, Hol9, Hol10, Hol11, Hol12, Hol13, Hol14, Hol15, Hol16)
nrow(TillNow)
length(unique(TillNow$V1))
TillNowU <- TillNow[unique(TillNow$V1),]

#merge to get sequences
TillNowU <- TillNowU[,c(1,3)]
colnames(TillNowU) <- c("Interbull", "STEV_ORIG_ZIVAL") 
pedigreCB <- pedigreCB[,c(1,5,6,9,10)]
# CB_seq_num <- merge(TillNowU, pedigreCB, by="STEV_ORIG_ZIVAL")
#CB_seq_num <- CB_seq_num[,c(2,3)]
write.table(CB_seq_num, "~/Genotipi/Genotipi_DATA/Genotipi_latest/Crnobela//InterbullID_seq_CB.csv", row.names=F, quote=F, sep=",")


#pridobi podatke o Holstein genotipiziranih
DataCB <- pedigreCB[(which(pedigreCB$ZIV_ID_SEQ %in% CB_seq_num$ZIV_ID_SEQ)),]
DataCB$DAT_ROJSTVO <- as.Date(DataCB$DAT_ROJSTVO, format="%d.%m.%Y")
DataCB$leto <- format(DataCB$DAT_ROJSTVO, format="%Y")


#pripravi še za vnos v govedo
colnames(TillNow) <- c("Interbull", "CHIP", "ZIVID")
GovedoInd <- merge(TillNow, CB_seq_num, by="Interbull")
GovedoInd <- GovedoInd[,c(5,2)]
date <- as.Date("20160728", format="%Y%m%d")
colnames(GovedoInd) <- c ("ZIV_ID_SEQ", "GEN_CHIP")
GovedoInd$V3 <- date
colnames(GovedoInd)[3] <- c ( "GEN_DATUM")
write.table(GovedoInd, "~/Genotipi/Genotipi_DATA/Holstein/GovedoInd.csv", quote=F, row.names=F, sep=",")
