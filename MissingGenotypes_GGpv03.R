inds1 <- read.csv("/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Genotipi/Top_PLINKIDBv03/inds.txt", header=FALSE)
inds2 <- read.csv("/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Genotipi/Top_PLINK/IDBv03/inds.txt", header=FALSE)
inds3 <- read.csv("/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/IDBv03/inds.txt", header=FALSE)
inds4 <- read.csv("~/inds.txt") #To je iz goveda
indsRjava <- read.csv("~/indsRjaveChip.txt") #To je iz goveda
indsR50K <- read.csv("~/inds50K.txt") #To je iz goveda


npp <- read.csv("~/Genotipi/Genotipi_DATA/VSE_GENOTIPIZACIJE_GNPP.csv")
sum(inds3$V1 %in% npp$SAMPLE.ID)
npp[npp$SAMPLE.ID %in% inds2$V1,]
nrow(npp)
table(npp$X)

sum(inds2$V1 %in% inds3$V1)
sum(inds1$V1 %in% inds3$V1)
inds4$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL <- as.character(inds4$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL)
inds3$V1 <- as.character(inds3$V1)
sum(inds4$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL %in% inds3$V1)
inds3[!(inds3$V1 %in% inds4$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL),]
length(inds3[!(inds3$V1 %in% inds4$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL),])
nrow(inds4)

top <- read.table("/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Genotipi/Top_PLINK/ALLIND.txt")
ab <- read.table("/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Genotipi/Genotipi_latest/Rjava/ALLIND.txt")
sum(top$V1 %in% ab$V1)
ab$V1[!(ab$V1 %in% top$V1)]
length(ab$V1[!(ab$V1 %in% top$V1)])
top$V1[!(top$V1 %in% ab$V1)]
miss <- as.character(top$V1[!(top$V1 %in% ab$V1)])
tabelaGovedo[tabelaGovedo$IID %in% miss,]
length(top$V1[!(top$V1 %in% ab$V1)])

sum(top$V1 %in% tabelaGovedo$IID)
sum(ab$V1 %in% tabelaGovedo$IID)

indsDate <- read.csv("/home/jana/Genotipi/Genotipi_DATA/Rjava_TEMP/ALL_INDDATE.txt")
sum(miss$IID %in% indsDate$X0)
miss[(miss$IID %in% indsDate$X0),]


setwd("/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Genotipi/Top_PLINK/")
nrow
tabelaGovedo <- data.frame()
for (chip in c("50Kv01", "50Kv02", "GGPv02", "GGPv03", "GGPv04", "HD", "HDv02", "IDBv03")) {
  chipImiss <- read.table(paste0(chip, "/PLINK_MERGED_", chip, ".imiss"), header=TRUE)
  chipImiss$Chip <- chip
  tabelaGovedo <- rbind(tabelaGovedo, chipImiss[,c("IID", "F_MISS", "Chip")])
}

nrow(indsRjava) == sum(indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL %in% tabelaGovedo$IID)
indsRjava[!(indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL %in% tabelaGovedo$IID),]
miss <- tabelaGovedo[!(tabelaGovedo$IID %in% indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL),]
table(tabelaGovedo$Chip[!(tabelaGovedo$IID %in% indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL)])
table(indsRjava[!(indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL %in% tabelaGovedo$IID),"GEN_CHIP"])

sum(indsR50K$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL %in% tabelaGovedo$IID)

table(indsRjava$GEN_CHIP[indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL %in% tabelaGovedo[tabelaGovedo$Chip=="IDBv03", "IID"]])
sum(tabelaGovedo$Chip=="IDBv03")
length(unique(tabelaGovedo$IID[tabelaGovedo$Chip=="IDBv03"]))



#Zdej pa naredi tabelo
#iz goveda - indsRjava
#iz vseh genotipov - tabelaGovedo
#datumi - indsDate

sum(tabelaGovedo$IID %in% indsRjava$ZIV.DRZ_ORIG_ZIVAL..ZIV.STEV_ORIG_ZIVAL)
colnames(indsRjava) <- c("IID", "SEQ", "Chip", "GEN_DATUM")
TABELA <- merge(indsRjava, tabelaGovedo, by=c("IID", "Chip"), all=TRUE)
TABELAmiss <- TABELA[is.na(TABELA$GEN_DATUM),]
TABELAfull <- TABELA[!(is.na(TABELA$GEN_DATUM)),]
colnames(indsDate) <- c("IID", "GenDate")
TABELA <- merge(TABELA, indsDate, by="IID")

indsDate <- unique(indsDate)
sum(TABELAmiss$IID %in% indsDate$IID)
DateMiss <- (indsDate[indsDate$IID %in% TABELAmiss$IID,])
DATE <- merge(TABELAmiss, DateMiss, by="IID", all.y=TRUE)

TABELA$GenDate <- NA
TABELA <- rbind(TABELAfull, DATE)
write.csv(TABELA, "~/GENTABELA.csv", quote=FALSE, row.names = FALSE)



tabelaGovedo[tabelaGovedo$IID %in% TABELA[is.na(TABELA$F_MISS),"IID"],]

indsRjava[indsRjava$IID=="SI04513638",]
indsOld <- read.csv("/home/jana/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_12042016/Matija_Rigler_12042016/ALL_INDDATE.txt")
indsOld[!(indsOld$X0 %in% TABELA$IID),]
TABELA[TABELA$GEN_DATUM=='05.04.2015',]
