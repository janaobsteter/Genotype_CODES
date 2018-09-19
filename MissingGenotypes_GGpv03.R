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



##############################################
#combine tables for all genotyped BSW animals
##########################################
gen1 <- read.csv("~/Genotipi/Genotipi_CODES/TABELA_GenotipiziraneZivali_10082018.csv")
gen2 <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_23082018/08082018GovedoInd.csv")
gen3 <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_23082018/10082018GovedoInd.csv")
gen4 <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_04092018/14082018GovedoInd.csv", as.is=TRUE)
FMISSv03  <- read.table("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/GGPv03/PLINK_MERGED_GGPv03.imiss", header=TRUE)

gen2$GenoDate <- "08082018"
gen2$GEN_DATUM <- format(as.Date(paste(substr(gen2$GenoDate, 1, 2), substr(gen2$GenoDate, 3,4), substr(gen2$GenoDate, 5,8), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
gen3$GEN_DATUM <- format(as.Date(paste(substr(gen3$GenoDate, 1, 2), substr(gen3$GenoDate, 3,4), substr(gen3$GenoDate, 5,8), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
gen4$GEN_DATUM <- format(as.Date(paste(substr(gen4$GenoDate, 1, 2), substr(gen4$GenoDate, 3,4), substr(gen4$GenoDate, 5,8), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
colnames(gen1)[1] <- "ID"

GEN <- rbind(gen1, gen2[,c(1,4,7,6)])
GEN <- rbind(GEN, gen3[,c(1,4,7,6)])
GEN <- rbind(GEN, gen4[,c(1,4,7,6)])



#nekaj ni v redu z datumi
gendate <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/ALL_INDDATE.txt")
table(gendate$Chip)
gendate$Chip <- as.character(gendate$Chip)
genDate <- unique(gendate[gendate$Chip %in% c("50Kv01", "50Kv02", "GGPv02", "GGPv03", "GGPv04", "HD", "HDv02", "IDBv03"),1:3])
genDate <- genDate[!(genDate$Datum=="11072017"),]
nrow(genDate)
nrow(unique(genDate))
nrow(unique(genDate[genDate$Chip=="IDBv03",]))
nrow((genDate[genDate$Chip=="IDBv03",]))


genDate$Datum[genDate$Datum=="3112017"] <- "03112017"
genDate$Datum[genDate$Datum=="8082018"] <- "08082018"
table(genDate$Datum)

govedo <- read.csv("~/Genotipi/Genotipi_CODES/Govedo_GenotipiziraneZivali_04092018")
colnames(govedo)[1] <- "ID"
colnames(genDate)[1] <- "ID"

#those missing from govedo table
m <- genDate[!(genDate$ID %in% govedo$ID), ]
write.csv(m, "Tabela_missing.csv")
m <- read.csv("Tabela_missing.csv")[,-1]
m$Datum[m$Datum=="8082018"] <- "08082018"

#popravi datume
mi <- m[m$Chip=="IDBv03",]
mg <- m[m$Chip=="GGPv04",]
mi$GEN_DATUM <- format(as.Date(paste(substr(mi$Datum, 1, 2), substr(mi$Datum, 3,4), substr(mi$Datum, 5,8), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
mg$GEN_DATUM <- format(as.Date(paste(substr(mg$Datum,7,8), substr(mg$Datum, 5,6), substr(mg$Datum, 1,4), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
M <- rbind(mi, mg)

#tiste, ki ne manjkajo - popravi datume
ma <- genDate[genDate$ID %in% govedo$ID, ]
table(ma$Datum, ma$Chip)
mai <- ma[ma$Chip=="IDBv03",]
mai$GEN_DATUM <- format(as.Date(paste(substr(mai$Datum, 1, 2), substr(mai$Datum, 3,4), substr(mai$Datum, 5,8), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
mag <- ma[ma$Chip!="IDBv03",]
mag$GEN_DATUM <- format(as.Date(paste(substr(mag$Datum,7,8), substr(mag$Datum, 5,6), substr(mag$Datum, 1,4), sep="-"), format="%d-%m-%Y"), "%d.%m.%Y")
MA <- rbind(mai, mag)

#to so vse, ki sem jih pobrala 
MA[!(MA$ID %in% govedo$ID),]#to je check

g <- govedo[!(govedo$ID %in% genDate$ID),]
nrow(g)
table(g$GEN_CHIP)
head(g)


g <- g[,c(1,3,4)]
colnames(g) <- colnames(Merge)


h <- govedo[govedo$ID %in% MA$ID,]
h <- h[,c(1,3,4)]
mtest <- MA[,c(1,3,4)]
colnames(h) <- colnames(mtest)
Merge <- unique(merge(mtest, h, by=c("ID", "Chip", "GEN_DATUM"), all.x=TRUE))
nrow(Merge)

Mf <- M[,c(1,3,4)]

FINAL <- rbind(Merge, g)
FINAL <- rbind(FINAL, Mf)
FINAL <- unique(FINAL)
table(FINAL$Chip)

#preveri, koliko imaš duplicated ID - DATUM
FINALdva <- FINAL[,c("ID", "Chip")]
sum(duplicated(FINALdva))
ddd <- FINALdva[which(duplicated(FINALdva)),]
ddd1 <- FINALdva[which(duplicated(FINALdva, fromLast = TRUE)),]
DDD <- rbind(ddd, ddd1)
DDD[order(DDD$ID),]
DDDIDs <- unique(DDD$ID)

FMISS <- read.csv("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/FMISS_All.csv")
FMISS$IID[FMISS$IID=="SI24951463MAGNET"] <- "SI24951463"
nrow(FMISS[FMISS$IID %in% FINAL$ID,])
nrow(FMISS[!(FMISS$IID %in% FINAL$ID),])
FINAL[!(FINAL$ID %in% FMISS$IID),]
table(FINAL$Chip[!(FINAL$ID %in% FMISS$IID)])
# 
FINAL$GEN_DATUM[FINAL$GEN_DATUM=="31.12.17"] <- "31.12.2017"

FINALbrezDupl <- FINAL[!(FINAL$ID %in% DDDIDs),]
sum(duplicated(FINALbrezDupl$ID))
d1 <- FINALbrezDupl[which(duplicated(FINALbrezDupl$ID)),]
d2 <- FINALbrezDupl[which(duplicated(FINALbrezDupl$ID, fromLast = TRUE)),]
D <- rbind(d1, d2)
D[(order(D$ID)),]
D[which(duplicated(c("ID", "Chip"))),]

#to je zdj prazno
#FINALd1 <- FINAL[duplicated(FINAL[,c("ID", "Chip")], fromLast = TRUE),]
#FINALd2 <- FINAL[duplicated(FINAL[,c("ID", "Chip")], fromLast = FALSE),]

#FINAL <- FINAL[!(duplicated(FINAL[,c("ID", "Chip")], fromLast = TRUE)),]
#FINAL <- FINAL[!(duplicated(FINAL[,c("ID", "Chip")], fromLast = FALSE)),]

colnames(FMISS)[1] <- "ID"

FINALF <- unique(merge(FINALbrezDupl, FMISS, by=c("ID", "Chip"), all.x=TRUE))
sum(duplicated(FINALF[,c("ID","Chip")]))

write.table(FINALF, "~/Genotipi/Genotipi_DATA/Genotipi_latest/GenotipiziraneZivali_PlusFMISS_11092018.csv", quote=FALSE, row.names = FALSE)
#FINALF <- read.csv("~/Genotipi/Genotipi_DATA/Genotipi_latest/GenotipiziraneZivali_PlusFMISS_04092018.csv", sep=" ")

#Te morajo biti dopolnjene na roke! kER JE BILA ISTA ŽIVAL dvakrat genotipzirana na istem čipu! Ne morem potegniti iz merge files
FINALdupl <- rbind(FINALd1, FINALd2)
FINALdupl$GEN_DATUM[FINALdupl$GEN_DATUM=="31.12.2017"] <- "03.11.2017"
FINALdupl <- FINALdupl[order(FINALdupl$ID),]
FINALdupl$FMISS <- NA
write.table(FINALdupl, "~/Genotipi/Genotipi_DATA/Genotipi_latest/GenotipiziraneZivali_ManjkaFMISS_04092018.csv", quote=FALSE, row.names = FALSE)
