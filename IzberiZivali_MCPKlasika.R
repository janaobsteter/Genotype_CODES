#odbira živali za vzorčenje mleka F4F
library(ROracle)

#genetski trendi pri rjavi pasmi na kodo lastnosti 164 - mleko kg

drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")


con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)
PVmleko <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL,pv.VREDNOST_12_PV, pv.DAT_OCENA_PV
FROM govedo.zivali ziv,
  ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
and pv.sif_izv_vr_last=1
AND pv.SIFRA_LAST         =164
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")

PVbelj <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL,pv.VREDNOST_12_PV, extract(year from pv.DAT_OCENA_PV) leto_obracuna
FROM govedo.zivali ziv,
  ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
and pv.sif_izv_vr_last=1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
AND pv.SIFRA_LAST         =168
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")

PVmasc <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL,pv.VREDNOST_12_PV, extract(year from pv.DAT_OCENA_PV) leto_obracuna
FROM govedo.zivali ziv,
  ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
and pv.sif_izv_vr_last=1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
AND pv.SIFRA_LAST         =165
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")

mleko <- fetch(dbSendQuery(con, PVmleko))
belj <- fetch(dbSendQuery(con, PVbelj))
masc <- fetch(dbSendQuery(con, PVmasc))
masc <- merge(mleko, pvs, by="ID", all.y=T)

#naredi summary in najdi kmetije z največjo varianco
mlekoA <- aggregate(mleko$VREDNOST_12_PV ~mleko$CRE_SIFRA_CREDA, FUN=function(x) {sd(x)})
colnames(mlekoA) <- c("Creda", "MlekoSD")

beljA <- aggregate(belj$VREDNOST_12_PV ~belj$CRE_SIFRA_CREDA, FUN=function(x) {sd(x)})
colnames(beljA) <- c("Creda", "BeljSD")

mascA <- aggregate(masc$VREDNOST_12_PV ~masc$CRE_SIFRA_CREDA, FUN=function(x) {sd(x)})
colnames(mascA) <- c("Creda", "MascSD")

SD <- Reduce(function(x, y) merge(x, y, all=TRUE), list(mlekoA, beljA, mascA))

#poglej še po zaporedni laktaciji in stadiju laktacije
ZapLakt <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL, count(distinct tel.TEL_ID_SEQ) StLakt
FROM govedo.zivali ziv, GOVEDO.TELITVE tel  
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =tel.tel_ziv_id_seq
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')
group by ziv.STEV_ORIG_ZIVAL, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL")

zapLaktAll <- fetch(dbSendQuery(con, ZapLakt))
zapLakt <- aggregate(zapLaktAll$STLAKT ~zapLaktAll$CRE_SIFRA_CREDA, FUN=function(x) {var(x)})
colnames(zapLakt) <- c("Creda", "ZaplaktSD")

SD <- merge(SD, zapLakt, by="Creda")

#poglej stadij laktacije
DIMs <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA Creda, ziv.SIF_SPOL, tel.DAT_TELITEV, lak.DAT_PRESUSITEV
FROM govedo.zivali ziv, GOVEDO.TELITVE tel  , GOVEDO.LAKTACIJE lak,
(SELECT ziv.ZIV_ID_SEQ ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL, count(distinct tel.TEL_ID_SEQ) maxTel
FROM govedo.zivali ziv, GOVEDO.TELITVE tel  
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =tel.tel_ziv_id_seq
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')
group by ziv.ZIV_ID_SEQ, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL) stTel
WHERE ziv.SP1_SIFRA_PASMA = 1
and tel.TEL_ZIV_ID_SEQ = ziv.ZIV_ID_SEQ
AND tel.TEL_ZIV_ID_SEQ        =stTel.ID
and tel.ZAP_TELITEV=stTel.maxTel
and lak.LAK_TEL_ID_SEQ=tel.TEL_ID_SEQ
and lak.DAT_PRESUSITEV is NULL
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")


DIM <- fetch(dbSendQuery(con, DIMs))
DIM$DAT_TELITEV <- as.Date(DIM$DAT_TELITEV, format="%y-%m-%d")
DIM$DIM <- difftime(Sys.Date(), DIM$DAT_TELITEV)
DIM$DIM <- as.numeric(DIM$DIM)
summary(DIM$DIM)
#DIM <- DIM[-which(DIM$DIM > 720),] # izloči tiste nad 720 dni


DIMA <- aggregate(DIM$DIM ~ DIM$CREDA, FUN=function(x) {sd(x)})
colnames(DIMA) <- c("Creda", "DIMSD")
SD <- merge(SD, DIMA, by="Creda")

#Kazeini
kapaCSN <- read.csv('/home/jana/Documents/F4F/MlecniProteini/KappaCaseinGenotype_python.csv')
crede <- paste0("SELECT ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL
FROM govedo.zivali ziv
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")

Crede <- fetch(dbSendQuery(con, crede))
kapaCSN <- merge(kapaCSN, Crede, by="ID")
library(dplyr)
library(AggregateR)
kapaByHerd <- Aggregate(as.data.frame(kapaCSN$SKUPEN), by=kapaCSN$CRE_SIFRA_CREDA)[,c(1:7)]
colnames(kapaByHerd) <- c("Creda", "AA", "AB", "AC", "BB", "BC", "BE")
SD <- merge(SD, kapaByHerd, by="Creda")  


#Pridobi še mlečnosti, izplen beljakovin in maščob

prireja <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID,ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL celID, ziv.IME_ZIVAL, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL, round(avg(lak.KG_MLEKO_305),0) mlekoAvg,round( avg(lak.KG_BELJAK_305),0) beljAvg, round(avg(lak.KG_MAST_305),0) mascAvg
FROM govedo.zivali ziv, GOVEDO.TELITVE tel , GOVEDO.LAKTACIJE lak
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =tel.tel_ziv_id_seq
and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')
group by ziv.STEV_ORIG_ZIVAL, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL, ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL, ziv.IME_ZIVAL")


Prireja <- fetch(dbSendQuery(con, prireja))
Prireja$OdsBelj <- round(Prireja$BELJAVG / Prireja$MLEKOAVG * 100,2)
Prireja$OdsMsc <- round(Prireja$MASCAVG / Prireja$MLEKOAVG * 100, 2)


MlekoA <- aggregate(Prireja$MLEKOAVG ~ Prireja$CRE_SIFRA_CREDA, FUN=function(x) {sd(x)})
colnames(MlekoA) <- c("Creda", "MLEKOSD")
SD <- merge(SD, MlekoA, by="Creda")
BeljA <- aggregate(Prireja$OdsBelj ~ Prireja$CRE_SIFRA_CREDA, FUN=function(x) {sd(x)})
colnames(BeljA) <- c("Creda", "BeljOdsSD")
SD <- merge(SD, BeljA, by="Creda")
#SEDAJ POBERI TISTE Z NAJVEČJIMI SD_JI i VARIANCAMI za vse kriterije
#) po PV
PVMlekoTop <- SD[order(-SD$MlekoSD),"Creda"][1:10]
PVBeljTop <- SD[order(-SD$BeljSD),"Creda"][1:10]
PVMascTop <- SD[order(-SD$MascSD),"Creda"][1:10]
Inter <- intersect(intersect(PVMlekoTop, PVBeljTop), PVMascTop)

ZapLaktTop <- SD[order(-SD$ZaplaktSD),"Creda"]#[1:10]
ZapLaktTop <- data.frame(ORderZapLakt = SD[order(-SD$ZaplaktSD),"Creda"])#[1:10]
DIMTop <- SD[order(-SD$DIMSD),"Creda"][1:10]
intersect(PVBeljTop, ZapLaktTop)

DIMOrder <- data.frame(ORderDIM =SD[order(-SD$DIMSD),"Creda"])
MlekoOrder <- data.frame(OrderMleko =SD[order(-SD$MLEKOSD),"Creda"])
BeljOrder <- data.frame(OrderBeljOds =SD[order(-SD$BeljOdsSD),"Creda"])

Orders <- cbind(ZapLaktTop, MlekoOrder, DIMOrder,BeljOrder)
write.table(Orders, '/home/jana/Genotipi/Genotipi_CODES/F4FKraveCrede_ORDERPrirejaDIM.csv', quote=FALSE, row.names=FALSE)

PVs <- merge(mleko[,c(1,4)], belj[,c(1,4)], by="ID")
PVs <- merge(PVs, masc[,c(1,4)], by="ID")
colnames(PVs) <- c("ID", "PVmleko", "PVbelj", "PVmasc")
kazeini <- kapaCSN[,c(1,8)]
 kazeini$ID <- gsub("SI", "", kazeini$ID)
#Zaporedna laktacija
Izbira1 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[1],]
Izbira1 <- merge(Izbira1, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira1 <- merge(Izbira1, PVs, by="ID", all.x=TRUE)
Izbira1 <- merge(Izbira1, kazeini, by="ID", all.x=TRUE)
Izbira1 <- merge(Izbira1, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)


Izbira2 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[2],]
Izbira2 <- merge(Izbira2, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira2 <- merge(Izbira2, PVs, by="ID", all.x=TRUE)
Izbira2 <- merge(Izbira2, kazeini, by="ID", all.x=TRUE)
Izbira2 <- merge(Izbira2, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)


Izbira3 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[3],]
Izbira3 <- merge(Izbira3, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira3 <- merge(Izbira3, PVs, by="ID", all.x=TRUE)
Izbira3 <- merge(Izbira3, kazeini, by="ID", all.x=TRUE)
Izbira3 <- merge(Izbira3, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)

Izbira4 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[4],]
Izbira4 <- merge(Izbira4, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira4 <- merge(Izbira4, PVs, by="ID", all.x=TRUE)
Izbira4 <- merge(Izbira4, kazeini, by="ID", all.x=TRUE)
Izbira4 <- merge(Izbira4, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)


Izbira5 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[5],]
Izbira5 <- merge(Izbira5, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira5 <- merge(Izbira5, PVs, by="ID", all.x=TRUE)
Izbira5 <- merge(Izbira5, kazeini, by="ID", all.x=TRUE)
Izbira5 <- merge(Izbira5, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)

Izbira6 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[6],]
Izbira6 <- merge(Izbira6, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira6 <- merge(Izbira6, PVs, by="ID", all.x=TRUE)
Izbira6 <- merge(Izbira6, kazeini, by="ID", all.x=TRUE)
Izbira6 <- merge(Izbira6, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)


Izbira7 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[7],]
Izbira7 <- merge(Izbira7, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira7 <- merge(Izbira7, PVs, by="ID", all.x=TRUE)
Izbira7 <- merge(Izbira7, kazeini, by="ID", all.x=TRUE)
Izbira7 <- merge(Izbira7, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)


Izbira8 <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==ZapLaktTop[8],]
Izbira8 <- merge(Izbira8, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
Izbira8 <- merge(Izbira8, PVs, by="ID", all.x=TRUE)
Izbira8 <- merge(Izbira8, kazeini, by="ID", all.x=TRUE)
Izbira8 <- merge(Izbira8, Prireja[,c(1,4,5,6)], by="ID", all.x=TRUE)


Top8Cred <- rbind(Izbira1, Izbira2, Izbira3, Izbira4, Izbira5, Izbira6, Izbira7, Izbira8)
write.csv(Top8Cred, '/home/jana/Documents/F4F/Top8Cred_MCPklasika.csv', quote=FALSE, row.names=FALSE)
as.vector(unique(Top8Cred$CRE_SIFRA_CREDA))

AllCrede <- matrix(ncol=14)
colnames(AllCrede) <- colnames(credaPV(86))
for (creda in unique(Crede$CRE_SIFRA_CREDA)) {
  AllCrede <- rbind(AllCrede, credaPV(creda))
}

write.table(AllCrede, '/home/jana/Genotipi/Genotipi_CODES/F4FKraveCrede_PrirejaDIM.csv', quote=FALSE, row.names=FALSE)
credaPV <- function (creda) {
  CredaPV <- zapLaktAll[zapLaktAll$CRE_SIFRA_CREDA==creda,]
  CredaPV <- merge(CredaPV, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
  CredaPV <- merge(CredaPV, PVs, by="ID", all.x=TRUE)
  CredaPV <- merge(CredaPV, kazeini, by="ID", all.x=TRUE)
  CredaPV <- merge(CredaPV, Prireja[,c(1,2,3,6,7,8,9,10)], by="ID", all.x=TRUE)
  CredaPV <- CredaPV[,c(10,11,4,5,6,7,8,9,12,13,14,15,16)]
  colnames(CredaPV) <- c("ID", "IME", "Zaporedna lakt", "DIM", "PVMleko", "PVBelj", "PVMasc", "Kapa-kazein", "Povprecna st. lakt.", "Povprecne belj 305", "Povprecne masc 305", "Ods Belj", "Ods Masc")
  return(CredaPV)
}
write.table(credaPV(4950), '/home/jana/Documents/F4F/KlasikaAnaliza/Kuhar.csv', quote=FALSE, row.names=FALSE, sep=",")


CredaPV <- zapLaktAll
CredaPV <- merge(CredaPV, DIM[,c("ID", "DIM")], by="ID", all.x=TRUE)
CredaPV <- merge(CredaPV, PVs, by="ID", all.x=TRUE)
CredaPV <- merge(CredaPV, kazeini, by="ID", all.x=TRUE)
CredaPV <- merge(CredaPV, Prireja[,c(1,2,3,6,7,8,9,10)], by="ID", all.x=TRUE)
CredaPV <- CredaPV[,c(10,11,4,5,6,7,8,9,12,13,14,15,16)]
colnames(CredaPV) <- c("ID", "IME", "Zaporedna lakt", "DIM", "PVMleko", "PVBelj", "PVMasc", "Kapa-kazein", "Povprecna st. lakt.", "Povprecne belj 305", "Povprecne masc 305", "Ods Belj", "Ods Masc")
write.table(CredaPV, '/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Prireja_Krave.csv', quote=FALSE, row.names=FALSE)

