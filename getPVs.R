library(ROracle)

#genetski trendi pri rjavi pasmi na kodo lastnosti 164 - mleko kg

drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")


con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

poizvedba <- paste("SELECT ziv.ZIV_ID_SEQ, ziv.SIF_SPOL,pv.VREDNOST_12_PV, extract(year from pv.DAT_OCENA_PV)
FROM govedo.zivali ziv,
ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
AND pv.SIFRA_LAST         =164
AND ziv.DRZ_ORIG_ZIVAL    ='SI'")

tabela <- fetch(dbSendQuery(con,poizvedba))
means <- aggregate(tabela$VREDNOST_12_PV ~ tabela$`EXTRACT(YEARFROMPV.DAT_OCENA_PV)`, FUN=mean)



########################################################################3
#gPV
mleko <- read.table("/home/jana/gPV_F4F/Mleko_Pv.txt")
mlekoSI <- read.table("/home/jana/gPV_F4F/SI_mlekoPV.txt")
mlekoB <- read.table("/home/jana/gPV_F4F/PVBeljakovine.txt")
VisinaKriza <- read.table("/home/jana/gPV_F4F/PVVisinaKriza.txt")
SkocniSklep <- read.table("/home/jana/gPV_F4F/PVSkocniSklep.txt")
Dolgozivost <- read.table("/home/jana/gPV_F4F/PVDolgozivost.txt")

ids <- c()
for (name in row.names(mleko)) {
  ids <- c(ids, substring(name,12, 19))
}

mleko <- Dolgozivost
mleko <- SkocniSklep
mleko <- mlekoB
mlekoSI <- mlekoSI[,c(8,9,10,11)]
mleko <- mleko[,c(1,4,10,11,12,13)]
VisinaKriza <- VisinaKriza[,c(1,3,9,10,11,12)]
mleko <- mleko[,c(1,3,9,10,11,12)]
mleko <- mleko[,c(1,2,9,10,11)]
mleko <- mleko[,c(1,9,10,11,12)]
mleko <- VisinaKriza
colnames(mleko) <- c("Lastnost", "Enota", "gPV", "gPV_12", "Rang", "Rang%")
colnames(mleko) <- c("Lastnost", "Enota", "gPV_12", "Rang", "Rang%")
colnames(mleko) <- c("Lastnost",  "gPV", "gPV_12", "Rang", "Rang%")
colnames(mlekoSI) <- c("gPV", "gPV_12", "Rang", "Rang%")
mleko$ID <- ids
mlekoSI$ID <- ids


library(ROracle)
drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")
con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

poizvedba <- paste("select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03'")

F4F <- fetch(dbSendQuery(con,poizvedba))

#ali imaš ID-je
length(intersect(ids, F4F$STEV_ORIG_ZIVAL)) #IMAŠ VSE!
length(intersect(mleko$ID, F4F$STEV_ORIG_ZIVAL)) #IMAŠ VSE!

PVp <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL,pv.VREDNOST_12_PV, pv.DAT_OCENA_PV
FROM govedo.zivali ziv,
  ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
AND pv.SIFRA_LAST         =169
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")

PVps <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL,pv.VREDNOST_12_PV, extract(year from pv.DAT_OCENA_PV) leto_obracuna
FROM govedo.zivali ziv,
  ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
AND pv.SIFRA_LAST         =100
AND ziv.STEV_ORIG_ZIVAL  in (select ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and gen.GEN_CHIP='IDBv03')")

pvs <- fetch(dbSendQuery(con, PVp))
pvSI <- fetch(dbSendQuery(con, PVps))
mlekoPV <- merge(mleko, pvs, by="ID", all.y=T)
mlekoPV_SI <- merge(mlekoSI, pvSI, by="ID", all.y=T)

cor(mlekoPV$gPV_12, mlekoPV$VREDNOST_12_PV)
plot(mlekoPV$gPV_12 ~ mlekoPV$VREDNOST_12_PV)
cor(mlekoPV_SI$gPV_12, mlekoPV_SI$VREDNOST_12_PV)
plot(mlekoPV_SI$gPV_12 ~ mlekoPV_SI$VREDNOST_12_PV)

mlekoPV_B <- mlekoPV
write.table(mlekoPV_B, "PVzivali_beljakovine.csv", quote=FALSE, row.names = FALSE)
write.table(mlekoPV, "PVzivali_mlekoKg.csv", quote=FALSE, row.names = FALSE)
hist(mlekoPV$gPV_12)
hist(mlekoPV_B$VREDNOST_12_PV)

#beljakovine
par(mfrow = c(1,2))
hist(mlekoPV_B$VREDNOST_12_PV, xlab = "Klasične PV 12", ylab = "Frekvenca", main =  "Beljakovine [kg]")
hist(mlekoPV_B$gPV_12, xlab = "Genomske PV 12", ylab = "Frekvenca", main =  "Beljakovine [kg]")
#mlekoKg
hist(mlekoPV$VREDNOST_12_PV, xlab = "Klasične PV 12", ylab = "Frekvenca", main =  "Beljakovine [kg]")
hist(mlekoPV$gPV_12, xlab = "Genomske PV 12", ylab = "Frekvenca", main =  "Beljakovine [kg]")
#Visina kriza
hist(mlekoPV$VREDNOST_12_PV, xlab = "Klasične PV 12", ylab = "Frekvenca", main =  "Višina križa [cm]")
hist(mlekoPV$gPV_12, xlab = "Genomske PV 12", ylab = "Frekvenca", main =  "Višina križa  [cm]")
#Skočni sklep
hist(mlekoPV$VREDNOST_12_PV, xlab = "Klasične PV 12", ylab = "Frekvenca", main =  "Skočni sklep")
hist(mlekoPV$gPV_12, xlab = "Genomske PV 12", ylab = "Frekvenca", main =  "Skočni sklep")


#Anželak
Anz <- mlekoPV[mlekoPV$CRE_SIFRA_CREDA==143,]
AnzSI <- mlekoPV_SI[mlekoPV_SI$CRE_SIFRA_CREDA==143,]
Anz <- Anz[order(Anz$VREDNOST_12_PV),]
Anz$RangPV <- 1:23
Anz$RangPV_ID <- Anz$ID
Anz <- Anz[order(Anz$gPV_12),]
Anz$RanggPV <- 1:23
Anz$RanggPV_ID <- Anz$ID

cor(Anz$VREDNOST_12_PV, Anz$gPV_12)
cor(Anz$RanggPV, Anz$RangPV)
cor(AnzSI$VREDNOST_12_PV, AnzSI$gPV_12)
cor(Anz$RanggPV, Anz$RangPV)
plot(Anz$RanggPV, Anz$RangPV)
plot(Anz$gPV_12, Anz$VREDNOST_12_PV)
plot(AnzSI$gPV_12, AnzSI$VREDNOST_12_PV)
#Anz$RanggPV <- as.numeric(Anz$RanggPV)
#Anz$RangPV <- as.numeric(Anz$RangPV)
plotRanks(Anz$RangPV, Anz$RanggPV)

plotRanks <- function(a, b, labels.offset=0.1, arrow.len=0.1)
{
  old.par <- par(mar=c(1,1,1,1))
  
  # Find the length of the vectors
  len.1 <- length(a)
  len.2 <- length(b)
  
  # Plot two columns of equidistant points
  plot(rep(1, len.1), 1:len.1, pch=20, cex=0.8, 
       xlim=c(0, 3), ylim=c(0, max(len.1, len.2)),
       axes=FALSE, xlab="", ylab="") # Remove axes and labels
  points(rep(2, len.2), 1:len.2, pch=20, cex=0.8)
  
  # Put labels next to each observation
  text(rep(1-labels.offset, len.1), 1:len.1, a)
  text(rep(2+labels.offset, len.2), 1:len.2, b)
  
  # Now we need to map where the elements of a are in b
  # We use the match function for this job
  a.to.b <- match(a, b)
  
  # Now we can draw arrows from the first column to the second
  arrows(rep(1.02, len.1), 1:len.1, rep(1.98, len.2), a.to.b, 
         length=arrow.len, angle=20)
  par(old.par)
}
