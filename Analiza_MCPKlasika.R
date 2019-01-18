outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

MCP <- read.csv("/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Skupni_rezultati_popravljeni_11122017.csv", header=TRUE)
#dodaj SI pred ID-je. Če se začne s 4, dodaj 0 na začetek
IDs <- c()
for (id in MCP$ID) {
  name <- ifelse( nchar(id) == 7, paste0("SI0", id), paste0("SI", id)) #substr(id, 1, 1)==4 &&
  IDs <- c(IDs, name)
}
MCP$ID <- IDs

MCP2 <- read.csv("/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Rezultati april-maj 2018.csv", header=TRUE)
MCP2$Rejec <- as.character(MCP2$Rejec)
MCP2$Rejec[which(MCP2$Rejec=="Leskovšek")] <- "Vasle"
MCP2$Rejec <- as.factor(MCP2$Rejec)
MCP$Season <- "Zima"
MCP2$Season <- "Pomlad"
MCP2$KappaKazein <- NA
summary(MCP$r)
summary(MCP2$r)
par(mfrow=c(2,1))
hist(MCP$a30)
hist(MCP2$a30)
hist(as.numeric(MCP$r.s.))
hist(as.numeric(MCP2$r.s.))


MCPA <- rbind(MCP, MCP2)
KontroleReje <- read.csv('/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Rejci_kontrole.csv', header=TRUE)
KontroleReje$Datum <- as.Date(KontroleReje$Datum, format="%d.%m.%Y")
MCPA <- merge(MCPA, KontroleReje, by=c("Rejec", "Season"))
hist(MCPA$r.s.)
hist(MCPA$a30.mm.)

MCPA$r.s.[MCPA$r.s. == 0] <- NA
MCPA$a30.mm.[MCPA$a30.mm. == 0] <- NA
MCPA$k20.s.[MCPA$k20.s. == 0] <- NA
MCPA$STEV_ORIG <- gsub("SI", "", MCPA$ID)

#reda in kappa kazein
MCPA$IME[MCPA$IME =="Sivka II"] <- "Sivka"
MCPA$ID[MCPA$ID == "SI97158890"] <- "SI94158890"
MCPA$STEV_ORIG[MCPA$STEV_ORIG == 97158890] <- 94158890
MCPA$ID[MCPA$ID == "SI63817455"] <- "SI63817445"
MCPA$STEV_ORIG[MCPA$STEV_ORIG == 63817455] <- 63817445
kk <- read.csv("~/Documents/F4F/MlecniProteini/KappaCaseinGenotype_python.csv")
length(intersect(MCPA$ID, kk$ID))
length(unique(MCPA$ID))
length(unique(kk$ID))
unique(MCPA$ID[!(MCPA$ID %in% kk$ID)])
length(unique(MCPA$ID[!(MCPA$ID %in% kk$ID)]))
MCPA <- merge(MCPA, kk[,c("ID", "SKUPEN")], by="ID", all.x=TRUE)
breed <- unique(read.csv("~/Documents/F4F/AnimalBreed.csv", colClasses = c(STEV_ORIG = "character")))
MCPA$STEV_ORIG[!(MCPA$STEV_ORIG %in% breed$STEV_ORIG)]
outersect(MCPA$STEV_ORIG, breed$STEV_ORIG)
MCPA <- merge(MCPA, breed, by="STEV_ORIG", all.x = TRUE)

tel <- read.csv("~/Documents/F4F/Rezultati_MCPKlasiak/Telitve.csv")
IDs <- c()

for (id in tel$STEV_ORIG) {
  name <- ifelse( nchar(id) == 7, paste0("0", id), id) #substr(id, 1, 1)==4 &&
  IDs <- c(IDs, name)
}
tel$STEV_ORIG <- as.character(IDs)

length(unique(tel$STEV_ORIG))
length(unique(MCPA$STEV_ORIG))
tel$DAT_TELITEV <- as.Date(tel$DAT_TELITEV, format="%d-%m-%y")

nrow(MCPA)
length(unique(MCPA$STEV_ORIG[(MCPA$STEV_ORIG %in% tel$STEV_ORIG)]))
unique(MCPA$STEV_ORIG[!(MCPA$STEV_ORIG %in% tel$STEV_ORIG)])
MCPA$STEV_ORIG[!(MCPA$STEV_ORIG %in% tel$STEV_ORIG)]

tel <- merge(MCPA, tel, by="STEV_ORIG")
tel <- tel[tel$Datum > tel$DAT_TELITEV,]
tel$dnPoLakt <- difftime(tel$Datum, tel$DAT_TELITEV, unit="days")
tel <- unique(tel)

#izberi samo ustrezne laktacije
tel$DatID <- paste0(tel$STEV_ORIG, tel$Datum)
TEL <- data.frame()
for (datID in unique(tel$DatID)) {
  tmp <- tel[tel$DatID==datID,]
  TEL <- rbind(TEL, tmp[order(tmp$dnPoLakt),][1,])
}

MCPA <- merge(MCPA, TEL[,c("STEV_ORIG", "Datum", "DAT_TELITEV", "dnPoLakt")], by=c("STEV_ORIG", "Datum"))
MCPA$StLakt <- ifelse(MCPA$dnPoLakt< 100, 1, ifelse(100 <= MCPA$dnPoLakt & MCPA$dnPoLakt < 200, 2, ifelse(200 <= MCPA$dnPoLakt & MCPA$dnPoLakt < 300, 3, 4)))
aggregate(MCPA$a30.mm. ~ MCPA$StLakt, FUN="mean")

write.csv(unique(MCPA), "~/Documents/F4F/Rezultati_MCPKlasiak/SkupniRezultati_Klasika.csv", quote=FALSE, row.names=FALSE)


#analiza vplivov
MCPA$Season <- as.factor(MCPA$Season)
MCPA$StLakt <- as.factor(MCPA$StLakt)
MCPA$Cas <- as.factor(MCPA$Cas)
MCPA$KappaKazein <- as.factor(MCPA$KappaKazein)
MCPA$SP1_SIFRA_PASMA <- as.factor(MCPA$SP1_SIFRA_PASMA)
MCPA$a30.mm. <- as.numeric(MCPA$a30.mm.)
MCPA$Cas[MCPA$Cas=="zjvečer"] <- "zvečer"
MCPA$Cas[MCPA$Cas=="zvecer"] <- "zvečer"
MCPA$Cas[MCPA$Cas=="/"] <-NA
MCPA$KappaKazein[MCPA$KappaKazein==""] <-NA
model1 <- lm(MCPA$a30.mm. ~ MCPA$Season  + MCPA$Cas + MCPA$proteini + MCPA$Mascoba + MCPA$SP1_SIFRA_PASMA + MCPA$StLakt)
plot(model1)
summary(model1)
anova(model1)

intersect(KontroleReje$Rejec, MCP$Rejec)
intersect(KontroleReje$Rejec, MCPA$Rejec)
length(intersect(KontroleReje$Rejec, MCPA$Rejec))
length(intersect(KontroleReje$Rejec, MCP$Rejec))
length(unique(KontroleReje$Rejec))
length(unique(MCP$Rejec))
unique(MCP$Rejec)
outersect(KontroleReje$Rejec, MCP$Rejec)


MCP <- merge(MCP, KontroleReje, by="Rejec", all.x=TRUE)
MCP$Datum <- as.Date(MCP$Datum, format="%d.%m.%Y")
if (length(unique(MCP$ID)) == nrow(MCP)) {print("No duplicates in the data")} 

#primerjaj pomlad-zimo
library(tidyr)
MCPA$MonthKl <- format(MCPA$Datum, "%m-%y")
zima <- MCPA[MCPA$Season=="Zima", ]
colnames(zima)[3:23] <- paste0(colnames(zima)[3:23], "_Z")
pomlad <- MCPA[MCPA$Season=="Pomlad", ]
colnames(pomlad)[3:23] <- paste0(colnames(pomlad)[3:23], "_P")
compare <- merge(zima, pomlad, by=c("ID", "STEV_ORIG"))
cor(compare$a30.mm._P, compare$a30.mm._Z, use = "pairwise.complete.obs")
cor(compare$k20.s._P, compare$k20.s._Z, use = "pairwise.complete.obs")
cor(compare$r.s._Z, compare$r.s._P, use = "pairwise.complete.obs")
plot(compare$a30.mm._Z, compare$a30.mm._P)

# #repeatability
# library(heritability)
# repa30 <- lm(MCPA$a30.mm. ~ MCPA$ID )
# repeatability(repa30)



#Tukaj pridobi podatke o teh živalih iz baze
drv <- dbDriver("Oracle")
# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")
con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

#Zaporedna Latkacija
ZapLakt <- paste0("SELECT ziv.STEV_ORIG_ZIVAL ID, ziv.CRE_SIFRA_CREDA, 
count(distinct tel.TEL_ID_SEQ) StLakt FROM govedo.zivali ziv, GOVEDO.TELITVE tel  
where ziv.ZIV_ID_SEQ(+)        =tel.tel_ziv_id_seq
AND ziv.STEV_ORIG_ZIVAL  in ('", paste(MCP$ID,collapse = "','"), "') group by ziv.STEV_ORIG_ZIVAL, ziv.CRE_SIFRA_CREDA, ziv.SIF_SPOL")

ZapL <- fetch(dbSendQuery(con,ZapLakt))
KraveLakt <- read.csv('/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Prireja_Krave.csv')
#KraveLakt <- CredaPV
KraveLakt$ID <- gsub("SI", "", KraveLakt$ID)

#DIM
DIMs <- paste0("SELECT distinct  ziv.ZIV_ID_SEQ, ziv.stev_orig_zival,
  ziv.SP1_SIFRA_PASMA pasma,
  extract (year from ziv.DAT_ROJSTVO) datRoj,
  mleko.mlekoPV,
  beljKg.beljKgPV,
  mascKg.mascKgPV,
  beljOds.beljOdsPV,
  mascOds.mascOdsPV,
  COUNT(DISTINCT tel.TEL_ID_SEQ) ZapTel,
  MAX(tel.DAT_TELITEV) DatTel
FROM  zivali ziv, GENOTIPIZIRANE_ZIVALI gen,
 telitve tel,
 (  select pv.PV_ZIV_ID_SEQ ZIV_ID_SEQ, pv.VREDNOST_12_PV mlekoPV from ARHIV.PLEMENSKE_VREDNOSTI pv where pv.SIFRA_LAST=164 and pv.SIF_IZV_VR_LAST=4) mleko,
 (  select pv.PV_ZIV_ID_SEQ ZIV_ID_SEQ, pv.VREDNOST_12_PV beljKGPV from ARHIV.PLEMENSKE_VREDNOSTI pv where pv.SIFRA_LAST=168 and pv.SIF_IZV_VR_LAST=4 ) beljKg,
 (  select pv.PV_ZIV_ID_SEQ ZIV_ID_SEQ, pv.VREDNOST_12_PV mascKgPV from ARHIV.PLEMENSKE_VREDNOSTI pv where pv.SIFRA_LAST=166 and pv.SIF_IZV_VR_LAST=4) mascKg,
 (  select pv.PV_ZIV_ID_SEQ ZIV_ID_SEQ, pv.VREDNOST_12_PV beljOdsPV from ARHIV.PLEMENSKE_VREDNOSTI pv where pv.SIFRA_LAST=167 and pv.SIF_IZV_VR_LAST=4) beljOds,
 (  select pv.PV_ZIV_ID_SEQ ZIV_ID_SEQ, pv.VREDNOST_12_PV mascOdsPV from ARHIV.PLEMENSKE_VREDNOSTI pv where pv.SIFRA_LAST=165 and pv.SIF_IZV_VR_LAST=4) mascOds
 WHERE ziv.ZIV_ID_SEQ    =tel.TEL_ZIV_ID_SEQ
AND ziv.stev_orig_zival in  ('", paste(MCP$ID,collapse = "','"), "') 
and mleko.ZIV_ID_SEQ(+) = ziv.ZIV_ID_SEQ
and beljKg.ZIV_ID_SEQ(+) = ziv.ZIV_ID_SEQ
and mascKg.ZIV_ID_SEQ(+) = ziv.ZIV_ID_SEQ
and beljOds.ZIV_ID_SEQ(+) = ziv.ZIV_ID_SEQ
and mascOds.ZIV_ID_SEQ(+) = ziv.ZIV_ID_SEQ
GROUP BY ziv.ZIV_ID_SEQ,
  ziv.SP1_SIFRA_PASMA, ziv.stev_orig_zival,
  ziv.DAT_ROJSTVO, mleko.mlekoPV, beljKg.beljKgPV,
  mascKg.mascKgPV,
  beljOds.beljOdsPV,
  mascOds.mascOdsPV")


DIM <- fetch(dbSendQuery(con,DIMs))
colnames(DIM)[2] <- "ID"
#NO MORE DUPLICATED IN 28112017 VERZIJI!
#duplicated IDs MCP Rodica
#Dupl <- MCP[MCP$ID %in% MCP$ID[duplicated(MCP$ID)],]
#Dupl <- Dupl[order(Dupl$ID),]
#write.table(Dupl, "/home/jana/Documents/F4F/Rezultati_MCPKlasiak/DuplicatedRecords.csv", quote=FALSE, row.names=FALSE)



#Tukaj združi MCP (Rodica) podatke s podatki o zaporedni laktaciji in prireji
Data <- merge(MCP, ZapL, by="ID")
Data <- merge(Data, DIM, by="ID")
#sedaj zračunaj DIM
Data$DATTEL <- as.Date(Data$DATTEL, format="%Y-%d-%m")
Data$Datum <- as.Date(Data$Datum, format="%Y-%d-%m")
Data$DIM <- round(as.numeric(gsub(" days", "", difftime(Data$Datum, Data$DATTEL)), 1))
#tu dodaj še plemenske vrednosti

#Dodaj še genotipe
genotipi <- read.csv("/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/AB/IDBv03/Genotipi_F4F/MonogenicGenotypes_Table.csv")
genotipi$ID <- gsub("SI", "", genotipi$ID)
length(intersect(Data$ID, genotipi$ID))
Data <- merge(Data, genotipi, by="ID", all.x=TRUE)
Data <- Data[Data$PASMA %in% 1:3,]
Data$PASMA <- gsub(3, "Holstein", Data$PASMA)
Data$PASMA <- gsub(1, "Rjava", Data$PASMA)
Data$PASMA <- gsub(2, "Lisasta", Data$PASMA)


#pridobi celotne ID-je
IDs <- paste0("SELECT ziv.drz_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL ID, ziv.STEV_ORIG_ZIVAL
FROM govedo.zivali ziv
where ziv.STEV_ORIG_ZIVAL  in ('", paste(MCP$ID,collapse = "','"), "')")

IDd <- fetch(dbSendQuery(con,IDs))
colnames(IDd) <- c("CompleteID", "ID")
length(intersect(IDd$ID, Data$ID)) == nrow(Data)

Data <- read.table('/home/jana/Documents/F4F/Rezultati_MCPKlasiak/TabelaRezultati_11122017.csv', sep="\t", header=TRUE)
Data <- merge(Data, IDd, by="ID")

write.table(Data, '/home/jana/Documents/F4F/Rezultati_MCPKlasiak/TabelaRezultati_11122017.csv', quote=FALSE, sep="\t", row.names=FALSE)
write.table(Data, '/home/jana/Genotipi/Genotipi_CODES//Rezultati_MCPKlasiak/TabelaRezultati_11122017.csv', quote=FALSE, sep="\t", row.names=FALSE)
###################################################################
###################################################################
nrow(MCP)
length(unique(MCP$ID))
length(intersect(MCP$ID, KraveLakt$ID))

Data <- merge(MCP, KraveLakt, by="ID")

Data$Rejec <- as.character(Data$Rejec)
Data$k20 <- as.numeric(Data$k20)
Data$r <- as.numeric(Data$r)
Data$a30 <- as.numeric(Data$a30)
Data$Mascoba <- as.numeric(Data$Mascoba)
Data1$Mascoba <- as.numeric(Data1$Mascoba)
Data1 <- Data[-(is.na(Data$Mascoba)),]
library(lme4)
model <- lmer(r ~ 1 + (1 | Rejec) + Mascoba + proteini, data=Data1, na.action=na.omit)
model <- lmer(r ~ 1 + (1 | Rejec) + Mascoba, data=Data1, na.action=na.omit )
model <- lmer(Data$r ~ 1+ (1 | Data$Rejec) )
summary(model)
