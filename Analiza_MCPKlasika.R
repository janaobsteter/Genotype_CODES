outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

MCP <- read.csv("/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Zbirnik rezultatov_28112017.csv", header=TRUE, colClasses = 'character')
KontroleReje <- read.csv('/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Rejci_kontrole.csv')
intersect(KontroleReje$Rejec, MCP$Rejec)
outersect(KontroleReje$Rejec, MCP$Rejec)


MCP <- merge(MCP, KontroleReje, by="Rejec", all.x=TRUE)
MCP$Datum <- as.Date(MCP$Datum, format="%d.%m.%Y")
if (length(unique(MCP$ID)) == nrow(MCP)) {print("No duplicates in the data")} 


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
#KraveLakt <- read.csv('/home/jana/Documents/F4F/Rezultati_MCPKlasiak/Prireja_Krave.csv')
#KraveLakt <- CredaPV
KraveLakt$ID <- gsub("SI", "", KraveLakt$ID)

#DIM
DIMs <- paste0("SELECT distinct  ziv.ZIV_ID_SEQ, ziv.stev_orig_zival,
  ziv.SP1_SIFRA_PASMA pasma,
  ziv.DAT_ROJSTVO,
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

write.table(Data, '/home/jana/Documents/F4F/Rezultati_MCPKlasiak/TabelaRezultati_01122017.csv', quote=FALSE, sep="\t", row.names=FALSE)
write.table(Data, '/home/jana/Genotipi/Genotipi_CODES//Rezultati_MCPKlasiak/TabelaRezultati_01122017.csv', quote=FALSE, sep="\t", row.names=FALSE)
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
