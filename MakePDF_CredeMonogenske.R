#script for producing PDFs for each table
#F4F projekt - odbria živali za genotipizacijo
#vsaka čreda posebna tabela
#primarni seznam - izbrane živali in seznam B - seznam nadomestnih živali v primeru, da je ena iz primarnega seznama izločena
setwd('/home/jana/Documents/F4F/OdbiraZivali/DodatnaOdbira_11052017/')
library(Hmisc)
library(knitr)


##  make my data
rejci <- read.csv('~/Documents/F4F/OdbiraZivali/RejciImena_15032017.csv')
colnames(rejci)[1] <- "Rejec"

library(ROracle)

#genetski trendi pri rjavi pasmi na kodo lastnosti 164 - mleko kg

drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")


con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

poizvedba <- paste("SELECT DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL ID_ZIVALI,
  ziv.ZIV_ID_SEQ,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL,
ziv.cre_SIFRA_CREDA
FROM zivali ziv, janao.genotipizirane_zivali gen
WHERE gen.ziv_id_seq = ziv.ziv_id_seq and gen.gen_chip='IDBv03'")
genTabela <- fetch(dbSendQuery(con,poizvedba))
#seznamA
#vseRjavePDF <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017_PDF.csv",header=T) 
vseRjavePDF1 <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_909_15032017_PDF_zeIzloceneRejeZoran.csv",header=T) 
#nakdadno odbrane, da se pokoristijo sredstva -11.5.2017
vseRjavePDF2 <- read.csv("~/Documents/F4F/OdbiraZivali/DodatnaOdbira_11052017/RjaveKrave_145_Dodatne_11052017_PDF.csv",header=T) 
vseRjavePDF <- rbind(vseRjavePDF1, vseRjavePDF2)
nrow(vseRjavePDF)

#Preveri spremembo črede
genCrede <- genTabela[,c("ID_ZIVALI", "CRE_SIFRA_CREDA")]
beforeCrede <- vseRjavePDF[,c("ID_ZIVALI", "CRE_SIFRA_CREDA")]
colnames(beforeCrede)[2] <- "CREDA_PREJ"
crede <- merge(genCrede, beforeCrede, by="ID_ZIVALI", all.x=T)
crede <- unique(crede)
crede[which(crede$CRE_SIFRA_CREDA != crede$CREDA_PREJ),]

######
#dobi SNPe za lastnosti
#######
library(ggplot2)
setwd("/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/")
tsPed <- read.table("MGTraits_11072017.ped")
tsMap <- read.table("./MGTraits_11072017.map")
ts <- tsPed[,-c(3,4,5,6)]

conv <- read.csv("./Trait_conversion_map_IDBV3_3.csv")

names <- conv[,c(1,3)]
colnames(tsMap)[2] <- "Name"


keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
namesSNP <- keeping.order(tsMap, merge, y=names, by = "Name")
namesSNP$TraitNo <- 1:40
traitNo <- rep(1:40, each=2)
new <- data.frame(t(data.frame(TraitNo = as.factor(c(0,0, traitNo)))))
colnames(new) <- colnames(ts)
ts <- rbind(new, ts) #tukaj oštevilči SNPe

Tts <- as.data.frame(t(ts))
#za vsako lastnost naredi tabelo

snpTable <- data.frame(FID=tsPed$V1, ID=tsPed$V2)
for (tNo in c(10, 11, 13, 14, 15, 32)) {
  snpTable_new <- as.data.frame(t(Tts[Tts$TraitNo %in% c(0,tNo),]))
  colnames(snpTable_new) <- c("FID", "ID", "A1", "A2")
  snpTable_new <- snpTable_new[rownames(snpTable_new) != "TraitNo",]
  snpTable_new$Genotype <- paste0(snpTable_new$A1, snpTable_new$A2)
  colnames(snpTable_new)[3:5] <- c(paste0("A1_", tNo), paste0("A2_", tNo), paste0("Genotype_",tNo) )
  snpTable <- merge(snpTable, snpTable_new, by=c("FID", "ID"))
  alleles <- c(snpTable$A1, snpTable$A2)
  #ggplot(as.data.frame(alleles)) +  geom_bar(aes(x = alleles, fill = as.factor(alleles)), position = "dodge", stat = "count") + ggtitle(namesSNP$Full.Trait.Name[namesSNP$TraitNo==tNo])
}

######################################
#sutvari tabele po čredah
#############################
colnames(crede)[1] <- "ID"
snpi <- merge(snpTable, crede, by="ID")

for (c in unique(snpi$CRE_SIFRA_CREDA)) {
  credaTable <- subset(snpi, snpi$CRE_SIFRA_CREDA==c)
  credaGen <- credaTable[,c(1,seq(5, ncol(credaTable), by=3))]
  colnames(credaGen) <- c("IDg", "Weaverg", "Arahnomelijag", "ABCG2g", "KappaCasein Bg", "KappaCasein Eg", "SMAg")
  credaGen$Weaver <- ifelse( credaGen$Weaverg == "BB", "Zdrava", ifelse( credaGen$Weaverg == "AB"| credaGen$Weaverg=="BA",  "Prenašalka", NA) )
  credaGen$Arahnomelija <- ifelse( credaGen$Arahnomelijag == "AA", "Zdrava", ifelse( credaGen$Arahnomelijag == "AB"| credaGen$Arahnomelijag =="BA",  "Prenašalka", NA) )
  credaGen$ABCG2 <- ifelse( credaGen$ABCG2g == "AA", "Dve kopiji alela A", ifelse( credaGen$ABCG2g == "AB"| credaGen$ABCG2g == "BA",  "Ena kopija alela A", NA) )
  credaGen$SMA <- ifelse( credaGen$SMAg == "BB", "Zdrava", ifelse( credaGen$SMAg == "AB"| credaGen$SMAg == "BA",  "Prenašalka", NA) )
}
