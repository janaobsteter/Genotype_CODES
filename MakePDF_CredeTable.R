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

#seznamA
#vseRjavePDF <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017_PDF.csv",header=T) 
vseRjavePDF <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_909_15032017_PDF_zeIzloceneRejeZoran.csv",header=T) 
#nakdadno odbrane, da se pokoristijo sredstva -11.5.2017
vseRjavePDF <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_145_Dodatne_11052017_PDF.csv",header=T) 
nrow(vseRjavePDF)

#IZLOČI ŠE NEMŠKE KRAVE!!!
drzave <-substr(vseRjavePDF$ID_ZIVALI, start=0, stop=2)
tuje <- which(drzave != "SI")
length(tuje)
vseRjavePDF <- vseRjavePDF[-tuje,]
nrow(vseRjavePDF)
length(unique(vseRjavePDF$CRE_SIFRA_CREDA))
which(substr(vseRjavePDF$ID_ZIVALI, start=0, stop=2)!="SI")

vseRjavePDF$DAT_ROJSTVO <- as.Date(vseRjavePDF$DAT_ROJSTVO, format='%d.%m.%Y')
vseRjavePDF$CRE_SIFRA_CREDA <- as.character(vseRjavePDF$CRE_SIFRA_CREDA)
nrow(vseRjavePDF)

#ŠE ENRKAT PREVERI; ČE JE KATERA ŽE GENOTIPIZRANA - azurno stanje
genRJ <- read.csv("~/Genotipi/Rjave_Gen_23032017.csv")
genRJ <- read.csv("~/Genotipi/Rjave_Gen_10052017.csv")
zeGen <- which(vseRjavePDF$ID_ZIVALI %in% genRJ$GEN.GEN_DRZ..GEN.GEN_ORIG_STEV)
if (length(zeGen) != 0) {
  vseRjavePDF <- vseRjavePDF[-zeGen,]
}
nrow(vseRjavePDF)

#seznamB
seznamB <- read.csv("~/Documents/F4F/OdbiraZivali/seznamB_62_Dodatne_11052017.csv", header=T)
seznamB$REL <- as.character(seznamB$REL)
seznamB$REL[seznamB$REL=="MatHSTrojka"] <- "MatiPS"
seznamB$REL[seznamB$REL=="MatHS"] <- "MatiPS"
seznamB$REL[seznamB$REL=="PatHS"] <- "OcePS"
seznamB$DAT_ROJSTVO <- as.Date(seznamB$DAT_ROJSTVO, format='%d.%m.%Y')
seznamB$Datum <- format(seznamB$DAT_ROJSTVO, format='%Y%m%d')
seznamB$DAT_ROJSTVO <- format(seznamB$DAT_ROJSTVO, format='%d.%m.%Y')
seznamB$Datum <- as.numeric(seznamB$Datum)

#IZLOČI ŠE NEMŠKE KRAVE!!!
nrow(seznamB)
drzave <-substr(seznamB$ID_NadomestnaZival, start=0, stop=2)
tuje <- which(drzave != "SI")
length(tuje)
if (length(tuje) != 0) {
  seznamB <- seznamB[-tuje,]
}
nrow(seznamB)
which(substr(seznamB$ID_NadomestnaZival, start=0, stop=2)!="SI")

#stPocredah10 <- as.data.frame(table(vseRjave10$CRE_SIFRA_CREDA)) #tabela število krav po čredah
stPocredah <- as.data.frame(table(vseRjavePDF$CRE_SIFRA_CREDA)) #tabela število krav po čredah
write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/StPoCredah.csv")
write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/StPoCredah_Dodatne.csv")
#write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/CredeInSteviloZivali.csv", quote=F, row.names=F)
colnames(stPocredah) <- c('CRE_SIFRA_CREDA', 'ST')

#izberi živali v čredah (glede na število živali za genotipizacijo)
#stPocredah <- read.csv("~/Documents/F4F/OdbiraZivali/SteviloPoCredah_15032017.csv")
#Črede več kot 10 živali
stPocredah1 <- stPocredah[stPocredah$ST>= 10,]
sum(stPocredah$ST)
sum(stPocredah1$ST)
#TUKAJ JE KRITERIJ ZA SORTIRANJE
stPocredah <- stPocredah[order(-stPocredah$ST),]

#CredeLJ in NM
credeLJNM <- c(86,2475,9690,10253)
zivaliLJNM <- vseRjavePDF[which(vseRjavePDF$CRE_SIFRA_CREDA %in% credeLJNM),]
nrow(zivaliLJNM)
#stPocredah <- stPocredah[-which(stPocredah$CRE_SIFRA_CREDA %in% credeLJNM),]



#15.3.2017 --> 558 rjavih kravih v čredah večjih kot 10 krav!¨
#poglej R in PV
qplot(stPocredah$meanR, geom='histogram', bins=100)
qplot(stPocredah$meanSSI, geom='histogram', bins=100)
qplot(stPocredah$Index.mean.x, geom='histogram', bins=100)


#NAKDADNO IZBRANE ČREDE - DA SE DOPOLNI OBLJUBLJENO ŠTEVILO ŽIVALI: 11.5.2017




#funkcija za odbiro določenega števila živali izmed 1125 živali -v čredah z več kot 10 živalmi (brez polsester)
#začne z največjimi čredami
izberiCrede <- function (st) {
  sum <- 0
  crede <- c()
  row <- 1
  while (sum < st) {
    sum <- sum + stPocredah$ST[row]
    crede <- c(crede, as.character(stPocredah$CRE_SIFRA_CREDA[row]))
    row <- row +1
  }
  print(sum)
  return(crede)
}

izbraneCrede <- izberiCrede(556)
length(izbraneCrede) #tukaj so izbrane kmetije glede na velikost 
#DODATNE ŽIVALI 11.5.2017
izbraneCrede <- c(2992, 3497, 4725,5825,9624,13087,8946)
length(intersect(c(86,2475,9690,10253), vseRjavePDF$CRE_SIFRA_CREDA)) #LJ in NM kmetije
izbraneCrede <- unique(c(izbraneCrede, c(86,2475,9690,10253))) #tukaj združi s štirimi NM in LJ kmetijami - 3NM in 1 LJ - te štiri MORAJO biti vključene ne glede na št krav(ampak 3 že tako ali tako po velikosti)
length(izbraneCrede)
length(intersect(vseRjavePDF$CRE_SIFRA_CREDA, izbraneCrede))

#Izbrane živali
length(which(vseRjavePDF$CRE_SIFRA_CREDA %in% izbraneCrede))
vseRjavePDFIzb <- vseRjavePDF[(which(vseRjavePDF$CRE_SIFRA_CREDA %in% izbraneCrede)),]
nrow(vseRjavePDFIzb)
seznamBIzb <- seznamB[which(seznamB$CRE_SIFRA_CREDA %in% izbraneCrede),]
nrow(seznamBIzb)


#SEZNAMI ZA KONTOLORJE
setwd("~/Documents/F4F/OdbiraZivali/DodatnaOdbira_11052017/")
## knitr loop
for (creda in izbraneCrede){ #unique(vseRjavePDF$CRE_SIFRA_CREDA)
  knit2pdf("~/Genotipi/Genotipi_CODES/MakePDF_CredeTable.Rnw", output=paste0('PDF_', creda, '.tex'))
}

#SEZNAMI ZA REJCE
setwd("~/Documents/F4F/OdbiraZivali/DodatnaOdbira_11052017/")
## knitr loop
for (creda in izbraneCrede){ #unique(vseRjavePDF$CRE_SIFRA_CREDA)
  knit2pdf("~/Genotipi/Genotipi_CODES/MakePDF_CredeTable_REJCI.Rnw", output=paste0('PDF_Rejci_', creda, '.tex'))
}


#sekvence za Andrejo
#seznamA
vseRJSeq <- rjKrave10[c(1,5)]
vseRjavePDFIzb <- merge(vseRJSeq, vseRjavePDFIzb, by="ID_ZIVALI", all.y=T)
vseRjaveIzb <- data.frame(ZIV_ID_SEQ=vseRjavePDFIzb[,2], SEZNAM="A")
vseRjaveIzb <- data.frame(ZIV_ID_SEQ=vseRjavePDFIzb[,c(2,3)], SEZNAM="A")
nrow(vseRjavePDFIzb)

#seznamB
vseRJSeq <- rjKrave10[c(1,5)]
colnames(seznamBIzb)[3] <- "ID_ZIVALI"
seznamBIzb <- merge(vseRJSeq, seznamBIzb, by="ID_ZIVALI", all.y=T)
seznamBIzb <- data.frame(ZIV_ID_SEQ=seznamBIzb[,2], SEZNAM="B")
seznamBIzb <- data.frame(ZIV_ID_SEQ=seznamBIzb[,c(2,3)], SEZNAM="B")
nrow(seznamBIzb)

celSeznam <- rbind(vseRjaveIzb, seznamBIzb)
colnames(celSeznam)
nrow(celSeznam)
write.csv(celSeznam, "~/Documents/F4F/OdbiraZivali/CelSeznamAplusB_15032017.csv", quote=F, row.names=F)
write.csv(celSeznam, "~/Documents/F4F/OdbiraZivali/DodatnaOdbira_11052017/CelSeznamAplusB_Dodaten_11052017.csv", quote=F, row.names=F)



##################
#problemi z imeni, čreda 7177
prob <- vseRjavePDF[which(vseRjavePDF$CRE_SIFRA_CREDA==7177),]

izbraniRejci <- as.data.frame(rejci[rejci$F4F_CRE_SIFRA_CREDA %in% izbraneCrede,])
colnames(izbraniRejci)[2] <- "CRE_SIFRA_CREDA"
izbraniRejci <- merge(izbraniRejci, stPocredah, by="CRE_SIFRA_CREDA")
colnames(izbraniRejci)[2] <- "SteviloZivali"
write.table(izbraniRejci[,c(2,3)], "~/Documents/F4F/OdbiraZivali/IzbraneReje_15032017.csv", quote=F, row.names=F, sep="\t")
write.table(izbraniRejci, "~/Documents/F4F/OdbiraZivali/IzbraneReje_inStZivali_15032017.csv", quote=F, row.names=F, sep="\t")
write.table(izbraniRejci, "~/Documents/F4F/OdbiraZivali/IzbraneReje_inStZivali_DodatnaIzbira_11052017.csv", quote=F, row.names=F, sep="\t")




