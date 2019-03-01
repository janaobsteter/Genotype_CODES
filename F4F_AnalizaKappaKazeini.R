reje1 <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_909_15032017.csv")
reje1 <- reje1[,c("CRE_SIFRA_CREDA", "ID_ZIVALI")]
reje3 <- read.csv("~/Documents/F4F/OdbiraZivali/seznamB_15032017.csv")
reje3 <- reje3[,c("CRE_SIFRA_CREDA", "ID_NadomestnaZival")]
colnames(reje3) <- c("CRE_SIFRA_CREDA", "ID_ZIVALI")
reje2 <- read.csv("~/Documents/F4F/OdbiraZivali/DodatnaOdbira_11052017/RjaveKrave_145_Dodatne_11052017_PDF.csv")
reje2 <- reje2[,c("CRE_SIFRA_CREDA", "ID_ZIVALI")]
reje4 <- read.csv("~/Documents/F4F/AnzelakDodatneZivali_12.csv", header=FALSE)
reje4$CRE_SIFRA_CREDA <- 143
colnames(reje4)[1] <- "ID_ZIVALI"
reje <- rbind(reje1, reje2)
reje <- rbind(reje, reje3)
reje <- rbind(reje, reje4)

kapaCSN <- read.csv('/home/jana/Documents/F4F/MlecniProteini/KappaCaseinGenotype_python.csv', header=TRUE)[,c(2,8)]
colnames(kapaCSN) <- c("ID", "KapaCSN")
kapaCSN$KapaCSN <- gsub("BB", "B/B", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("AA", "A/A", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("AB", "A/B", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("AC", "A/C", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("BC", "B/C", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("BE", "B/E", kapaCSN$KapaCSN)
colnames(kapaCSN) <- c("ID", "KapaCSN")

nrow(kapaCSN)
length(which(kapaCSN$ID %in% reje$ID_ZIVALI))
colnames(kapaCSN)[1] <- "ID_ZIVALI"
kapaCSN <- merge(kapaCSN, reje, by="ID_ZIVALI", all.x=TRUE)

rejci <- read.csv("~/Documents/F4F/Rezultati_MCPKlasiak///Rejci_SIfre.csv", header=TRUE)
colnames(rejci) <- c("Ime", "CRE_SIFRA_CREDA", "Kontrola")
rejci <- rejci[,1:2]
length(unique(kapaCSN$CRE_SIFRA_CREDA))

kapaCSN <- unique(merge(kapaCSN, rejci, by="CRE_SIFRA_CREDA", all.x=TRUE))

library(plyr)
A <- t(table(kapaCSN$KapaCSN, kapaCSN$Ime) )
write.table(A, "~/Documents/F4F/MlecniProteini/SteviloKazeiniPoCredah.csv", quote=FALSE, sep=";")
Ap <- round((A / rowSums(A) * 100), 0)
write.table(Ap, "~/Documents/F4F/MlecniProteini/ProcentiKazeiniPoCredah.csv", quote=FALSE, sep=";")  
