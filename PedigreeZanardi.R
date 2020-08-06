
#ZIV_ID, OCE_ID, MATI_ID, DATUMROJ (YYYYMMDD), SPOL (M/F)
#PAZI; da je oče pred mamo!!!

pedigreCB <- read.csv("~/Genotipi/Genotipi_DATA/Crnobele_pedigre.csv", header=T)
pedigreCB <- pedigreCB[,c(4,5,6,7,8)]
pedigreCB$DAT_ROJSTVO <- as.Date(pedigreCB$DAT_ROJSTVO, format="%d.%m.%Y", origin="01-01-1920")
pedigreCB$DAT_ROJSTVO <- format(pedigreCB$DAT_ROJSTVO, format="%Y%m%d")
pedigreCB$SIF_SPOL[which(pedigreCB$SIF_SPOL==1)] <- "M"
pedigreCB$SIF_SPOL[which(pedigreCB$SIF_SPOL==2)] <- "F"
pedigreCB <- pedigreCB[-c(which(pedigreCB$SIF_SPOL==3)),]
pedigreCB$MATIST[which(pedigreCB$MATIST=="")] <- "UUUUUUUU"
pedigreCB$OCEST[which(pedigreCB$OCEST=="")] <- "UUUUUUUU"
write.table(pedigreCB, "~/Genotipi/Genotipi_DATA/CBPedigre_ZANARDI", sep=";", row.names=F, col.names = F, quote=F, na="UUUUUUUUUUUUUU")



#pedigreRJ <- read.csv("~/Documents/Rjave_pedigre.csv", header=T)
pedigreRJ <- read.csv("~/Genotipi/Rjava_pedigree.csv", header=T)
#pedigreRJ <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_pedigree11012017.csv", header=T)
pedigreRJ <- pedigreRJ[,c(4,6,5,7,8)]
pedigreRJ$DAT_ROJSTVO <- as.Date(pedigreRJ$DAT_ROJSTVO, format="%d-%m-%y", origin="01-01-1920")
pedigreRJ$DAT_ROJSTVO <- format(pedigreRJ$DAT_ROJSTVO, format="%Y%m%d")

########################################
#fix pedigre

missingFather <- data.frame(ZIVAL=unique(pedigreRJ$OCEST[!(pedigreRJ$OCEST %in% pedigreRJ$ZIVAL)]),
                            OCEST=NA, MATIST=NA, DAT_ROJSTVO=NA, SIF_SPOL=1)
missingMother <- data.frame(ZIVAL=unique(pedigreRJ$MATIST[!(pedigreRJ$MATIST %in% pedigreRJ$ZIVAL)]),
                            OCEST=NA, MATIST=NA, DAT_ROJSTVO=NA, SIF_SPOL=2)
nrow(missingFather)
nrow(missingMother)

pedigreRJ <- rbind(missingFather, missingMother, pedigreRJ)

########################################

pedigreRJ$SIF_SPOL[which(pedigreRJ$SIF_SPOL==1)] <- "M"
pedigreRJ$SIF_SPOL[which(pedigreRJ$SIF_SPOL==2)] <- "F"
pedigreRJ <- pedigreRJ[-c(which(pedigreRJ$SIF_SPOL==3)),]
pedigreRJ$MATIST[which(pedigreRJ$MATIST=="")] <- "UUUUUUUU" #tudi če vrže vn opozorilo, je to OK, ker nardi NA in potem te NA zapišeš v UUUUUUUUUUUUUUU
pedigreRJ$OCEST[which(pedigreRJ$OCEST=="")] <- "UUUUUUUU"
pedigreRJ <- pedigreRJ[order(pedigreRJ$DAT_ROJSTVO),]

write.table(pedigreRJ, "~/Genotipi/Genotipi_DATA/RJPedigre_ZANARDI", sep=";", row.names=F, col.names = F, quote=F, na="UUUUUUUUUUUUUU")



#######################################
#če nočeš celega pedigreja, ampak samo za genotipizirane živali
#####################################
gen_zivali <- read.csv("~/Genotipi/Genotipi_DATA/Gen_zivali11012017.csv")
colnames(gen_zivali) <- "ID"
pedigreRJ <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_pedigree11012017.csv", header=T)
#pridobi pedigre le za gen zivali
#popravi - fixPedigree 
length(intersect(pedigreRJ$ZIVAL, gen_zivali$ID))
pedigreRJ_F <- pedigreRJ[which(pedigreRJ$ZIVAL %in% gen_zivali$ID),]
library(pedigree)
library(pedantics)
#pedigree <- genotyped_ped [,c(1,2,3)]
#zirhtaj pedigre uzgoraj (datum, spol ...)
pedigreRJ_F <- pedigreRJ_F[,c(1:3)]
pedigreRJ_Fix <- fixPedigree (Ped=pedigreRJ_F) 
colnames(pedigreRJ_Fix) <- c("ZIVAL", "dam", "sire")
pedigre_ADD <- pedigreRJ[,c(1,4,5)]
pedigre_Z <- merge(pedigreRJ_Fix, pedigre_ADD, by="ZIVAL")
pedigre_Z <- pedigre_Z[,c(1,3,2,4,5)] #OČE PRED MAMO!
#SOOOORTTTT!!!"!!!!!!!
pedigre_Z <- pedigre_Z[order(pedigre_Z$DAT_ROJSTVO),]
pedigre_Z$dam[which(pedigre_Z$dam=="")] <- "UUUUUUUU" #tudi če vrže vn opozorilo, je to OK, ker nardi NA in potem te NA zapišeš v UUUUUUUUUUUUUUU
pedigre_Z$sire[which(pedigre_Z$sire=="")] <- "UUUUUUUU"
write.table(pedigre_Z, "~/Genotipi/Genotipi_DATA/RJPedigre_ZANARDI", sep=";", row.names=F, col.names = F, quote=F, na="UUUUUUUUUUUUUU")

#za blupf90
pedigreZ <- read.csv("~/Genotipi/Genotipi_DATA/RJPedigre_ZANARDI", sep=";", header=FALSE)
colnames(pedigreZ) <- c("ZIVAL", "OCE", "MATI")
pedigreZ$OCE <- gsub("UUUUUUUUUUUUUU", 0, pedigreZ$OCE)
pedigreZ$MATI <- gsub("UUUUUUUUUUUUUU", 0, pedigreZ$MATI)
write.table(pedigreZ[,c(1,2,3)], "/home/jana/Genotipi/Genotipi_WORK/PedigreIMPUTED_blup.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
############################################

pedigreLS <- read.csv("~/Genotipi/Genotipi_DATA/Lisaste/Pedigree.csv", header=T)
pedigreLS <- pedigreLS[,c(4,5,6,7,8)]
pedigreLS$DAT_ROJSTVO <- as.Date(pedigreLS$DAT_ROJSTVO, format="%d.%m.%Y", origin="01-01-1920")
pedigreLS$DAT_ROJSTVO <- format(pedigreLS$DAT_ROJSTVO, format="%Y%m%d")
pedigreLS$SIF_SPOL[which(pedigreLS$SIF_SPOL==1)] <- "M"
pedigreLS$SIF_SPOL[which(pedigreLS$SIF_SPOL==2)] <- "F"
pedigreLS <- pedigreLS[-c(which(pedigreLS$SIF_SPOL==3)),]
#pedigreLS$MATIST[which(is.na(pedigreLS$MATI))] <- "UUUUUUUU"
#pedigreLS$OCEST[which(pedigreLS$OCE=="")] <- "UUUUUUUU"
write.table(pedigreLS, "~/Genotipi/Genotipi_DATA/Lisaste/LSPedigre_ZANARDI", sep=";", row.names=F, col.names = F, quote=F, na="UUUUUUUUUUUUUU")
