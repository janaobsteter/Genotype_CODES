#pedigree vseh rjavih živali
rjPed <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_pedigree10022017.csv")
rjPed$DATUM <- as.Date(rjPed$DAT_ROJSTVO, format="%d.%m.%Y")
rjPed4 <- rjPed[,c(1,2,3,9)] # obdrži le id, dam, sire, datum rojstva
rjPed4[,1:3] <- lapply(rjPed4[,1:3], as.numeric)
rjPed3 <- rjPed4[,1:3]

library(pedigree)
library(MCMCglmm)
library(pedantics)
library(MasterBayes)


rjKrave10 <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav.csv")
rjKrave10 <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKRave_MC_15032017.csv")
#Tukaj maš nekaj podvojenih z različnimi PV
rjKrave10 <- rjKrave10[!duplicated(rjKrave10$ZIV_ID_SEQ, fromFirst=T),]
#rjKrave10 <- read.csv("~/Documents/F4F/OdbiraZivali/VseRjKrave_TudiManjKot10.csv")
#rjKrave <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVse.csv") #to so vse krave, ki pridejo v poštev: oddajajo mlekarni, žive, v laktaciji
length(unique(rjKrave10$ZIV_ID_SEQ))
rjKrave10[which(rjKrave10$ZIV_ID_SEQ %in% rjKrave10$ZIV_ID_SEQ[duplicated(rjKrave10$ZIV_ID_SEQ)]),]
length(intersect(rjPed$ZIV_ID_SEQ, rjKrave10$ZIV_ID_SEQ))
length(which(rjKrave10$ZIV_ID_SEQ %in% rjPed$ZIV_ID_SEQ))

"
#poglej, kje se po SSI nahajajo krave v najvačjih čredah
Nad20Crede <- VseCredeDATA$CREDA[which(VseCredeDATA$RJKRAVE > 20)]
rjKrave20 <- rjKrave10[which(rjKrave10$CRE_SIFRA_CREDA %in% Nad20Crede),]

VseKR <- qplot(rjKrave10$VREDNOST_12_PV, geom='histogram', bins=100)
Nad20Kr <- qplot(rjKrave20$VREDNOST_12_PV, geom='histogram', bins=100)
multiplot(VseKR, Nad20Kr)
"
#prune the pedigree only for cows you wish to calculate relatedness for
CowsOrd <- MasterBayes::orderPed(rjPed3, time_born = rjPed4$DATUM)
prunedCows <- prunePed(CowsOrd, keep=rjKrave10$ZIV_ID_SEQ, make.base=F)
#prunedCowsBase <- prunePed(rjPed3, keep=rjKrave$ZIV_ID_SEQ, make.base=T)
colnames(prunedCows) <- c("id", "dam", "sire")
length(intersect(prunedCows$id, rjKrave10$ZIV_ID_SEQ))
prunedCows[,1:3] <- lapply(prunedCows[,1:3], as.numeric)


#fix the pedigree
prunedCowsFixed <- pedantics::fixPedigree(Ped=prunedCows)
#prunedCowsFixed[which(!(prunedCowsFixed$sire %in% prunedCowsFixed$id)),]
#length(which(!(prunedCowsFixed$dam %in% prunedCowsFixed$dam)))
#length(which(!(prunedCowsFixed$sire %in% prunedCowsFixed$sire)))

#calculate relatedness
#prunedCowsFixed <- MasterBayes::orderPed(prunedCowsFixed)
#set pedigree object
library(pedigreemm)
pedRjCows <- pedigreemm::pedigree(sire=as.numeric(as.character(prunedCowsFixed$sire)),
                     dam=as.numeric(as.character(prunedCowsFixed$dam)),
                     label=as.numeric(as.character((prunedCowsFixed$id))))

## --- Izracun koef. sorodstva ---

L <- relfactor(ped=pedRjCows)
object.size(L)
## 8251704 bytes
A <- crossprod(L)
object.size(A)
#205540192 bytes
rownames(A) <- colnames(A) <- prunedCowsFixed$id
## 7337

## --- Izbor ---

#podatki = rjave krave, za katere nas zanim sorodstvo - 3386 živali
podatki <- rjKrave10[which(rjKrave10$ZIV_ID_SEQ %in% rjPed$ZIV_ID_SEQ),]
podatki <- podatki[which(podatki$ZIV_ID_SEQ %in% prunedCowsFixed$id),]
## Samo cistopasemske --> že prej samo rjave


## Sestevek koef. sorodstva za dane zivali
A <- A[as.character(podatki$ZIV_ID_SEQ), as.character(podatki$ZIV_ID_SEQ)]
podatki$R <- rowSums(A)


#attach R to krave
R10 <- unique(podatki[,c('ZIV_ID_SEQ', 'R')])
rjKrave10R <- merge(rjKrave10, R10, by='ZIV_ID_SEQ', all.x=T)
#rjKrave10RnoNA <- rjKrave10R[-which(is.na(rjKrave10R$VREDNOST_12_PV)),]
#write.csv(rjKrave10R, "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusR.csv", row.names=F)
#write.csv(rjKrave10RnoNA, "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusRnoNA.csv", row.names=F)

## Porazdelitev vsote koef. sorod.
library(ggplot2)
library(doBy)
qplot(podatki$R, geom='histogram', bins=100)



################################
##################################
#samo vrži vn polsestre od vseh rjavih krav
require(Rmisc)
#rjPV <- qplot(rjKrave10R$VREDNOST_12_PV, bins=100)
#rjR <- qplot(rjKrave10R$R,  bins=100)
#multiplot(rjPV, rjR)

#izloči krave, ki nimajo znanega oceta
brezOceta <- rjKrave10R[which(is.na(rjKrave10R$OCE)),]
nrow(brezOceta)
rjKrave10R <- rjKrave10R[-which(is.na(rjKrave10R$OCE)),]
nrow(rjKrave10R)

#izloči patHS v isti čredi
rjKrave10R$together <- paste(rjKrave10R$CRE_SIFRA_CREDA, rjKrave10R$OCE, sep="")
rjKrave10R <- rjKrave10R[order(rjKrave10R$together),]
duplSis <- unique(rjKrave10R$together[(duplicated(rjKrave10R$together))])
length(duplSis)
length(intersect(duplSis, rjKrave10R$together))
halfSis <- rjKrave10R[which(rjKrave10R$together %in% duplSis),]
nrow(halfSis)

uniqueKrave <- rjKrave10R[which(!(rjKrave10R$together %in% duplSis)),]
nrow(uniqueKrave)

#tukaj imaš sedaj 536 vrstic v polsestrah, ostat pa ti more 233 živali
#zdj izberi sestre glede na manj sorodne
halfSis <- halfSis[order(halfSis$together, halfSis$R),]
leastR <- halfSis[!(duplicated(halfSis$together, fromFirst=T)),]
seznamBpatHS <- halfSis[(duplicated(halfSis$together, fromFirst=T)),]
seznamApatHS <- leastR[,c('ZIV_ID_SEQ','together')]
colnames(seznamApatHS)[1] <- 'IzbranaZival'
nadomPatHS <- merge(seznamBpatHS, seznamApatHS, by='together', all.x=T)
nadomPatHSPDF <- nadomPatHS[,c('CRE_SIFRA_CREDA', 'IzbranaZival', 'ID_ZIVALI', 'IME_ZIVAL', 'DAT_ROJSTVO'),]
nadomPatHSPDF$REL <- 'PatHS'

#to so vse rjaveKrave MC  brez polsester po očetu (aktivne krave v laktaciji)
vseRjavebrezHS <- rbind(uniqueKrave, leastR)
nrow(vseRjavebrezHS) #1247, 1313
length(intersect(nadomPatHS$IzbranaZival, vseRjavebrezHS$ZIV_ID_SEQ)) #preveri, če so Izbrane živali v seznamu B tudi v novem seznam odbranih živali (#233)


#poskusi drugače --> s tem v bistvu izločiš vse matere - hči in polsestre po mamini strani
matere <- vseRjavebrezHS[,c('ZIV_ID_SEQ','CRE_SIFRA_CREDA')]
colnames(matere) <- c("MATI", "CRE_SIFRA_CREDA")
hcere <- vseRjavebrezHS[,c('ZIV_ID_SEQ','CRE_SIFRA_CREDA', 'MATI')]
colnames(hcere) <- c("HCI", "CRE_SIFRA_CREDA", "MATI")
mh <- merge(matere, hcere, by=c("MATI", "CRE_SIFRA_CREDA")) #identificiraj matere in hčere v isti čredi (tukaj tudi trojke)
nrow(mh)
mh <- mh[order(mh$CRE_SIFRA_CREDA),]

#izloči matere
seznamBmatere <- vseRjavebrezHS[which(vseRjavebrezHS$ZIV_ID_SEQ %in% mhR$MATI),] #vse matere na seznam B
vseRjavebrezHS <- vseRjavebrezHS[-which(vseRjavebrezHS$ZIV_ID_SEQ %in% mh$MATI),]
nrow(vseRjavebrezHS)

Rhcer <- vseRjavebrezHS[,c('ZIV_ID_SEQ', 'R')] #tukaj so koef. sorodnostni hčera
colnames(Rhcer)[1] <- 'HCI'
mhR <- merge(mh, Rhcer, by='HCI', all.x=T)#dodaj R-je
nrow(mhR)
mhR$maticreda <- paste(mhR$CRE_SIFRA_CREDA,mhR$MATI,  sep="") #skupen čreda matere in id, da dobiš eno spremenljivko --> kje se ponovi
boljseHcere <- mhR[!(duplicated(mhR$maticreda, fromFirst=T)),] # hčere z višjim koeficientom sorodstva
slabseHcere <- mhR[which(!(mhR$HCI %in% boljseHcere$HCI)),]
nrow(slabseHcere)

#izloči slabše hčere
vseRjavebrezHS <- vseRjavebrezHS[-(which(vseRjavebrezHS$ZIV_ID_SEQ %in% slabseHcere$HCI)),]
nrow(vseRjavebrezHS)

slabseHcere <- slabseHcere[,c(1,3,5)]
colnames(slabseHcere)[1] <- c('NadomHS3')
boljseHcere <- boljseHcere[,c(1,5)]
colnames(boljseHcere)[1] <- 'IzbranaZival'
nadomMatHSTrojka <- merge(boljseHcere, slabseHcere, by='maticreda')
nadomMatHSTrojka <- nadomMatHSTrojka[order(nadomMatHSTrojka$IzbranaZival),]


#pripravi še seznam, kjer povežeš izločeno nater z izbrano živaljo
mhR1 <- mhR[,c(1,2)]
colnames(mhR1)[2] <- "ZIV_ID_SEQ"
nadomMater <- merge(mhR1, rjKrave10R, by='ZIV_ID_SEQ')
nadomMaterPDF <- nadomMater[,c('CRE_SIFRA_CREDA', 'HCI', 'ID_ZIVALI', 'IME_ZIVAL', 'DAT_ROJSTVO')] #tukaj so nadomestne živali matere - dodati moraš še podatke o materi: IME + DATUM ROJSTVA
nadomMaterPDF$REL <- 'MATI'
colnames(nadomMaterPDF)[2] <- 'IzbranaZival'

nadomMatHSTrojka <- merge(boljseHcere, slabseHcere, by='maticreda')
nadomMatHSTrojka <- nadomMatHSTrojka[order(nadomMatHSTrojka$IzbranaZival),]
colnames(nadomMatHSTrojka)[3] <- 'ZIV_ID_SEQ'
nadomMatHSTrojka <- merge(nadomMatHSTrojka, rjKrave10R, by=c('ZIV_ID_SEQ', 'CRE_SIFRA_CREDA'), all.x=T)
nadomMatHSTrojkaPDF <- nadomMatHSTrojka[,c('CRE_SIFRA_CREDA', 'IzbranaZival', 'ID_ZIVALI', 'IME_ZIVAL', 'DAT_ROJSTVO'),] #tukaj so zdj nadomestne živali živali pod ID_ZIVALI
nadomMatHSTrojkaPDF$REL <- 'MatHSTrojka'

#koliko teh je polsester po mami
vseRjavebrezHS$togetherM <- paste(vseRjavebrezHS$CRE_SIFRA_CREDA, vseRjavebrezHS$MATI, sep="")
length(unique(vseRjavebrezHS$togetherM))
duplMati <- vseRjavebrezHS$togetherM[duplicated(vseRjavebrezHS$togetherM)]
length(duplMati)
vseRjaveMdupl <- vseRjavebrezHS[which(vseRjavebrezHS$togetherM %in% duplMati),] # to so vse maternal polsestre, polovico jih obrži
nrow(vseRjaveMdupl)
vseRjaveMdupl <- vseRjaveMdupl[order(vseRjaveMdupl$togetherM, vseRjaveMdupl$R),]
seznamBmatHS <- vseRjaveMdupl[duplicated(vseRjaveMdupl$togetherM, fromLast=T),] #157 maternal HS

seznamAmatHS <- vseRjaveMdupl[duplicated(vseRjaveMdupl$togetherM, fromFirst=T),c('ZIV_ID_SEQ','togetherM')] #157 maternal HS
colnames(seznamAmatHS)[1] <- 'IzbranaZival'
nadomMatHS <- merge(seznamAmatHS, seznamBmatHS, by='togetherM')
nadomMatHSPDF <- nadomMatHS[,c('CRE_SIFRA_CREDA', 'IzbranaZival', 'ID_ZIVALI', 'IME_ZIVAL', 'DAT_ROJSTVO'),]
nadomMatHSPDF$REL <- 'MatHS'

#izloči iz seznama živali maternal HS in matere 
vseRjavebrezHS <- vseRjavebrezHS[-which(vseRjavebrezHS$ZIV_ID_SEQ %in% seznamBmatHS$ZIV_ID_SEQ),] #943
nrow(vseRjavebrezHS)  

#vseRjave <- vseRjavebrezHS[-(which(vseRjavebrezHS$ZIV_ID_SEQ %in% seznamBmatere$ZIV_ID_SEQ)),] # izloči matere, PREIMENUJE TABELO!  vseRjave
#vseRjave <- vseRjave[-(which(vseRjave$ZIV_ID_SEQ %in% seznamBhcere$ZIV_ID_SEQ)),] # tu gre za primere HS, kjer sta obe hčeri kot tudi mati v isti čredi
#vseRjave <- vseRjave[-(which(vseRjave$ZIV_ID_SEQ %in% seznamBmatHS$ZIV_ID_SEQ)),] # tu gre za primere HS, kjer sta obe hčeri kot tudi mati v isti čredi

#izloči še potencialne že genotipizirane
zenskeGen <- read.csv("~/Genotipi/ZenskeGen_15022017.csv")
length(which(vseRjavebrezHS$ZIV_ID_SEQ %in% zenskeGen$ZIV_ID_SEQ)) #takšnih živali je 6
vseRjave <- vseRjavebrezHS[-which(vseRjavebrezHS$ZIV_ID_SEQ %in% zenskeGen$ZIV_ID_SEQ),]
nrow(vseRjave)

#merge with animal names
#plusNames <- read.csv('/home/jana/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav.csv')
#plusNames <- plusNames[,c(1,2)]
#plusNames <- unique(plusNames)
#vseRjave <- merge(plusNames, vseRjave, by="ZIV_ID_SEQ", all.y=T)


#ustvari še B seznam
#matere
seznamB <- unique(rbind(nadomMaterPDF, nadomMatHSTrojkaPDF, nadomMatHSPDF, nadomPatHSPDF))
seznamB <- seznamB[order(seznamB$CRE_SIFRA_CREDA, seznamB$IzbranaZival),]
#dodaj ID izbrane živali
idZivali <- rjKrave10[,c('ZIV_ID_SEQ','ID_ZIVALI')]
colnames(idZivali) <- c('IzbranaZival', 'IzbranaZivalID')
seznamB <- unique(merge(idZivali, seznamB, by='IzbranaZival', all.y=T))
seznamB <- seznamB[,c(3,2,4:ncol(seznamB))]
colnames(seznamB)[2:3] <- c('ID_IzbranaZival', 'ID_NadomestnaZival')

#sedaj iz seznamaB izloči vrstice za tiste živali, ki so bile izločene kasneje, ampak so še vedno kot IZbrane na seznamuB
kasnejeIzl <- intersect(seznamB$ID_IzbranaZival, seznamB$ID_NadomestnaZival)
seznamB1 <- seznamB[-which(seznamB$ID_IzbranaZival %in% kasnejeIzl),]
length(unique(seznamB1$ID_IzbranaZival)) #373 -- more se ujemat s spodnjo vrstico
length(intersect(seznamB1$ID_IzbranaZival, vseRjave$ID_ZIVALI)) #373
seznamB <- seznamB1 #498 vrstic



#ZORAN KRamer izločil te črede
izlCrede <- as.character(c(3989,30968,10135,31629,15727,2264,31656,3110,11610))#te reje je izločil Zoran Kramer
vseRjave <- vseRjave[-which(vseRjave$CRE_SIFRA_CREDA %in% izlCrede),]
nrow(vseRjave)
seznamB <- seznamB[-which(seznamB$CRE_SIFRA_CREDA %in% izlCrede),]



#pusti samo stolpce za pdf sezname - št črede, id živali, ime živali, datum rojstva
vseRjavePDF <- vseRjave[,c('CRE_SIFRA_CREDA', 'ID_ZIVALI', 'IME_ZIVAL', 'DAT_ROJSTVO')]
  
#write.table(vseRjave, "~/Documents/F4F/OdbiraZivali/RjaveKrave_928_15022017.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
write.csv(vseRjave, "~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
write.csv(vseRjave, "~/Documents/F4F/OdbiraZivali/RjaveKrave_909_15032017.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
#write.csv(vseRjavePDF, "~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017_PDF.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
write.csv(vseRjavePDF, "~/Documents/F4F/OdbiraZivali/RjaveKrave_909_15032017_PDF_zeIzloceneRejeZoran.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
#vseRjave <- read.table("~/Documents/F4F/OdbiraZivali/RjaveKrave_928_15022017.csv", header=T) #TO SO ODBRANE ŽIVALI v čredah nad 10 - štart!!!
vseRjave <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017_PDF.csv", header=T) #TO SO ODBRANE ŽIVALI v čredah nad 10 - štart!!!
vseRjave <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017.csv", header=T) #TO SO ODBRANE ŽIVALI v čredah nad 10 - štart!!!
vseRjave <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_909_15032017.csv", header=T) #TO SO ODBRANE ŽIVALI v čredah nad 10 - štart!!!
#tukaj pa zapiši tiste, kje štartaš z vsemi čredami,ne samo pod 10
#write.table(vseRjave, "~/Documents/F4F/OdbiraZivali/RjaveKraveVse_1288_20022017.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
#vseRjave <- read.table("~/Documents/F4F/OdbiraZivali/RjaveKraveVse_1288_20022017.csv", header=T) #TO SO ODBRANE ŽIVALI v vseh čredah!!!
#write.table(seznamB, "~/Documents/F4F/seznamBVse_20022017.csv", row.names=F, quote=F)
write.csv(seznamB, "~/Documents/F4F/OdbiraZivali/seznamB_22022017.csv", row.names=F, quote=F)
write.csv(seznamB, "~/Documents/F4F/OdbiraZivali/seznamB_15032017.csv", row.names=F, quote=F)
seznamB <- read.table("~/Documents/F4F/OdbiraZivali/seznamB_22022017.csv", header=T)

reje <- as.data.frame(unique(vseRjave$CRE_SIFRA_CREDA))
write.csv(reje, "~/Documents/F4F/RejeSifre_15032017.csv")
################################################################################################################################################
################################################################################################################################################

#preveri, če je še vedno kaj polsester ali mati-hči parov
sum(duplicated(vseRjave$together)) #očetovske polsestre
sum(duplicated(vseRjave$togetherM)) #maternalne polsestre
a <- vseRjave[,c('ZIV_ID_SEQ', 'CRE_SIFRA_CREDA')]
colnames(a) <- c('MATI', 'CREDA')
b <- vseRjave[,c('ZIV_ID_SEQ', 'MATI', 'CRE_SIFRA_CREDA')]
colnames(b) <- c('HCI', 'MATI', 'CREDA')
c <- merge(a, b, by=c('MATI', 'CREDA')) #MORA BITI PRAZEN
nrow(c)


################################################################################################################################################
################################################################################################################################################

avgSSIRj <- read.csv("~/Documents/F4F/OdbiraZivali/VseCredeMC_avgSSI_RJAVE.csv")
avgSSIRj <- aggregate(vseRjave$VREDNOST_12_PV, by=c(list(vseRjave$CRE_SIFRA_CREDA)), FUN=function (x) (mean(x,na.rm=T)))
avgSSIRj1 <- aggregate(vseRjave$VREDNOST_12_PV, by=c(list(vseRjave$CRE_SIFRA_CREDA)), FUN=function (x) (sd(x, na.rm=T)))
avgSSIRj <- merge(avgSSIRj, avgSSIRj1, by='Group.1')
colnames(avgSSIRj) <- c("CRE_SIFRA_CREDA", "meanSSI", "sdSSI")

avgR <- aggregate(vseRjave$R, by=c(list(vseRjave$CRE_SIFRA_CREDA)), FUN=function (x) (mean(x,na.rm=T)))
sdR <- aggregate(vseRjave$R, by=c(list(vseRjave$CRE_SIFRA_CREDA)), FUN=function (x) (sd(x,na.rm=T)))
avgR <- merge(avgR, sdR,by='Group.1' )
colnames(avgR) <- c("CRE_SIFRA_CREDA", "meanR", "sdR")
avgData <- merge(avgSSIRj, avgR, by='CRE_SIFRA_CREDA')

vseRjave <- merge(vseRjave, avgData, by='CRE_SIFRA_CREDA', all.x=T)
vseRjave$deltaSSI <- (abs(vseRjave$meanSSI - vseRjave$VREDNOST_12_PV))/vseRjave$sdSSI #the larger the better
vseRjave$deltaR <- -((vseRjave$R - vseRjave$meanR)/vseRjave$sdR) #the smaller the better
vseRjave$Index <- 0.5*vseRjave$deltaSSI + 0.5*vseRjave$deltaR

#ODBERI ŽIVALI!!!!

#število po čredi
#stPocredah10 <- as.data.frame(table(vseRjave10$CRE_SIFRA_CREDA)) #tabela število krav po čredah
stPocredah <- as.data.frame(table(vseRjave$CRE_SIFRA_CREDA)) #tabela število krav po čredah
#write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/CredeInSteviloZivali.csv", quote=F, row.names=F)
colnames(stPocredah) <- c('CRE_SIFRA_CREDA', 'ST')

#zivaliData <- zivaliData[order(zivaliData$Index),] # tu so podvojene zivali
#zivaliData <- zivaliData[!(duplicated(zivaliData$ZIV_ID_SEQ, fromFirst=T)),] #izberi tiste z višjim indeksom (od podvojenih vrstic)
#length(intersect(zivaliData$ZIV_ID_SEQ, vseRjave$ZIV_ID_SEQ))
#selIndex <- zivaliData[,c('ZIV_ID_SEQ', 'Index')] #dodaj indeks rjavim kravam, kjer možno
#vseRjave <- merge(vseRjave, selIndex, by='ZIV_ID_SEQ', all.x=T)
#qplot(vseRjave$Index, geom='histogram', bins=500)

#TO SO PODATKI IZ IZBRANIH 821 ŽIVALI!!! in 78 ČRED!
#vseRjaveIndex <- vseRjave[which(!(is.na(vseRjave$Index))),]

#avgData
rejciIndex <- aggregate(vseRjave$Index, by=list(vseRjave$CRE_SIFRA_CREDA), function(x) mean(x, na.rm=TRUE))
colnames(rejciIndex) <- c("CRE_SIFRA_CREDA", "Index.mean")
avgData <- merge(avgData, rejciIndex, by='CRE_SIFRA_CREDA')
stPocredah <- merge(stPocredah, avgData, by="CRE_SIFRA_CREDA")

nrow(stPocredah[which(stPocredah$ST >= 10),])
sum(stPocredah$ST[which(stPocredah$ST >= 8)])
sum(stPocredah$ST[which(stPocredah$ST >= 10)])
qplot(stPocredah$Freq, geom='histogram', bins=100)

write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/SteviloPoCredah_Koncen.csv", quote=F, row.names = F)
write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/SteviloPoCredah_15032017.csv", quote=F, row.names = F)
#po indeksu in število po čredi
#izberi samo črede z več kot 10 rjavimi kravami
stPocredahVse <- stPocredah
stPocredah <- stPocredah[stPocredah$ST >=10,]
#stPocredah <- merge(stPocredah, rejciR, by="CRE_SIFRA_CREDA", all.x=T)

stPocredah <- stPocredah[order(stPocredah$Index.mean, stPocredah$ST, decreasing=T),]


#credeVecKot10 <- stPocredah[which(stPocredah$Freq >= 10),] #krave v čredah z več kot 10 kravami
#credeVecKot10 <- credeVecKot10[order(-credeVecKot10$Freq),] #order po številu živali

####################################################################################################
####################################################################################################
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


###################################################
###################################################
#izbrane samo po velikost
stPocredah <- stPocredah[order(stPocredah$ST, decreasing=T),]
izbraneCrede <- izberiCrede(280)
IzbraneVelikost <- vseRjave[which(vseRjave$CRE_SIFRA_CREDA %in% izbraneCrede),]
mean(IzbraneVelikost$R)
mean(IzbraneVelikost$Index, na.rm=T)
mean(IzbraneVelikost$VREDNOST_12_PV, na.rm=T)
sd(IzbraneVelikost$VREDNOST_12_PV, na.rm=T)
std.error(IzbraneVelikost$VREDNOST_12_PV, na.rm=T)


#PLOT po živalih
RvseRJ <- qplot(vseRjave$R, geom='histogram', bins=100)
RizRJ <- qplot(IzbraneVelikost$R, geom='histogram', bins=100)
multiplot(RvseRJ, RizRJ)

#PLOT po credah
VCrede <- avgData[which(avgData$CRE_SIFRA_CREDA %in% izbraneCrede),]
IzCrede <- qplot(VCrede$meanR, geom='histogram', bins=50) + xlim(c(45, 67))
vseCrede <- qplot(avgData$meanR, geom='histogram', bins=100) + xlim(c(45,67))
multiplot(IzCrede, vseCrede)


#############################################################################################
#izbrane po R-u in velikost
#############################################################################################
stPocredah <- stPocredah[order(stPocredah$meanR, -stPocredah$ST),]

izbraneCrede <- izberiCrede(280)
IzbraneR <- vseRjave[which(vseRjave$CRE_SIFRA_CREDA %in% izbraneCrede),]
mean(IzbraneR$R)
mean(IzbraneR$Index, na.rm=T)
mean(IzbraneR$VREDNOST_12_PV, na.rm=T)
sd(IzbraneR$VREDNOST_12_PV, na.rm=T)
std.error(IzbraneR$VREDNOST_12_PV, na.rm=T)


#PLOT po živalih
RvseRJ <- qplot(vseRjave$R, geom='histogram', bins=100)
RizRJ <- qplot(IzbraneR$R, geom='histogram', bins=100)
multiplot(RvseRJ, RizRJ)

#PLOT po credah
RCrede <- avgData[which(avgData$CRE_SIFRA_CREDA %in% izbraneCrede),]
IzCrede <- qplot(RCrede$meanR, geom='histogram', bins=100) + xlim(c(45, 67))
vseCrede <- qplot(avgData$meanR, geom='histogram', bins=100) + xlim(c(45,67))
multiplot(IzCrede, vseCrede)


#############################################################################################
#izbrane po Index-u in velikost
#############################################################################################
stPocredah <- stPocredah[order(-stPocredah$Index.mean, -stPocredah$ST),]

izbraneCrede <- izberiCrede(280)
IzbraneI <- vseRjave[which(vseRjave$CRE_SIFRA_CREDA %in% izbraneCrede),]
mean(IzbraneI$R)
mean(IzbraneI$Index, na.rm=T)
mean(IzbraneI$VREDNOST_12_PV, na.rm=T)
sd(IzbraneI$VREDNOST_12_PV, na.rm=T)
std.error(IzbraneI$VREDNOST_12_PV, na.rm=T)


#PLOT po živalih
IvseRJ <- qplot(vseRjave$R, geom='histogram', bins=100)
IizRJ <- qplot(IzbraneI$R, geom='histogram', bins=100)
multiplot(IvseRJ, IizRJ)

#PLOT po credah
ICrede <- avgData[which(avgData$CRE_SIFRA_CREDA %in% izbraneCrede),]
IzCrede <- qplot(ICrede$meanR, geom='histogram', bins=100) + xlim(c(45, 67))
vseCrede <- qplot(avgData$meanR, geom='histogram', bins=100) + xlim(c(45,67))
multiplot(IzCrede, vseCrede)


#############################################################################################
#############################################################################################
#pokritost po očetih
poOce <- as.data.frame(table(vseRjave$OCE))
poOce <- as.data.frame(table(IzbraneR$OCE))
poOce <- poOce[order(-poOce$Freq),]
qplot(poOce$Freq, geom='histogram', bins=100)

"""
vseRjaveICrede <- vseRjave[which(vseRjave$CRE_SIFRA_CREDA %in% izbraneCrede),]
mean(vseRjaveICrede$R)
mean(stPocredah$R.mean)
mean(stPocredahVse$R.mean)
ICrede <- stPocredah[which(stPocredah$CRE_SIFRA_CREDA %in% izbraneCrede),]
RCrede <- rejciR[which(rejciR$CRE_SIFRA_CREDA %in% izbraneCrede),]
IzCrede <- qplot(RCrede$R.mean, geom='histogram', bins=50) + xlim(c(45, 67))
vseCrede <- qplot(stPocredahVse$R.mean, geom='histogram', bins=100) + xlim(c(45,67))
multiplot(IzCrede, vseCrede)

PVCrede <- avgSSIRj[which(avgSSIRj$CREDA %in% izbraneCrede),]
IzCredePV <- qplot(PVCrede$AVGSSI, geom='histogram', bins=50) + xlim(c(80, 140))
vseCredePV <- qplot(CredeSSI$AVGSSI, geom='histogram', bins=50)+ xlim(c(80, 140))

multiplot(IzCredePV, vseCredePV)



#PREVERI R-je in PV-je izbranih živali
#vse reje mlekarne Celeia
#samoKrave, ne telice
podatkiKrave <- podatki[which(podatki$ZIV_ID_SEQ %in% vseRjave$ZIV_ID_SEQ),]
rejciR <- summaryBy(R ~ CRE_SIFRA_CREDA, FUN=descStat, data=podatkiKrave)
qplot(rejciR$R.mean, geom='histogram', bins=100)

#reje mlekarne Celeia z več kot 10 rjavimi kravami
rejci10R <- summaryBy(R ~ CRE_SIFRA_CREDA, FUN=descStat, data=podatki10)
qplot(rejci10R$R.mean, geom='histogram', bins=100)

#rejciUnder120 <- rejciR[which(rejciR$R.mean < 120),]
#rejci10Under120 <- rejci10R[which(rejci10R$CRE_SIFRA_CREDA %in% rejciUnder120$CRE_SIFRA_CREDA),]

avgSSIRj <- read.csv("~/Documents/F4F/OdbiraZivali/VseCredeMC_avgSSI_RJAVE.csv")
#merge R and SSI
colnames(avgSSIRj)[1] <- "CRE_SIFRA_CREDA"
credeData <- merge(rejciR, avgSSIRj, by="CRE_SIFRA_CREDA")



###TO STORI; KO ŽE IMAŠ KONČNI SEZNAM vseRjave!!!!!!!!!!!!!!!!!
#zdaj združi še z vsemi kravami
zivaliData <- merge(rjKrave10R, credeData, by="CRE_SIFRA_CREDA", all.x=T)
zivaliData$deltaSSI <- (abs(zivaliData$AVGSSI - zivaliData$VREDNOST_12_PV))/zivaliData$SSISD #the larger the better
zivaliData$deltaR <- -((zivaliData$R - zivaliData$R.mean)/zivaliData$R.sd) #the smaller the better
zivaliData$Index <- 0.5*zivaliData$deltaSSI + 0.5*zivaliData$deltaR

zivaliData$CRE_SIFRA_CREDA <- as.factor(zivaliData$CRE_SIFRA_CREDA)
qplot(zivaliData$Index, geom='histogram', fill=zivaliData$CRE_SIFRA_CREDA, bins=100)
#zivaliData <- zivaliData[-(which(is.na(zivaliData$Index))),]
nrow(zivaliData[which(zivaliData$Index > 0.5),])


write.csv(zivaliData, "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusRnoNA.csv", row.names=F)



zivaliData <- read.csv( "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusRnoNA.csv")




#pokritost po očetih
poOce <- as.data.frame(table(vseRjave$OCE))
poOce <- poOce[order(-poOce$Freq),]
qplot(poOce$Freq, geom='histogram', bins=100)




"""


