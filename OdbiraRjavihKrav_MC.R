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
#rjKrave10 <- read.csv("~/Documents/F4F/OdbiraZivali/VseRjKrave_TudiManjKot10.csv")
rjKrave <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVse.csv") #to so vse krave, ki pridejo v poštev: oddajajo mlekarni, žive, v laktaciji
length(intersect(rjPed$ZIV_ID_SEQ, rjKrave$ZIV_ID_SEQ))


#poglej, kje se po SSI nahajajo krave v najvačjih čredah
Nad20Crede <- VseCredeDATA$CREDA[which(VseCredeDATA$RJKRAVE > 20)]
rjKrave20 <- rjKrave10[which(rjKrave10$CRE_SIFRA_CREDA %in% Nad20Crede),]

VseKR <- qplot(rjKrave10$VREDNOST_12_PV, geom='histogram', bins=100)
Nad20Kr <- qplot(rjKrave20$VREDNOST_12_PV, geom='histogram', bins=100)
multiplot(VseKR, Nad20Kr)

#prune the pedigree only for cows you wish to calculate relatedness for
CowsOrd <- MasterBayes::orderPed(rjPed3, time_born = rjPed4$DATUM)
prunedCows <- prunePed(CowsOrd, keep=rjKrave$ZIV_ID_SEQ, make.base=F)
  #prunedCowsBase <- prunePed(rjPed3, keep=rjKrave$ZIV_ID_SEQ, make.base=T)
length(intersect(prunedCows$id, rjKrave$ZIV_ID_SEQ))
colnames(prunedCows) <- c("id", "dam", "sire")
prunedCows[,1:3] <- lapply(prunedCows[,1:3], as.numeric)


#fix the pedigree
prunedCowsFixed <- pedantics::fixPedigree(Ped=prunedCows)
#prunedCowsFixed[which(!(prunedCowsFixed$sire %in% prunedCowsFixed$id)),]
#length(which(!(prunedCowsFixed$dam %in% prunedCowsFixed$dam)))
#length(which(!(prunedCowsFixed$sire %in% prunedCowsFixed$sire)))

#calculate relatedness
#prunedCowsFixed <- MasterBayes::orderPed(prunedCowsFixed)
#set pedigree object
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
podatki <- rjKrave[which(rjKrave$ZIV_ID_SEQ %in% rjPed$ZIV_ID_SEQ),]
podatki <- podatki[which(podatki$ZIV_ID_SEQ %in% prunedCowsFixed$id),]
## Samo cistopasemske --> že prej samo rjave


## Sestevek koef. sorodstva za dane zivali
A <- A[as.character(podatki$ZIV_ID_SEQ), as.character(podatki$ZIV_ID_SEQ)]
podatki$R <- rowSums(A)


#attach R to krave
R10 <- podatki[,c(1,12)]
rjKrave10R <- merge(rjKrave10, R10, by='ZIV_ID_SEQ', all.x=T)
#rjKrave10RnoNA <- rjKrave10R[-which(is.na(rjKrave10R$VREDNOST_12_PV)),]
#write.csv(rjKrave10R, "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusR.csv", row.names=F)
#write.csv(rjKrave10RnoNA, "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusRnoNA.csv", row.names=F)

## Porazdelitev vsote koef. sorod.
library(ggplot2)
qplot(podatki$R, geom='histogram', bins=100)

"""
#vse reje mlekarne Celeia
#samoKrave, ne telice
podatkiKrave <- podatki[which(podatki$ZIV_ID_SEQ %in% rjKrave10R$ZIV_ID_SEQ),]
rejciR <- summaryBy(R ~ CRE_SIFRA_CREDA, FUN=descStat, data=podatkiKrave)
qplot(rejciR$R.mean, geom='histogram', bins=100)

#reje mlekarne Celeia z več kot 10 rjavimi kravami
rejci10R <- summaryBy(R ~ CRE_SIFRA_CREDA, FUN=descStat, data=podatki10)
qplot(rejci10R$R.mean, geom='histogram', bins=100)

rejciUnder120 <- rejciR[which(rejciR$R.mean < 120),]
rejci10Under120 <- rejci10R[which(rejci10R$CRE_SIFRA_CREDA %in% rejciUnder120$CRE_SIFRA_CREDA),]

avgSSIRj <- read.csv("~/Documents/F4F/OdbiraZivali/VseCredeMC_avgSSI_RJAVE.csv")
#merge R and SSI
colnames(avgSSIRj)[1] <- "CRE_SIFRA_CREDA"
credeData <- merge(rejciR, avgSSIRj, by="CRE_SIFRA_CREDA")


#zdaj združi še z vsemi kravami
zivaliData <- merge(rjKrave10R, credeData, by="CRE_SIFRA_CREDA", all.x=T)
zivaliData$deltaSSI <- (abs(zivaliData$AVGSSI - zivaliData$VREDNOST_12_PV))/zivaliData$SSISD #the larger the better
zivaliData$deltaR <- -((zivaliData$R - zivaliData$R.mean)/zivaliData$R.sd) #the smaller the better
zivaliData$Index <- 0.5*zivaliData$deltaSSI + 0.5*zivaliData$deltaR

zivaliData$CRE_SIFRA_CREDA <- as.factor(zivaliData$CRE_SIFRA_CREDA)
qplot(zivaliData$Index, geom='histogram', fill=zivaliData$CRE_SIFRA_CREDA, bins=100)
zivaliData <- zivaliData[-(which(is.na(zivaliData$Index))),]
nrow(zivaliData[which(zivaliData$Index > 0.5),])


#write.csv(zivaliData, "~/Documents/F4F/OdbiraZivali/RjaveKrave_Telitve_RejeVecKot10RjKrav_plusRnoNA.csv", row.names=F)
"""



################################
##################################
#samo vrži vn polsestre od vseh rjavih krav
require(Rmisc)
rjPV <- qplot(rjKrave10R$VREDNOST_12_PV, bins=100)
rjR <- qplot(rjKrave10R$R,  bins=100)
multiplot(rjPV, rjR)

#izloči krave, ki nimajo znanega oceta
rjKrave10R <- rjKrave10R[-which(is.na(rjKrave10R$OCE)),]
rjKrave10R$together <- paste(rjKrave10R$CRE_SIFRA_CREDA, rjKrave10R$OCE, sep="")
rjKrave10R <- rjKrave10R[order(rjKrave10R$together),]
duplSis <- rjKrave10R$together[which(duplicated(rjKrave10R$together))]
halfSis <- rjKrave10R[which(rjKrave10R$together %in% duplSis),]




uniqueKrave <- rjKrave10R[which(!(rjKrave10R$together %in% duplSis)),]
nrow(uniqueKrave)

#zdj izberi sestre glede na manj sorodne
halfSis <- halfSis[order(halfSis$together, halfSis$R),]
leastR <- halfSis[!(duplicated(halfSis$together, fromFirst=T)),]
seznamBpatHS <- halfSis[(duplicated(halfSis$together, fromFirst=T)),]

#to so vse rjaveKrave MC  brez polsester po očetu (aktivne krave v laktaciji)
vseRjavebrezHS <- rbind(uniqueKrave, leastR)
nrow(vseRjavebrezHS)


#poskusi drugače --> s tem v bistvu izločiš vse matere - hči in polsestre po mamini strani
matere <- vseRjavebrezHS[,c(1,3)]
colnames(matere) <- c("MATI", "CRE_SIFRA_CREDA")
hcere <- vseRjavebrezHS[,c(1,3,6)]
colnames(hcere) <- c("HCI", "CRE_SIFRA_CREDA", "MATI")
mh <- merge(matere, hcere, by=c("MATI", "CRE_SIFRA_CREDA")) #identificiraj matere in hčere v isti čredi (tukaj tudi trojke)
nrow(mh)
mh <- mh[order(mh$CRE_SIFRA_CREDA),]
Rhcer <- vseRjavebrezHS[,c('ZIV_ID_SEQ', 'R')] #tukaj so koef. sorodnostni hčera
colnames(Rhcer)[1] <- 'HCI'
mhR <- merge(mh, Rhcer, by='HCI')#dodaj R-je
nrow(mhR)
mhR$maticreda <- paste(mhR$CRE_SIFRA_CREDA,mhR$MATI,  sep="") #skupen čreda matere in id, da dobiš eno spremenljivko --> kje se ponovi
mhR <- mhR[order(mhR$maticreda, mhR$R),]
slabseHcere <- mhR[duplicated(mhR$maticreda, fromLast=T),] # hčere z višjim koeficientom sorodstva
nrow(slabseHcere)

#izloči matere in slabše hčere
seznamBmatere <- vseRjavebrezHS[which(vseRjavebrezHS$ZIV_ID_SEQ %in% mhR$MATI),] #vse matere na seznam B
seznamBhcere <- vseRjavebrezHS[which(vseRjavebrezHS$ZIV_ID_SEQ %in% slabseHcere$HCI),] #polsestre z višjim R v isti čredi (trojke) na seznam B



#koliko teh je polsester po mami
vseRjavebrezHS$togetherM <- paste(vseRjavebrezHS$CRE_SIFRA_CREDA, vseRjavebrezHS$MATI, sep="")
duplMati <- vseRjavebrezHS$togetherM[(duplicated(vseRjavebrezHS$togetherM))]
vseRjaveMdupl <- vseRjavebrezHS[which(vseRjavebrezHS$togetherM %in% duplMati),] # to so vse maternal polsestre, polovico jih obrži
vseRjaveMdupl <- vseRjaveMdupl[order(vseRjaveMdupl$togetherM, vseRjaveMdupl$R),]
seznamBmatHS <- vseRjaveMdupl[duplicated(vseRjaveMdupl$togetherM, fromLast=T),] #157 maternal HS

#izloči iz seznama živali maternal HS in matere 
vseRjave <- vseRjavebrezHS[-(which(vseRjavebrezHS$ZIV_ID_SEQ %in% seznamBmatere$ZIV_ID_SEQ)),] # izloči matere, PREIMENUJE TABELO!  vseRjave
vseRjave <- vseRjave[-(which(vseRjave$ZIV_ID_SEQ %in% seznamBhcere$ZIV_ID_SEQ)),] # tu gre za primere HS, kjer sta obe hčeri kot tudi mati v isti čredi
vseRjave <- vseRjave[-(which(vseRjave$ZIV_ID_SEQ %in% seznamBmatHS$ZIV_ID_SEQ)),] # tu gre za primere HS, kjer sta obe hčeri kot tudi mati v isti čredi

#izloči še potencialne že genotipizirane
zenskeGen <- read.csv("~/Genotipi/ZenskeGen_15022017.csv")
length(which(vseRjave$ZIV_ID_SEQ %in% zenskeGen$ZIV_ID_SEQ)) #takšnih živali je 6
vseRjave <- vseRjave[-which(vseRjave$ZIV_ID_SEQ %in% zenskeGen$ZIV_ID_SEQ),]

#write.table(vseRjave, "~/Documents/F4F/OdbiraZivali/RjaveKrave_928_15022017.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
vseRjave10 <- read.table("~/Documents/F4F/OdbiraZivali/RjaveKrave_928_15022017.csv", header=T) #TO SO ODBRANE ŽIVALI v čredah nad 10 - štart!!!
#tukaj pa zapiši tiste, kje štartaš z vsemi čredami,ne samo pod 10
#write.table(vseRjave, "~/Documents/F4F/OdbiraZivali/RjaveKraveVse_1288_20022017.csv", row.names=F, quote=F) #TO SO ODBRANE ŽIVALI!!!
vseRjave <- read.table("~/Documents/F4F/OdbiraZivali/RjaveKraveVse_1288_20022017.csv", header=T) #TO SO ODBRANE ŽIVALI v vseh čredah!!!


#pokritost po očetih
poOce <- as.data.frame(table(vseRjave$OCE))
poOce <- poOce[order(-poOce$Freq),]
qplot(poOce$Freq, geom='histogram', bins=100)

#ustvari še B seznam
#matere
seznamBmatere$rel <- "MATI"
seznamBhcere$rel <- "matSESTRA_trojka"
seznamBmatHS$rel <- "matSESTRA"
seznamBpatHS$togetherM <- paste(seznamBpatHS$CRE_SIFRA_CREDA, seznamBpatHS$MATI, sep="") 
seznamBpatHS$rel <- "patSESTRA"

seznamB <- rbind(seznamBmatere, seznamBhcere, seznamBmatHS, seznamBpatHS)
seznamB <- seznamB[order(seznamB$CRE_SIFRA_CREDA),]
#write.table(seznamB, "~/Documents/F4F/seznamBVse_20022017.csv", row.names=F, quote=F)


#preveri, če je še vedno kaj polsester ali mati-hči parov
sum(duplicated(vseRjave$together)) #očetovske polsestre
sum(duplicated(vseRjave$togetherM)) #maternalne polsestre
a <- vseRjave[,c('ZIV_ID_SEQ', 'CRE_SIFRA_CREDA')]
colnames(a) <- c('MATI', 'CREDA')
b <- vseRjave[,c('ZIV_ID_SEQ', 'MATI', 'CRE_SIFRA_CREDA')]
colnames(b) <- c('HCI', 'MATI', 'CREDA')
c <- merge(a, b, by=c('MATI', 'CREDA')) #MORA BITI PRAZEN


#parametri odbranih živali
selR <- qplot(vseRjave$R, geom='histogram', bins=100) #preveri R in PV teh krav
selPV <- qplot(vseRjave$VREDNOST_12_PV, geom='histogram', bins=100)
allPV <- qplot(podatki$VREDNOST_12_PV, geom='histogram', bins=100)
allR <- qplot(podatki$R, geom='histogram', bins=100)

multiplot(allPV, selPV)
multiplot(allR, selR)

zivaliData <- zivaliData[order(zivaliData$Index),] # tu so podvojene zivali
zivaliData <- zivaliData[!(duplicated(zivaliData$ZIV_ID_SEQ, fromFirst=T)),] #izberi tiste z višjim indeksom (od podvojenih vrstic)
selIndex <- zivaliData[,c('ZIV_ID_SEQ', 'Index')] #dodaj indeks rjavim kravam, kjer možno
vseRjaveI <- merge(vseRjave, selIndex, by='ZIV_ID_SEQ', all.x=T)
qplot(vseRjaveI$Index, geom='histogram', bins=500)
length(which(vseRjaveI$Index <= 0)) #159 jih ima indeks pod 0
length(which(vseRjaveI$Index >= 0))

#število po čredi
stPocredah10 <- as.data.frame(table(vseRjave10$CRE_SIFRA_CREDA)) #tabela število krav po čredah
stPocredah <- as.data.frame(table(vseRjave$CRE_SIFRA_CREDA)) #tabela število krav po čredah
#write.csv(stPocredah, "~/Documents/F4F/OdbiraZivali/CredeInSteviloZivali.csv", quote=F, row.names=F)
stPocredah <- stPocredah[order(-stPocredah$Freq),]
nrow(stPocredah10[which(stPocredah10$Freq >= 10),])
nrow(stPocredah[which(stPocredah$Freq >= 10),])
sum(stPocredah$Freq[which(stPocredah$Freq >= 8)])
sum(stPocredah$Freq[which(stPocredah$Freq >= 10)])
sum(stPocredah10$Freq[which(stPocredah10$Freq >= 10)])
qplot(stPocredah$Freq, geom='histogram', bins=100)


#po indeksu in število po čredi
vseRjaveI <- vseRjaveI[order(-vseRjaveI$Index),]
stPocredah <- as.data.frame(table(vseRjaveI$CRE_SIFRA_CREDA)) #tabela število krav po čredah
stPocredah <- stPocredah[order(-stPocredah$Freq),]
length(which(stPocredah$Freq >= 10))
credeVecKot10 <- stPocredah[which(stPocredah$Freq >= 10),] #krave v čredah z več kot 10 kravami
credeVecKot10 <- credeVecKot10[order(-credeVecKot10$Freq),] #order po številu živali

#funkcija za odbiro določenega števila živali izmed 1125 živali -v čredah z več kot 10 živalmi (brez polsester)
#začne z največjimi čredami
izberiCrede <- function (st) {
  sum <- 0
  crede <- c()
  row <- 1
  while (sum < st) {
    sum <- sum + credeVecKot10$Freq[row]
    crede <- c(crede, as.character(credeVecKot10$Var1[row]))
    row <- row +1
  }
  return(c(sum, crede))
}

izberiCrede(800)
vseRjaveICrede <- vseRjaveI[which(vseRjaveI$CRE_SIFRA_CREDA %in% credeVecKot10$Var1),]


#dodatna odbira - izloči index pod 0 --> ostanejo index NA in nad 0
Pod0 <- vseRjaveICrede$ZIV_ID_SEQ [(which(vseRjaveICrede$Index <=0))]
Nad0 <- vseRjaveICrede[which(!(vseRjaveICrede$ZIV_ID_SEQ %in% Pod0)),]
length(unique(Nad0$CRE_SIFRA_CREDA))
NoveCrede <- as.data.frame(table(Nad0$CRE_SIFRA_CREDA) )
NoveCrede10 <- NoveCrede[which(NoveCrede$Freq > 10),]
sum(NoveCrede10$Freq)
min(NoveCrede$Freq)

######################################
#od živali, ki imajo PV
#živali z isto čredo in istim očetom
two <- zivaliData[,c(1,7)]
two$together <- paste(two$CRE_SIFRA_CREDA, two$OCE, sep="")
two <- two[order(two$together),]
zivaliData$together <- paste(zivaliData$CRE_SIFRA_CREDA, zivaliData$OCE, sep="")
duplSis <- zivaliData$together[duplicated(zivaliData$together)]
#polsestre v čredni
halfSis <- zivaliData[which(zivaliData$together %in% duplSis),] #156 unique čreda-oče
halfSis <- halfSis[order(halfSis$together),]

#živali brez polsester v redči
uniqueZivali <- zivaliData[which(!(zivaliData$together %in% duplSis)),]

#chose the best halfSis
halfSis <- halfSis[order(halfSis$CRE_SIFRA_CREDA, halfSis$OCE, halfSis$Index),]
bestIndex <- halfSis[!(duplicated(halfSis$together, fromLast=T)),]

zivaliDataBrezHS <- rbind(uniqueZivali, bestPV)
sel <- qplot(zivaliDataBrezHS$Index, geom='histogram', bins=100) + xlab("932 krav")
all <- qplot(zivaliDataBrezHS$Index, geom="histogram", bins=100) + xlab("Vse rjave krave")
allPV <- qplot(zivaliDataBrezHS$VREDNOST_12_PV, geom="histogram", bins=100)+ xlim(c(80, 150))
allR <- qplot(zivaliDataBrezHS$R, geom="histogram", bins=100)

qplot(zivaliData$Index, geom="histogram", bins=100)
StPoCredah <- as.data.frame(table(zivaliDataBrezHS$CRE_SIFRA_CREDA))
StPoCredah <- StPoCredah[order(-StPoCredah$Freq),]
sum(StPoCredah$Freq[which(StPoCredah$Freq >= 10)])
zivaliBrezHDNad10 <- zivaliDataBrezHS[which(zivaliDataBrezHS$CRE_SIFRA_CREDA %in% StPoCredah$Var1[which(StPoCredah$Freq >= 10)]),]
sel1 <- qplot(zivaliBrezHDNad10$Index, geom='histogram', bins=100) + xlab("555 krav")
selPV1 <- qplot(zivaliBrezHDNad10$VREDNOST_12_PV, geom='histogram', bins=100) + xlim(c(80, 150))
selR1 <- qplot(zivaliBrezHDNad10$R, geom='histogram', bins=100) 
multiplot(sel1,sel, all)
multiplot(sel1, all)
multiplot(selPV1, allPV)
multiplot(selR1, allR)

###########################################################################################
#to so vse živali v čredah, ki po vseh kriterijih še vedno imajo več kot 10 krav (črede od 555 krav)
zivaliBrezHS <- podatki10[which(podatki10$CRE_SIFRA_CREDA %in% zivaliBrezHDNad10$CRE_SIFRA_CREDA),]
qplot(zivaliBrezHS$R, geom = 'histogram')
sum(StPoCredah$Freq[which(StPoCredah$Freq < 10)])

Nad1.5 <- zivaliDataBrezHS[which(zivaliDataBrezHS$Index > 0.5),]
STpoCredah <- as.data.frame(table(Nad1.5$CRE_SIFRA_CREDA))
STpoCredah <- STpoCredah[order(-STpoCredah$Freq),]


#two <- zivaliDataUnique[,c(1,7)] #MUST BE 0
#two$dupl <- duplicated(two) | duplicated(two, fromLast=T)

zivaliDataDupl <- zivaliData[two$dupl,]
zivaliDataDupl <- zivaliDataDupl[order(zivaliDataDupl$CRE_SIFRA_CREDA),]

Nad1.5 <- zivaliData[which(zivaliData$Index > 0.5),]
STpoCredah <- as.data.frame(table(Nad1.5$CRE_SIFRA_CREDA))
STpoCredah <- STpoCredah[order(-STpoCredah$Freq),]


reje10 <- avgSSIRj[which(avgSSIRj$RJKRAVE>=10),]
reje10Under120 <- reje10[which(reje10$CREDA %in% rejciUnder120$CRE_SIFRA_CREDA),]



################################################################################
############################################################################
 #zapiši tabelo za vsako čredi v pdf

library(gridExtra)
pdf("VseZivali.pdf", height=11, width = 8.5)
grid.table(vseRjave)
dev.off()



\documentclass{article}
\usepackage{longtable}
\begin{document}

<<results='asis'>>=
  library(xtable)

df = data.frame(matrix(rnorm(400), nrow=100))
xt = xtable(df)
print(xt, 
      tabular.environment = "longtable",
      floating = FALSE
)
@
  \end{document}