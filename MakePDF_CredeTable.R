#script for producing PDFs for each table
#F4F projekt - odbria živali za genotipizacijo
#vsaka čreda posebna tabela
#primarni seznam - izbrane živali in seznam B - seznam nadomestnih živali v primeru, da je ena iz primarnega seznama izločena
setwd('/home/jana/Documents/F4F/OdbiraZivali/PripravaPDF')
library(Hmisc)
library(knitr)

##  make my data
vseRjavePDF <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017_PDF.csv",header=T) 
vseRjavePDF$DAT_ROJSTVO <- as.Date(vseRjavePDF$DAT_ROJSTVO, format='%d.%m.%Y')

seznamB <- read.csv("~/Documents/F4F/OdbiraZivali/seznamB_22022017.csv", header=T)
seznamB$DAT_ROJSTVO <- as.Date(seznamB$DAT_ROJSTVO, format='%d.%m.%Y')

rejci <- read.csv('~/Documents/F4F/OdbiraZivali/PripravaPDF/OdbraneCrede.csv')
colnames(rejci)[1] <- "Rejec"



##izloči reje, ki niso primerne
izlCrede <- c(10253,3989,30968,86,10135,2475,31629,9690,15727,2264,31656,3110,11610) #te reje je izločil Zoran Kramer
vseRjavePDF <- vseRjavePDF[-which(vseRjavePDF$CRE_SIFRA_CREDA %in% izlCrede),]
seznamB <- seznamB[-which(seznamB$CRE_SIFRA_CREDA %in% izlCrede),]

#izberi živali v čredah (glede na število živali za genotipizacijo)
stPocredah <- as.data.frame(table(vseRjavePDF$CRE_SIFRA_CREDA)) #tabela število krav po čredah
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
  print(sum)
  return(crede)
}



izbraneCrede <- izberiCrede(300)
vseRjavePDF <- vseRjavePDF[which(vseRjavePDF$CRE_SIFRA_CREDA %in% izbraneCrede),]
seznamB <- seznamB[which(seznamB$CRE_SIFRA_CREDA %in% izbraneCrede),]


## knitr loop
for (creda in as.character(32717)){ #unique(vseRjavePDF$CRE_SIFRA_CREDA)
  knit2pdf("MakePDFTable.Rnw", output=paste0('PDF_', creda, '.tex'))
}

