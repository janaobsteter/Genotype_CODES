#first argument is name of the ped file
#second argument is the working dir
args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
WorkingDir <- args[2]
######
#dobi SNPe za lastnosti
#######
library(ggplot2)
#setwd("/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/")
setwd(WorkingDir)
tsPed <- read.table(paste0(name, "_Monogenic.ped"))
tsMap <- read.table(paste0(name, "_Monogenic.map"))
ts <- tsPed[,-c(3,4,5,6)]

conv <- read.csv("~/Genotipi/Genotipi_CODES/Trait_conversion_map_IDBV3_3.csv")

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
#write.csv(namesSNP, "/home/jana/Genotipi/Genotipi_CODES/Monogenic_IDBv03.txt", quote=FALSE)
namesSNP$TraitNo <- 1:nrow(namesSNP)
traitNo <- rep(1:nrow(namesSNP), each=2)
new <- data.frame(t(data.frame(TraitNo = as.factor(c(0,0, traitNo)))))
colnames(new) <- colnames(ts)
ts <- rbind(new, ts) #tukaj oštevilči SNPe

Tts <- as.data.frame(t(ts))
#za vsako lastnost naredi tabelo

snpTable <- data.frame(FID=tsPed$V1, ID=tsPed$V2)
#10 = Weaver, 11 = Arachnomelia, 13 = ABCG2, 32 = SMA 
traitNames <- c()
namesSNP$Full.Trait.Name <- as.character(namesSNP$Full.Trait.Name)
for (tNo in c(10, 11, 13, 14, 15, 31)) {
  snpTable_new <- as.data.frame(t(Tts[Tts$TraitNo %in% c(0,tNo),]))
  colnames(snpTable_new) <- c("FID", "ID", "A1", "A2")
  snpTable_new <- snpTable_new[rownames(snpTable_new) != "TraitNo",]
  snpTable_new$Genotype <- paste0(snpTable_new$A1, snpTable_new$A2)
  colnames(snpTable_new)[3:5] <- c(paste0("A1_", tNo), paste0("A2_", tNo), paste0("Genotype_",tNo) )
  snpTable <- merge(snpTable, snpTable_new, by=c("FID", "ID"))
  alleles <- c(snpTable$A1, snpTable$A2)
  traitNames <- c(traitNames, paste0(namesSNP$Full.Trait.Name[namesSNP$TraitNo == tNo], "g") )
    #ggplot(as.data.frame(alleles)) +  geom_bar(aes(x = alleles, fill = as.factor(alleles)), position = "dodge", stat = "count") + ggtitle(namesSNP$Full.Trait.Name[namesSNP$TraitNo==tNo])
}


###################################
#preberi še kapa in beta casein
###################################
kapaCSN <- read.csv(paste0(name, '_KappaCaseinGenotype_python.csv'))
kapaCSN <- kapaCSN[,c("ID","SKUPEN")]
colnames(kapaCSN) <- c("ID", "KapaCSN")
#če hočeš dobiti genotipi s slashom za izpise
kapaCSN$KapaCSN <- gsub("BB", "B/B", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("AA", "A/A", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("AB", "A/B", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("AC", "A/C", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("BC", "B/C", kapaCSN$KapaCSN)
kapaCSN$KapaCSN <- gsub("BE", "B/E", kapaCSN$KapaCSN)
#tukaj še za beta kazein
#betaCSN <- read.csv(name, '_BetaCaseinGenotype_python.csv', header=TRUE)[,c(2,5)]
#colnames(betaCSN) <- c("ID", "BetaCSN")
######################################
#sutvari tabele po čredah
#############################

#write.csv(snpi, '/home/jana/Documents/F4F/TabelaMonogenskeBolezniSNPi.csv', quote=FALSE, row.names=FALSE)
#snpi <- read.csv('/home/jana/Documents/F4F/TabelaMonogenskeBolezniSNPi.csv')

GenTable <- snpTable[,c(2,seq(5, ncol(snpTable), by=3))]
traitNames <- gsub("Bovine Progressive Degenerative Myeloencephalopathy", "Weaver", traitNames)
traitNames <- gsub(" (Weaver)", "", traitNames)
traitNames <- gsub("Arachnomelia Syndrome_SUOX", "Arahnomelija", traitNames)
traitNames <- gsub("Spinal Muscular Atrophy", "SMA", traitNames)
traitNames <- c("ID", traitNames)
colnames(GenTable) <- traitNames#c("ID", "Weaverg", "Arahnomelijag", "ABCG2g", "KappaCasein Bg", "KappaCasein Eg", "SMAg")
GenTable$Weaver <- ifelse( GenTable$Weaverg == "BB", "Zdrava", ifelse( GenTable$Weaverg == "AB"| GenTable$Weaverg=="BA",  "Prenašalka", NA) )
GenTable$Arahnomelija <- ifelse( GenTable$Arahnomelijag == "AA", "Zdrava", ifelse( GenTable$Arahnomelijag == "AB"| GenTable$Arahnomelijag =="BA",  "Prenašalka", NA) )
GenTable$ABCG2 <- ifelse( GenTable$ABCG2g == "AA", "Dve kopiji alela A", ifelse( GenTable$ABCG2g == "AB"| GenTable$ABCG2g == "BA",  "Ena kopija alela A", NA) )
GenTable$SMA <- ifelse( GenTable$SMAg == "BB", "Zdrava", ifelse( GenTable$SMAg == "AB"| GenTable$SMAg == "BA",  "Prenašalka", NA) )
GenTable <- merge(GenTable, kapaCSN, by="ID", all.x=TRUE)
#GenTable <- merge(GenTable, betaCSN, by="ID", all.x=TRUE)
#GenTable <- GenTable[,c("1,8,9,11,10,12,13"ID)]

#zapiši to tabelo
write.csv(GenTable, paste0(WorkingDir, "MonogenicGenotypes_Table.csv"), quote=FALSE, row.names=FALSE)
print("Successfully created MonogenicGenotypes_Table.csv")