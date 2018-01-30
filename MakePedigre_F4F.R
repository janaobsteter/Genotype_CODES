ped <- read.csv("~/Documents/F4F/Rezultati_MCPKlasiak/RjavePedigree_09012018.csv", na=" ")
genoZ <- read.csv("~/Documents/F4F/GenotipiziraneF4F.txt", header=FALSE)
colnames(genoZ) <- "ID"
length(intersect(ped$ZIVAL, genoZ$ID))

#for which animals do you data - i.e. which are genotyped
dataG <- c(ped$ZIVAL %in% genoZ$ID)

#trimPed - tukaj izubereš, koliko generacij nazaj
library(pedigree)
Order <- orderPed(ped[,c(4,6,5)])
pedO <- ped[order(Order),]
pedT <- pedO[trimPed(pedO[,c(4,6,5)], data=dataG, ngenback = 5),]

pedT$ZIVAL <- as.character(pedT$ZIVAL)
pedT$MATIST <- as.character(pedT$MATIST)
pedT$OCEST <- as.character(pedT$OCEST)


pedT$MATIST[pedT$MATIST == ""] <- 0
pedT$ZIVAL[pedT$ZIVAL == ""] <- 0
pedT$OCEST[pedT$OCEST == ""] <- 0

write.table(pedT[,c(4,6,5)], "~/Documents/F4F/Rezultati_MCPKlasiak/PedigreeBlupf90.ped", na="0", quote=FALSE, row.names=FALSE, col.names=FALSE)



#Dat blupf90 file
Data <- read.csv("~/Documents/F4F/Rezultati_MCPKlasiak/TabelaRezultati_11122017.csv", sep="\t")
DataBlup1 <- Data[,c("CompleteID", "Rejec", "r", "k20", "a30", "Mascoba", "proteini", "Laktoza", "SSC", "CRE_SIFRA_CREDA", "STLAKT", "DIM", "KapaCSN")]
DataBlup1$KapaCSN <- gsub("/", "", DataBlup1$KapaCSN)
write.table(DataBlup1, "~/Documents/F4F/Rezultati_MCPKlasiak/Blupf90.dat",  quote=FALSE, row.names=FALSE, sep=" ", na="0", col.names=FALSE)
write.table(DataBlup1, "~/Documents/F4F/Rezultati_MCPKlasiak/Blupf90_NAMES.dat",  quote=FALSE, row.names=FALSE, sep=" ", na="0")
