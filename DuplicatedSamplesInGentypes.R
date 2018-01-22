map <- read.csv("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/SNPchimp_result_3782403083.csv")
index <- read.csv("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/BovineSNP50_Final_SNPs_54001.csv", sep=" ")

nrow(index)
nrow(map)
length(intersect(map$SNP_name, index$SNP.Name))
colnames(index) [1] <- "SNP_name"
Map <- merge(index, map, by="SNP_name")
Map$positionStop <- Map$position

#write 705 map format - ampak ti ne koristi, ker je file 704
write.table(Map[,c("SNP_name", "SNP.Index", "chromosome", "position", "ITB_index")], "~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/50Kv01_705.map",
            row.names=FALSE, col.names=FALSE, quote=FALSE)
#write PLINK map
MapPlink <- Map[,c("Chr", "SNP_name", "position", "positionStop" )]
MapPlink <- MapPlink[order(MapPlink$Chr),]
write.table(MapPlink,
            "~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_50Kv01//PLINKnew.map", quote=FALSE, row.names=FALSE, col.names=FALSE)

ped <- read.csv("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/file704.20091221.svn", sep=" ", header=FALSE)
pedO <- ped[order(ped$V3),]
write.table(pedO, "~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/file704.20091221_Ordered.svn", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pedO[,c(1)], "~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/file704.20091221_Ordered1.svn", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pedO[,c(2,3,4,5)], "~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/file704.20091221_Ordered2.svn", col.names=FALSE, row.names=FALSE, quote=FALSE)

ped <- read.csv("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/SNP_pod_RJ_a/file704.20091221.svn", sep=" ", header=FALSE)
ped$V10 <- as.numeric(ped$V10)
pedO = ped[order(ped$V10),]


ped$V9 <- as.numeric(ped$V9)
ped$V12[10:99] <- ped$V11[10:99]
ped$V11[10:99] <- ped$V10[10:99]
ped$V10[10:99] <- ped$V9[10:99]

###########################################
#to je analizra novoprispelih vzorcev
sample <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_15012018/we_bl_12012018_Sample_Map.txt", sep="\t")
snp <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_15012018/we_bl_12012018_SNP_Map.txt", sep="\t")
ped <- read.csv("~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_15012018/we_bl_12012018/we_bl_final_report_12012018.txt", sep="\t", skip=9)
ped$Sample.ID <- as.character(ped$Sample.ID)
length(unique(ped$Sample.ID))
length(intersect(ped$Sample.ID, sample$ID))
nrow(ped) / nrow(snp) == nrow(sample)
length(intersect(snp$Name, ped$SNP.Name)) == nrow(snp)

id1 <- ped[ped$Sample.ID == sample$ID[1],]
id1SNPs <- ped$SNP.Name[ped$Sample.ID == sample$ID[1]]
for (idNum in 2:40) {
  id2SNPs <- ped$SNP.Name[ped$Sample.ID == sample$ID[idNum]]
  if (identical(id1SNPs, id2SNPs) == FALSE)  {
    print(identical(id1SNPs, id2SNPs))
    print(as.character(sample$ID[idNum]))
  }
}

SpInd <- ped [ped$Sample.ID=="SI53595706",]
ID1 <- SpInd[duplicated(SpInd$SNP.Name, fromLast=FALSE),]
ID2 <- SpInd[(duplicated(SpInd$SNP.Name, fromLast=TRUE)),]
sum(ID1$GC.Score)
sum(ID2$GC.Score)
#any deviation from 0 in log R ratio are evidence for copy number change
sum(ID1$Log.R.Ratio, na.rm=TRUE)
sum(ID2$Log.R.Ratio, na.rm=TRUE)
sum(ID1$Allele1...AB=="-")
sum(ID2$Allele1...AB=="-")
#we take ID1

identical(id1SNPs, ID1$SNP.Name) #ZDJ STA ISTA
#izbriši tole podvojeno vn iz zgonjege ped fajla
ped1 <- ped[ped$Sample.ID != "SI53595706",]
ped1 <- rbind(ped1, ID1)
write.table(ped1, "~/Genotipi/Genotipi_DATA/Rjava_TEMP/Genotipi_15012018/we_bl_Final_report_12012018.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
