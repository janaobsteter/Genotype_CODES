#seq has to have two columns - first is animals ID (Interbull or SIxxx številka) and the second the sequence

#################################################
sifrant_SNP <- read.csv("~/Genotipi/Genotipi_CODES/Sifrant_SNP.csv", sep=",")
#write.table(sifrant_SNP, "~/Genotipi/Genotipi_CODES/Sifrant_SNP.csv", sep=",", quote=F, row.names=F)
seq <- read.csv("~/Genotipi/Genotipi_DATA/Holstein/InterbullID_seq.csv", sep=",", header=T)
Sample800_geno <- read.table("~/Genotipi/Genotipi_DATA/Holstein/Genotipi_09082016/50Kv02_chip800.ped")
Sample800_SNPs <- read.table("~/Genotipi/Genotipi_DATA/Holstein/Genotipi_09082016/50Kv02_chip800.map")


rownames(Sample800_geno) <- Sample800_geno$V2
Sample800_geno <- Sample800_geno[,-c(1,2,3,4,5,6)]


AllelNames <- rep(Sample800_SNPs$V2,each=2)
colnames(Sample800_geno) <- AllelNames
Sample800_geno$ID <- rownames(Sample800_geno)

#############################################################
#pripraviza Govedo vsak SNP alel svoja vrstica
#############################################################
SNP800 <- matrix(ncol=3, nrow=1) #create an empty table
colnames(SNP800) <- c("ID", "Allele", "SNPCode")

#put each column underneath the other - get the table with all MS in one column
for (i in 1:(ncol(Sample800_geno)-1)) {
  SNP_col <- data.frame(Sample800_geno[,c(ncol(Sample800_geno),i)])
  colnames(SNP_col) <- c("ID", "Allele")
  SNP_col$SNPCode <- colnames(Sample800_geno[i])
  SNP800 <- rbind(SNP800, SNP_col)
}

SNP800 <- SNP800[-1,] #remove first row
SNP800 <- SNP800[order(SNP800$ID),]

#sort SNP800 file to add codes
SNP800_sorted <- SNP800[order(SNP800[,1],SNP800[,3],SNP800[,2]),]

#NUMBER OF animals*1600 rows in SNP800 table (800 SNPs * 2 alleles)


liha <- seq(from=1, to=((length(unique(SNP800_sorted$SNPCode)))*2*(nrow(Sample800_geno))), by=2)
soda <- seq(from=2, to=((length(unique(SNP800_sorted$SNPCode)))*2*(nrow(Sample800_geno))), by=2)
SNP800_sorted$SNPCode <- as.character(SNP800_sorted$SNPCode)

#get codes from sifrant when adding genotyped to already existing table
for (i in liha) { 
  SNP <- SNP800_sorted$SNPCode[i]
  SNPallele <- paste (SNP, "_1", sep="")
  SNP800_sorted$SNPCode[i] <- SNPallele
}

for (i in soda) { 
  SNP <- SNP800_sorted$SNPCode[i]
  SNPallele <- paste (SNP, "_2", sep="")
  SNP800_sorted$SNPCode[i] <- SNPallele
}

colnames(SNP800_sorted) <- c("ID", "ZGPL_VREDNOST", "SGPL_NAZIV")


SNP800_codes <- as.data.frame(matrix(nrow=1, ncol=4))
colnames(SNP800_codes) <- c("SGPL_NAZIV", "ID", "ZGPL_VREDNOST", "SGPL_SIFRA")
#add codes of SNPs from sifrant
for (i in unique(SNP800_sorted$ID)) {
  Ind <- subset(SNP800_sorted, SNP800_sorted$ID==i)
  IndCode <- merge(Ind, sifrant_SNP, by="SGPL_NAZIV", all=T)
  IndCode$ID <- IndCode$ID[1]
  SNP800_codes <- rbind(SNP800_codes, IndCode)
}




SNP800_codes <- SNP800_codes[-1,]


#rearrange and remove MS_name
#add sequence
colnames(seq) [1] <- "ID"
#seq <- seq[,c(2,3)]
SNP800_sorted1 <- merge (SNP800_codes, seq, by="ID")
SNP800_sorted2 <- SNP800_sorted1[,c(5, 4, 3)]
colnames(SNP800_sorted2) <- c("ZGP_ZIV_ID_SEQ","ZGP_SGPL_SIFRA","ZGP_VREDNOST")
SNP800_sorted2$ZGP_VREDNOST <- gsub("0", "", SNP800_sorted2$ZGP_VREDNOST)


write.table(SNP800_sorted2, "~/Genotipi/Genotipi_DATA/Holstein/Parental800_forGovedo.csv", sep=",", row.names=F, quote=F, na="")
