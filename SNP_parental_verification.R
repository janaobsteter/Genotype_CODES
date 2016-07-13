setwd("~/Genotipi/Genotipi_DATA/Genotipi1_12042016/PLINK_genotypeFiles/")
#1) READ IN MAP FILES AND FIND THE SNPS COMMON TO ALL CHIPS
GGP_map <- read.csv("GGP/OUTPUT/PLINK_MERGED.map", sep=" ", header=F)
GGP3_map <- read.csv("GGPv03/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
GP3_map <- read.csv("GP3v02/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
GP4_map <- read.csv("GP4/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
HD_map <- read.csv("HD/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
K50_map <- read.csv("50K/OUTPUT/PLINK_MERGED.map", sep="\t", header=F)
HD138K_map <- read.csv("/home/janao/Genotipi/Genotipi_DATA/Genotipi03062016/PLINK_genotypeFiles/HD138K/OUTPUT/PLINK_MERGED.map", sep=" ", header=F)
ISAG200 <- read.csv("~/Genotipi/ParentalVerification_SNPSNP/ISAG_200SNPs.csv", header=F, sep=",") #200 proposed ISAG SNPs
ICBF <- read.table("~/Genotipi/ParentalVerification_SNPSNP/ICBF_Parentage_SNP_Selection.csv", sep=",", header=T) #McClure proposed 800 SNPs
ICBF800 <- subset(ICBF, ICBF$X800_Parentage_SNP==1)
K50_map1 <- read.csv("~/Genotipi/Genotipi_DATA/SNPChimp_Maps//Illumina50K_SNPChimp.csv", sep=",") #read in SNPChimP 50K SNP names
min <- read.csv("~/Genotipi/MS_Imputation/MinSNP_Set.txt", header=F)


#Check concordance between chips and min SNPset
ISAGcore <- ISAG200[(which(ISAG200$V2=="ISAGcore")),]
ISAGadd <- ISAG200[(which(ISAG200$V2=="ISAGadditional")),]

length(intersect(GGP_map$V2, ISAG200$V1))
length(intersect(GGP3_map$V2, ISAG200$V1))
length(intersect(GP4_map$V2, ISAG200$V1))
length(intersect(HD_map$V2, ISAG200$V1))
length(intersect(HD138K_map$V2, ISAG200$V1))
length(intersect(K50_map$V2, ISAG200$V1))

length(intersect(GGP_map$V2, ISAGcore$V1))
length(intersect(GGP3_map$V2, ISAGcore$V1))
length(intersect(GP4_map$V2, ISAGcore$V1))
length(intersect(HD_map$V2, ISAGcore$V1))
length(intersect(HD138K_map$V2, ISAGcore$V1))
length(intersect(K50_map$V2, ISAGcore$V1))

length(intersect(GGP_map$V2, ISAGadd$V1))
length(intersect(GGP3_map$V2, ISAGadd$V1))
length(intersect(GP4_map$V2, ISAGadd$V1))
length(intersect(HD_map$V2, ISAGadd$V1))
length(intersect(HD138K_map$V2, ISAGadd$V1))
length(intersect(K50_map$V2, ISAGadd$V1))

#check match between ICBF 800 and out chips
length(intersect(GGP_map$V2, ICBF800$SNPID))
length(intersect(GGP3_map$V2, ICBF800$SNPID))
length(intersect(GP4_map$V2, ICBF800$SNPID))
length(intersect(HD_map$V2, ICBF800$SNPID))
length(intersect(HD138K_map$V2, ICBF800$SNPID))
length(intersect(K50_map$V2, ICBF800$SNPID))
length(intersect(K50_map1$SNP_name, ICBF800$SNPID))

#distribution of MAF amongst ISAG 100/200 in our Brown-Swiss population
MAFs <- read.table("/home/janao/Genotipi/Genotipi1_12042016/StepImputation/StepAllInOne/ImputationOnto50K/MAF.frq", header=T)
MAFs$rank <- c(1:nrow(MAFs))
MAF_ISAG100 <- MAFs[(which(MAFs$SNP %in% ISAGcore$V1)),]
hist(MAF_ISAG100$MAF)
MAF_ISAG200 <- MAFs[(which(MAFs$SNP %in% ISAGadd$V1)),]
hist(MAF_ISAG200$MAF)



MAFrank <- MAFs[order(-MAFs$MAF),] # rank by decreasing MAF
MAFtop100 <- MAFrank[1:100,] #top 100 with highest MAF
MAFtop200 <- MAFrank[1:200,] # top 200 with highest MAF
length(intersect(MAFtop100$SNP, ISAGcore$V1))
length(intersect(MAFtop100$SNP, ISAGcore$V1))
length(intersect(MAFtop100$SNP, ISAGadd$V1))
length(intersect(MAFtop200$SNP, ISAGadd$V1))

#bind together ICBF names, chip names and rs numbers
names <- merge(ICBF800, K50_map1, by="rs")
#there are 904 unique SNP names (last column is the name on our chip) --> extract them, write to a txt file and extract using PLINK
Chip800 <- as.data.frame(names$SNP)
write.table(Chip800, "~/Genotipi/Chip800.txt", row.names=F, col.names=F, quote=F)


colnames(names)[ncol(names)] <- "SNP"
colnames(MAFs)[5] <- "MAFsample"
MAF <- merge (MAFs, names, by="SNP") # 799 rows
hist(MAF$MAFsample, labels=T, main=NULL, ylim=c(0,200))


#Check Min SNPSet for MS imputation names with Illumina and GeneSeek names
min <- read.csv("~/Genotipi/MS_Imputation/Mcclure_Supplementary_Tables/Mcclure_MinSNPset.csv", sep=",", header=T)
K50_UMD <- read.csv("~/Genotipi/Illumina50Kv02_UMD3.csv", sep=",")
length(intersect(K50_UMD$position, min$UMD3Position)) #52
length(intersect(K50_UMD$SNP_name, min$MarkerName)) #56
GP3_UMD <- read.csv("~/Genotipi/Genotipi_DATA/SNPChimp_Maps/GP3_UMD3.csv", sep=",")
length(intersect(GP3_UMD$position, min$UMD3Position)) #52
length(intersect(GP3_UMD$SNP_name, min$MarkerName)) #56
HD <- read.csv ("~/Genotipi/Genotipi_DATA/SNPChimp_Maps/Illumina700K.csv", sep=",")
length(intersect(HD$SNP_name, min$MarkerName)) #880
rs880 <- which(HD$SNP_name %in% min$MarkerName)
rs880 <- HD[c(rs880),]
colnames(min)[1] <- "SNP_name"
rs880 <- merge (rs880, min, by="SNP_name")
write.table(rs880, "~/Genotipi/MS_Imputation/MinSNPSet_rsNum.csv", row.names = F, quote = F)


#check intersect of min880 SNPs with rs numbers
length(intersect(GP3_UMD$rs, rs880$rs)) # 682
length(intersect(GGP$rs, rs880$rs)) # 607
length(intersect(HD$rs, rs880$rs)) # 751
length(intersect(HDv02$rs, rs880$rs)) # 840
length(intersect(K50_map1$rs, rs880$rs)) # 56
K50v01 <- read.csv ("~/Genotipi/Genotipi_DATA/SNPChimp_Maps/Illumina50Kv01.csv", sep=",")
length(intersect(K50v01$rs, rs880$rs)) # 57

#check 800 for parental verification with GeneSeek Chips using rs
length(intersect(ICBF800$rs, GP3_UMD$rs)) #800
GGP <- read.csv("~/Genotipi/Genotipi_DATA/SNPChimp_Maps/GGPv02_map.csv", sep=",")
length(intersect(ICBF800$rs, GGP$rs))
HD <- read.csv("~/Genotipi/Genotipi_DATA/SNPChimp_Maps/HD_map.csv", sep=",")
length(intersect(ICBF800$rs, HD$rs))
HDv02 <- read.csv("~/Genotipi/Genotipi_DATA/SNPChimp_Maps/HDv02_map.csv", sep=",")
length(intersect(ICBF800$rs, HDv02$rs))
length(intersect(ICBF800$rs, K50_map1$rs))

#get the SNP names on the chips for these 800 SNPs
GGP_800 <- (intersect(ICBF800$rs, GGP$rs))
GGP_names <- GGP[(which(GGP$rs %in% GGP_800)),]
HD_800 <- (intersect(ICBF800$rs, HD$rs))
HD_names <- HD[(which(HD$rs %in% HD_800)),]
HD2_800 <- (intersect(ICBF800$rs, HDv02$rs))
HD2_names <- HDv02[(which(HDv02$rs %in% HD2_800)),]

Names800 <- as.data.frame(GGP_names$SNP_name)
write.table(Names800, "~/Genotipi/ParentalVerification_SNPSNP/Names_800SNPs.txt", row.names=F, col.names = F, quote=F)

############################################################3
#read in merged files with 800 SNPs from samples
###########################################################
Sample800_geno <- read.table("/home/janao/Genotipi/ParentalVerification_SNPSNP/Chip800_Brown.ped")
Sample800_SNPs <- read.table("/home/janao/Genotipi/ParentalVerification_SNPSNP/Chip800_Brown.map")
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
SNP800_sorted <- SNP800[order(SNP800[,1],SNP800[,3],SNP800[,2]),]



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
#add codes of SNPs from sifratn
for (i in unique(SNP800_sorted$ID)) {
  Ind <- subset(SNP800_sorted, SNP800_sorted$ID==i)
  IndCode <- merge(Ind, sifrant_SNP, by="SGPL_NAZIV", all=T)
  IndCode$ID <- IndCode$ID[1]
  SNP800_codes <- rbind(SNP800_codes, IndCode)
}

SNP800_codes <- SNP800_codes[-1,]



#APPLICABLE WHEN CREATING SIFRANT
#add additional column with codes, MS within each animal sorted in the same way (order)
#Names are AllleNames, codes are 1:1600
SNP800_sorted$Code <- NA
SNP800_sorted$Code <- rep(1:((length(unique(SNP800_sorted$SNPCode)))*2), times=(nrow(Sample800_geno))) #codes for SNP 
# for (i in Ani_SNP) { #assign codes in order within an animal
#   for (l in 0:1599) {
#     SNP800_sorted$Code[i+l] <- Codes[l+1] 
#   }
# }

#rearrange and remove MS_name
#add sequence
seq <- read.csv("~/Genotipi/Rjave_seq_ID.csv", sep=",")
colnames(seq) [1] <- "ID"
SNP800_sorted1 <- merge (SNP800_sorted, seq, by="ID")
SNP800_sorted2 <- SNP800_sorted1[,c(5, 4, 2)]
colnames(SNP800_sorted2) <- c("ZGP_ZIV_ID_SEQ","ZGP_SGPL_SIFRA","ZGP_VREDNOST")
SNP800_sorted2$ZGP_VREDNOST <- gsub("0", "", SNP800_sorted2$ZGP_VREDNOST)

#create sifrant for SNPs - WHEN CREATIN SIFRANT
# sifrant_SNP <- SNP800_sorted1[1:1600,c(3,4)]
# liha <- seq(from=1, to=1600, by=2)
# soda <- seq(from=2, to=1600, by=2)
# sifrant_SNP$SNPCode <- as.character(sifrant_SNP$SNPCode)
# for (i in liha) { 
#   SNP <- sifrant_SNP$SNPCode[i]
#   SNPallele <- paste (SNP, "_1", sep="")
#   sifrant_SNP$SNPCode[i] <- SNPallele
# }
# 
# for (i in soda) { 
#   SNP <- sifrant_SNP$SNPCode[i]
#   SNPallele <- paste (SNP, "_2", sep="")
#   sifrant_SNP$SNPCode[i] <- SNPallele
# }


# colnames(sifrant_SNP) <- c("SGPL_NAZIV", "SGPL_SIFRA")

write.table(SNP800_sorted2, "~/Genotipi/ParentalVerification_SNPSNP/Govedo_SNP800.csv", sep=",", row.names=F, quote=F)
#SNP_sorted <- read.table("~/Genotipi/ParentalVerification_SNPSNP/Govedo_SNP800.csv", sep=" ", header=T)
# write.table(sifrant_SNP, "~/Genotipi/ParentalVerification_SNPSNP/Sifrant_SNP.csv", sep=",", quote=F, row.names=F)
