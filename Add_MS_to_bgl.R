######################################################
#add lab MS to the bgl file
#read in laboratory MS
LabMS <- read.csv("~/Genotipi/MS_Imputation/LabMS_14072016.csv", header=F) 
#read in animals that have SNPs ans MSs
LabInter <- read.csv ("~/Genotipi/MS_Imputation/MS_03062016/BeagleImputed1000it/Lab_and_Imp.txt", header=F)
#extract only the lab ms for the animals that also have SNPs
LabMS <- LabMS[(c(1,(which(LabMS$V1 %in% LabInter$V1)))),]
#read in a text file with MS by chromosome
MSChr <- read.table("~/Genotipi/MS_Imputation/MSByChr.txt", header=T)

#read in beagle file
#read in beagle format file of the genotypee
Chr1 <- read.table("ChromosomeSampleBeagleFile") # this is beagle format unphased file of the sample
Chr1phased <- read.table("ChromosomePhasedFile")
AniID <- as.data.frame(t(Chr1[1,])) #get Sample IDs - the ones that have both SNPs and Lab MS
colnames(AniID) <- "Col"
AniIDLab <- (which(AniID$Col %in% LabInter$V1)) #get a list of genotyped animals with lab MB
Chr1Lab <- Chr1[,c(1,2,AniIDLab)] #Extract corresponding animals with SNP and Lab MS

#get the min set SNP for the chromosome from the sample bgl file - extract required SNPs
Chr1Lab$V2 <- toupper(Chr1Lab$V2)
#length(intersect(Chr1Lab$V2, Chr1phased$V2)) # 71 - 70 SNPov + ID - MS, ki ga ni
Chr1Lab1 <- Chr1Lab[(which(Chr1Lab$V2 %in% Chr1phased$V2)),] #extracted SNPs and animals, SNPs from phased files extracted from sample file

#order animals in LabMS by the order in Chr1 bgl file
IDs <- as.data.frame(t(Chr1Lab1[1,]))
IDs <- as.data.frame(IDs[-c(1,2),])
IDs <- unique(IDs)
colnames(IDs) <- "ID"
colnames(LabMS)[1] <- "ID"
LabMS1 <- merge(IDs, LabMS, by="ID", all.y=T) #get the same order
LabMS2 <- as.data.frame(t(LabMS1))

#change the names of the MS for the two alleles to have the same name, without _1/2
MSnames <- as.vector(LabMS2[,1])
#colnames(MSnames) <- "Name"
soda <- seq(from=2, to=length(LabMS), by=2)
MSnames <- MSnames[soda]
MSnames <- rep(MSnames, each=2)
LabMS2$V1 <- c("ID",MSnames )


#prepare lab MS table in the same format as bgl file
#split the datasets for one to hold the first allele and the second to hold the other
liha <- seq(from=3, to=nrow(LabMS2), by=2)
soda <- seq(from=2, to=nrow(LabMS2), by=2)
LabMS_al1 <- LabMS2[c(1,soda),] #create two datasets
LabMS_al2 <- LabMS2[c(1,liha),]
TwoAl <- cbind(LabMS_al1, LabMS_al2)#cbind them back togetehr
TwoAl_ordered <- TwoAl[,order(TwoAl[1,])] # order by column name to get the same individuals next to eachtoher- i.e. the two alelles
TwoAl_ordered[,1] <- c("I", rep("M", 12)) #remove duplicated ID column

#add the corresponding MS line to the bgl file and order the marker as in marker and reference file
TwoAl_ordered <- TwoAl_ordered[-1,]
TwoAl_ordered <- TwoAl_ordered[c(1,2,6,3,4,7,11,10,5,9,8,12),]
TwoAl_ordered$V1.1[which(TwoAl_ordered$V1.1=="INRA23")] <- "INRA023" #change INRA23 to INRA023 to match the McClure tables
MSal <- TwoAl_ordered[(which(MSChr$Chr==ChrOfTheMS)),]
colnames(MSal) <- colnames(Chr1Lab1)
Chr1Ref <- rbind(Chr1Lab1, MSal)


#order the markers by the order in the marker file
marker <- read.csv("MarkerFile", sep="\t", header=F) #read marker file in
markers <- as.data.frame(marker$V1) #get the marker names
colnames(markers) <- "Marker"
id <- data.frame(Marker="ID")
markers <- rbind(id, markers)
colnames(Chr1Ref)[2] <- "Marker"
Chr1Ref <- merge(markers, Chr1Ref, by="Marker", sort=F) #merge to get the same order (markers first file and sort is false)
Chr1Ref <- Chr1Ref[,c(2,1,3:ncol(Chr1Ref))] #rearrange columns to get the Beagle format again


write.table(Chr1Ref, "FinalFilePathAndName", row.names=F, col.names=F, quote=F)




