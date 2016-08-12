##########################
#add phased to reference file
Ref1 <- read.table("ReferenceFile", stringsAsFactors = F)
Sample1 <- read.table("NewPhasedFile", stringsAsFactors = F)

#read in marker file to order both reference and sample file according to the marker file order
marker <- read.csv("MarkerFile", sep="\t", header=F) #read marker file in
markers <- as.data.frame(marker$V1) #get the marker names
colnames(markers) <- "Marker"
id <- data.frame(Marker="id")
markers <- rbind(id, markers)
colnames(Ref1) [2] <- "Marker"
RefMarker <- Ref1[(which(Ref1$Marker %in% markers$Marker)),] #extract minSNPset markers from Reference file
RefMarker <- merge(markers, RefMarker, by="Marker", sort=F, all.x=T) #merge with marker file and leave in all the marker SNPs
RefMarker <- RefMarker[(match(markers$Marker, RefMarker$Marker)),] # order reference file 

#merge sample with marker file and order according to the marker file
Sample1$V2[1] <- ("id")
colnames(Sample1)[2] <- "Marker"
SampleMarker <- merge(markers, Sample1, by="Marker", sort=F, all.x=T)
SampleMarker <- SampleMarker[(match(markers$Marker, SampleMarker$Marker)),] # order reference file 

#now you can merge pahsed sample with MS and reference together to produce a new reference
RefMarker <- RefMarker[,-2] #remove I; M column to avoid duplication in the merged file
SampleRef <- merge(SampleMarker, RefMarker, by="Marker", sort=F)
#rearrange columns to get beagle format
SampleRef <- SampleRef[,c(2,1,3:ncol(SampleRef))]
#for lines where sample had no values replace NA with M in the first column
SampleRef$V1[which(is.na(SampleRef$V1))] <- "M"

write.table(SampleRef, "AddedReference", quote=F, row.names=F, col.names = F, na="?")

