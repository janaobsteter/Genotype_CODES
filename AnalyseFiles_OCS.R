setwd("/home/jana/bin/AlphaMateLinux/TestEddie/")

criterion <- read.table("CRITERION.txt", sep=",")
pedA<- read.table("PedigreeNrm.txt")
indopt <- read.table("IndOpt.txt")
ped <- read.table("PEDIGREE.txt", sep=",")

length(pedA[,1])

length(intersect(ped$V1, pedA[,1]))
length(intersect(criterion$V1, pedA[,1]))
length(intersect(indopt$V1, pedA[,1]))
length(intersect(criterion$V1, indopt$V1))


ped <- read.csv("~/ExternalPedigreeTotal.txt")
ind <- read.csv("~/IndOpt.txt", header=FALSE)

table(ped$cat[ped$Indiv %in% ind$V1])
table(ped$sex[ped$Indiv %in% ind$V1])

ped  <- read.table("~/PedigreeAndGeneticValues_cat.txt", header=TRUE, sep=" ")
ext <- read.csv("~/ExternalPedigreeTotal.txt")
ext1 <- read.csv("~/ExternalPedigree.txt", header=FALSE)

extSex <- ext[,c("Indiv", "sex")]
colnames(extSex) <- c("Indiv", "sexExt")
pedExt <- merge(extSex, ped, by="Indiv")
table(pedExt$sex, pedExt$sexExt)
