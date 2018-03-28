ped <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/PedigreeAndGeneticValues_cat.txt", sep= " ", header=TRUE)

#pedG <- read.table("~/Gen10_Pedigree.txt", sep=" ", header=TRUE)

table(ped$age[ped$cat=="k"])
table(ped$age[ped$cat=="telF"])
table(ped$age[ped$cat=="pt"])

potomciM <- ped$Indiv[(ped$cat == "potomciNP") & (ped$sex == "M")]
krave5 <- ped$Indiv[(ped$cat=="k") & (ped$age %in% 1:5)]
ostaleZ <- ped$Indiv[ped$cat %in% c("pBM", "pt", "telF")]


write.table(ped[,c(2,3,4)], "/home/jana/bin/AlphaMateLinux/OCSSloPop/PEDIGREE.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep=" ")
write.table(sort(c(potomciM, krave5, ostaleZ)), "/home/jana/bin/AlphaMateLinux/OCSSloPop/IndOpt.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  

gender <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/Gender.txt")
write.table(gender[,c(2,3)],"/home/jana/bin/AlphaMateLinux/OCSSloPop/Gender.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


#output files
#non-zero contributors
parents <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/ContributionsModeOptTarget1.txt", header=TRUE)
