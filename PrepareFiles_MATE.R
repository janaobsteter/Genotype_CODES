#GENDDER NI OK; KER GA DELA ALPHASIM!!!!!!! TI GA DOLOČIŠ POSEBEJJ!!!!!!
ped <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/PedigreeAndGeneticValues_cat.txt", header = TRUE)
indopt <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/LimitMale_PreSelectFemale/IndOpt.txt")


gInd <- ped[ped$Indiv %in% indopt$V1,c("Indiv", "sex")]
gInd$sex1 <- ifelse (gInd$sex =="M", 1, 2)
write.table(gInd[c("Indiv", "sex1")], "/home/jana/bin/AlphaMateLinux/OCSSloPop/LimitMale_PreSelectFemale/GENDER.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")


EBV <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/renumbered_Solutions_40")
write.table(EBV[EBV$V2 %in% indopt$V1, c(2,3)], "/home/jana/bin/AlphaMateLinux/OCSSloPop/LimitMale_PreSelectFemale/CRITERION.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")

#GENDDER NI OK; KER GA DELA ALPHASIM!!!!!!! TI GA DOLOČIŠ POSEBEJJ!!!!!!
ped <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/PedigreeAndGeneticValues_cat.txt", header = TRUE)
indopt <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//IndOpt.txt")


gInd <- ped[ped$Indiv %in% indopt$V1,c("Indiv", "sex")]
gInd$sex1 <- ifelse (gInd$sex =="M", 1, 2)
write.table(gInd[c("Indiv", "sex1")], "/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//GENDER.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")


EBV <- read.table("/home/jana/bin/AlphaMateLinux/OCSSloPop/renumbered_Solutions_40")
write.table(EBV[EBV$V2 %in% indopt$V1, c(2,3)], "/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//CRITERION.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")

