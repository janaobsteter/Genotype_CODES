gen0 <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/GenPed_EBV0.txt")
ped0 <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", header=TRUE, sep=" ")
sol <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/renumbered_Solutions", header=FALSE, sep=" ")
colnames(sol) <- c("renID", "Indiv", "EBV")


pedSol <- merge(ped0, sol, by="Indiv", all.y=TRUE)
cor(pedSol$EBV, pedSol$gvNormUnres1)

pedn <- pedSol[pedSol$cat == "potomciNP",]
cor(pedn$EBV, pedn$gvNormUnres1)

write.table(ped0[,c(2,3,4)], "/home/jana/Documents/PhD/CompBio/TestingGBLUP/Blupf90.ped", sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
